# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

Changelog:
    20171116: fix scipy version 1.0
    20170310: add override params for pixmm
    20161219: Removed class variables
    20160815: Restructuring (clean up)
    20160803: Fix for RestrictROINormi13 where edge is both min and max
    20160802: sync with wad2.0
    20160701: Fixes for cropping, artefacts; fix ints
    20160205: Fix integer overflow on 32bit
    20160202: Finished uniformity
    20160201: Split Uniformity/Artefact detection off from QCMammo to enable recycling; starting from v20150522
"""
__version__ = '20171116'
__author__ = 'aschilham'

try:
    import pydicom as dicom
except ImportError:
    import dicom
import numpy as np
from math import pi
import copy
import scipy.ndimage as scind
import matplotlib.pyplot as plt
from PIL import Image # image from pillow is needed
from PIL import ImageDraw # imagedraw from pillow is needed, not pil
import scipy.misc
# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
if scipy_version[0] == 0:
    if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
        raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")

# First try if we are running wad1.0, since in wad2 libs are installed systemwide
try: 
    # try local folder
    import wadwrapper_lib
except ImportError:
    # try pyWADlib from plugin.py.zip
    try: 
        from pyWADLib import wadwrapper_lib

    except ImportError: 
        # wad1.0 solutions failed, try wad2.0
        try: 
            # try system package wad_qc
            from wad_qc.modulelibs import wadwrapper_lib
        except ImportError: 
            # use parent wad_qc folder, and add it to search path
            import sys
            # add root folder of WAD_QC to search path for modules
            _modpath = os.path.dirname(os.path.abspath(__file__))
            while(not os.path.basename(_modpath) == 'Modules'):
                _new_modpath = os.path.dirname(_modpath)
                if _new_modpath == _modpath:
                    raise
                _modpath = _new_modpath
            sys.path.append(os.path.dirname(_modpath))
            from wad_qc.modulelibs import wadwrapper_lib


class UnifStruct:
    def _max_pixel_value(self):
        # determine max allowed pixelvalue; as 2^bits_stored -1; note this value is not properly stored in dcmfileIn.BitsStored!
        dicomfields = [ ["0028,0101",  "Bits Stored"]]
        key = dicomfields[0][0]
        dicvalue = self.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
        if dicvalue == "":
            raise ValueError('Invalid dcmfile for max_pixel_value')

        return (2**dicvalue)-1

    def __init__ (self, dcmInfile, pixeldataIn):
        """
        As this is going to be used from different modules
        Just make sure the image is already inverted if needed!
        """
        self.verbose = False

        # input image
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn

        self.max_pixel_value = self._max_pixel_value() # max allowed value to store per pixel

        # for matlib plotting
        self.hasmadeplots = False

        # uniformity
        self.means = []
        self.stdevs = []
        self.unif_pct = -1
        self.snr_hol = -1
        self.unif_rois = [] # xy roi definitions # format: x0,wid, y0,hei

        # artefacts
        self.art_clusters = []
        self.art_image = None
        self.art_borderpx = []
        self.art_threshold = -1
        self.art_rois = [] # x0,y0,rad
        self.art_crop = [] # [xmin,xmax, ymin,ymax]

        # working with a smaller FOV
        self.unif_crop_frac = 1
        self.unif_crop_inoutoverin = -1
        self.unif_crop = []

        self.pixmm = None # allow hard override of pixmm, for example is ImagerPixelSpacing does not exist
        
    def pixDim(self, axis=0):
        if self.dcmInfile is None:
            return None
        if not self.pixmm is None:
            return self.pixmm
        try:
            pix_mm = self.dcmInfile.PixelSpacing[axis]
        except:
            pix_mm = self.dcmInfile.ImagerPixelSpacing[axis]
        return pix_mm

class Uniformity_QC:
    def __init__(self):
        self.qcversion = __version__
        self.artefactDetectorParameters() # default for non-L50 mammo

    def artefactDetectorParameters(self, UseStructure=False,
                                   bkscale=25, # dscale
                                   fgscale=5.0, # iscale
                                   threshold=15,
                                   ):
        self.artdet_UseStructure = UseStructure # use structure detector or LocalSNR
        # for Hologic Selenia there is a part of the image (12 px) that should be ignored
        self.artdet_bkscale = bkscale # if given, then subtract Gaussian(bkscale) from filtered image
        self.artdet_fgscale = fgscale # Gaussian scale of features
        self.artdet_threshold = threshold # threshold for artefact or not

    def RestrictROIUniformity(self, cs, borderpx=[0,0,0,0]):
        # attempt to crop to a box without over-exposed edges
        # define 3 horizontal rois: x= 20-40%, 40-60%, 60-80%; y=0-20%
        # calculate for each roi mean, sd; position (mean-n.sd)<=val<=(mean+n.sd); cut-off is max pos
        # check from center towards edge because edge might be min surrounded by max.
        roipts = []

        small_offset_px = int(3./cs.pixDim(0))# 0 # set > 0 to make found region smaller with this offset from all edges
        nsdev = .75
        nparts = 5
        widthpx  = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        xmin0_px = borderpx[0]
        xmax0_px = widthpx-borderpx[1]-1
        ymin0_px = borderpx[2]# cs.border_offset_px
        ymax0_px = heightpx-borderpx[3] -1 #cs.border_offset_px -1

        stepx = int(widthpx/nparts)
        stepy = int(heightpx/nparts)
        posx = range(stepx, stepx*nparts, stepx)
        posy = range(stepy, stepy*nparts, stepy)
        sx2 = int(stepx/2)
        sy2 = int(stepy/2)
        # y first
        ymin_px = ymin0_px
        ymax_px = ymax0_px
        for x0 in posx[0:-1]:
            yoffset0 = [0,0]
            # min y pos
            smallIm = cs.pixeldataIn[x0:(x0+stepx), 0:sy2 ]#stepy]
            while np.min(smallIm) == np.max(smallIm): # fix for XA images with a lot of empty content
                yoffset0[0] += sy2
                smallIm = cs.pixeldataIn[x0:(x0+stepx), yoffset0[0]:yoffset0[0]+sy2 ]#stepy]

            mean  = np.mean(smallIm)
            stdev = np.std(smallIm)
            thresh = nsdev*stdev
            hei = np.shape(smallIm)[1]
            for iy in reversed(range(hei)):
                val = np.mean(smallIm[:,iy])
                if iy == 0 or (val<mean-thresh or val > mean+thresh):
                    ymin_px = max(ymin_px, yoffset0[0]+iy+1+small_offset_px)
                    break
            # max y pos
            smallIm = cs.pixeldataIn[x0:(x0+stepx), (posy[-1]+sy2):(posy[-1]+stepy)]#posy[-1]:(posy[-1]+stepy)]
            while np.min(smallIm) == np.max(smallIm): # fix for XA images with a lot of empty content
                yoffset0[1] += sy2
                smallIm = cs.pixeldataIn[x0:(x0+stepx), (posy[-1]-yoffset0[1]+sy2):(posy[-1]+stepy-yoffset0[1]) ]#stepy]

            mean  = np.mean(smallIm)
            stdev = np.std(smallIm)
            thresh = nsdev*stdev
            hei = np.shape(smallIm)[1]
            for iy in range(hei):
                val = np.mean(smallIm[:,iy])
                if iy == hei-1 or (val<mean-thresh or val > mean+thresh):
                    ymax_px = min(ymax_px, (posy[-1]+sy2-yoffset0[1])+iy-1-small_offset_px)
                    break

        print('[restrict_y]',mean,stdev, ymin0_px,ymax0_px, thresh,ymin_px,ymax_px)

        xmin_px = xmin0_px
        xmax_px = xmax0_px
        for y0 in posy[0:-1]:
            xoffset0 = [0,0]
            # min x pos
            smallIm = cs.pixeldataIn[0:sx2, y0:(y0+stepy)]
            while np.min(smallIm) == np.max(smallIm): # fix for XA images with a lot of empty content
                xoffset0[0] += sy2
                smallIm = cs.pixeldataIn[xoffset0[0]:xoffset0[0]+sx2, y0:(y0+stepy)]

            mean  = np.mean(smallIm)
            stdev = np.std(smallIm)
            thresh = nsdev*stdev
            wid = np.shape(smallIm)[0]
            for ix in reversed(range(wid)):
                val = np.mean(smallIm[ix, :])
                if ix == 0 or (val<mean-thresh or val > mean+thresh):
                    xmin_px = max(xmin_px, xoffset0[0]+ix+1+small_offset_px)
                    break

            # max x pos
            smallIm = cs.pixeldataIn[(posx[-1]+sx2):(posx[-1]+stepx), y0:(y0+stepy)]
            while np.min(smallIm) == np.max(smallIm): # fix for XA images with a lot of empty content
                xoffset0[1] += sy2
                smallIm = cs.pixeldataIn[(posx[-1]-xoffset0[1]+sx2):(posx[-1]+stepx-xoffset0[1]), y0:(y0+stepy)]

            mean  = np.mean(smallIm)
            stdev = np.std(smallIm)
            thresh = nsdev*stdev
            wid = np.shape(smallIm)[0]
            for ix in range(wid):
                val = np.mean(smallIm[ix, :])
                if ix == wid-1 or (val<mean-thresh or val > mean+thresh):
                    xmax_px = min(xmax_px, posx[-1]+sx2-xoffset0[1]+ix-1-small_offset_px)
                    break

        print('[restrict_x]',mean,stdev, xmin0_px, xmax0_px, thresh,xmin_px,xmax_px)

        roipts = [
            [xmin_px,ymin_px],
            [xmax_px,ymin_px],
            [xmax_px,ymax_px],
            [xmin_px,ymax_px]]
        cs.unif_crop = [xmin_px,xmax_px, ymin_px,ymax_px]
        area_in = 1.*(xmax_px-xmin_px)*(ymax_px-ymin_px)
        area_all = np.prod(np.shape(cs.pixeldataIn))#1.*(xmax0_px-xmin0_px)/(ymax0_px-ymin0_px)
        cs.unif_crop_frac = area_in/area_all

        area_out = area_all-area_in
        mean_in = np.mean(cs.pixeldataIn[xmin_px:xmax_px:3, ymin_px:ymax_px:3])
        mean_all = np.mean(cs.pixeldataIn[xmin0_px:xmax0_px:3,ymin0_px:ymax0_px:3])
        mean_out = (area_all*mean_all-area_in*mean_in)/area_out
        cs.unif_crop_inoutoverin = (mean_in-mean_out)/mean_in
        return roipts

    def NeedsCroppingUnif(self, cs):
        """
        If needscropping2 says no, check if the max gradient is too large
        """
        
        scale = cs.phantommm2pix(2.) # 2 mm
        threshold = cs.max_pixel_value/scale # from min to max within 2 mm; could be artefact... so restrict to 1.5 cm from each edge
        
        gap = int(cs.phantommm2pix(15.)) # 15 mm
                  
        need = self.NeedsCropping2(cs.unif, mode='uniformity')
        if need == False:
            pdCopy = cs.pixeldataIn.astype(float)
            Lx = (scind.gaussian_filter(cs.pixeldataIn, scale/2, order=[1,0], mode='mirror'))
            lx_max = 0
            part = Lx[0:gap, :]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[-gap:-1, :]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[:, 0:gap]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[:, -gap:-1]
            lx_max = max([lx_max, part.max(), -part.min()])
            if lx_max > threshold:
                return True
            Lx = (scind.gaussian_filter(cs.pixeldataIn, scale/2, order=[0,1], mode='mirror'))
            part = Lx[0:gap, :]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[-gap:-1, :]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[:, 0:gap]
            lx_max = max([lx_max, part.max(), -part.min()])
            part = Lx[:, -gap:-1]
            lx_max = max([lx_max, part.max(), -part.min()])
            if lx_max > threshold:
                return True
            print('####y', lx_max, threshold)
            
        return need


    def Uniformity(self,cs,borderpx=[0,0,0,0], border_is_circle=False):
        """
        Calculation of non-uniformity

        Defined in EUREF_DMWG_Protocol_2003 as max deviation between avg pix value in five ROIs
        """
        error = True

        width_array = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        height_array = np.shape(cs.pixeldataIn)[1]

        pixel_spacing = cs.pixDim(0)## spacing in mm

        field_max_x_px = width_array-borderpx[1] ## total field width/height in px
        field_max_y_px = height_array-borderpx[3]
        field_min_x_px = borderpx[0]
        field_min_y_px = borderpx[2]

        if len(cs.unif_crop)>0:
            field_min_x_px = max([cs.unif_crop[0], field_min_x_px])
            field_max_x_px = min([cs.unif_crop[1], field_max_x_px])
            field_min_y_px = max([cs.unif_crop[2], field_min_y_px])
            field_max_y_px = min([cs.unif_crop[3], field_max_y_px])

        ## According to EUREF_DMWG_Protocol_2003 the areas need to be 4cm^2

        ## ROI definitions:

        ## 0    1
        ##
        ##   4  5
        ##
        ## 2    3

        ## or in circle
        ##   0
        ## 1 4 2,5
        ##   3

        width_roi_mm = 20.0 #width of roi in mm
        width_roi_px = int(width_roi_mm/pixel_spacing)  #width of roi in px
        height_roi_mm = 20.0 #height of roi in mm
        height_roi_px = int(height_roi_mm/pixel_spacing)  #height of roi in px

        cs.unif_rois = [] # format: x0,wid, yo,hei
        x0 = field_min_x_px
        y0 = field_min_y_px
        x1 = field_max_x_px-width_roi_px
        y1 = field_max_y_px-height_roi_px

        if border_is_circle:
            x01 = int((x0+x1)/2)
            y01 = int((y0+y1)/2)
            cs.unif_rois.append([x01,width_roi_px, y0,height_roi_px]) # 0
            cs.unif_rois.append([x0,width_roi_px, y01,height_roi_px]) # 1
            cs.unif_rois.append([x1,width_roi_px, y01,height_roi_px]) # 2
            cs.unif_rois.append([x01,width_roi_px, y1,height_roi_px]) # 3
        else:
            cs.unif_rois.append([x0,width_roi_px, y0,height_roi_px]) # 0
            cs.unif_rois.append([x1,width_roi_px, y0,height_roi_px]) # 1
            cs.unif_rois.append([x0,width_roi_px, y1,height_roi_px]) # 2
            cs.unif_rois.append([x1,width_roi_px, y1,height_roi_px]) # 3
        x0 = field_min_x_px+(field_max_x_px-field_min_x_px)/2
        y0 = field_min_y_px+(field_max_y_px-field_min_y_px)/2

        cs.unif_rois.append([int(x0-width_roi_px/2),width_roi_px, int(y0-height_roi_px/2),height_roi_px]) # 4 
        cs.unif_rois.append([int(x1),width_roi_px, int(y0-height_roi_px/2),height_roi_px]) # 5

        cs.means  = []
        cs.stdevs = []
        for r in cs.unif_rois:
            arr = cs.pixeldataIn[r[0]:(r[0]+r[1]),r[2]:(r[2]+r[3])]
            cs.means.append(np.mean(arr))
            cs.stdevs.append(np.std(arr))

        maxdev = abs(cs.means[0]-cs.means[4])
        for i in range(1,4):
            maxdev = max(maxdev,abs(cs.means[i]-cs.means[4]))

        cs.unif_pct = 100.*maxdev/cs.means[4]
        cs.snr_hol = cs.means[5]/cs.stdevs[5]
        if cs.verbose:
            print("[Uniformity] maxdev="+str(maxdev)+" unif="+str(cs.unif_pct)+" snr="+str(cs.snr_hol))
        error = False
        return error

#----------------------------------------------------------------------
    def Artefacts(self, cs, borderpx=[0,0,0,0], border_is_circle=False, uiobject=None):
        """
        """
        error = True
        doUseStructureAlways = False
        print('[Artefacts]',borderpx)
        # copy image, and cutoff black borders
        field_max_x_px = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        field_max_y_px = np.shape(cs.pixeldataIn)[1]
        field_min_x_px = 0
        field_min_y_px = 0

        if len(cs.unif_crop)>0:
            field_min_x_px = cs.unif_crop[0]
            field_max_x_px = cs.unif_crop[1]
            field_min_y_px = cs.unif_crop[2]
            field_max_y_px = cs.unif_crop[3]

        field_max_x_px -= borderpx[1] 
        field_max_y_px -= borderpx[3]
        field_min_x_px += borderpx[0]
        field_min_y_px += borderpx[2]

        cs.art_crop = [field_min_x_px,field_max_x_px, field_min_y_px,field_max_y_px]

        pdCopy = np.array(cs.pixeldataIn[field_min_x_px:field_max_x_px,field_min_y_px:field_max_y_px], dtype=np.float)
        # exclude outside max circle if so defined
        if border_is_circle:
            mw, mh = np.shape(pdCopy) # width/height in pixels
            cxy = [ mw/2., mh/2.]
            rxy = [ mw/2., mh/2.]
            Y,X = np.ogrid[:mh, :mw]
            dist = np.sqrt( ((X-cxy[0])/rxy[0])**2. + ((Y-cxy[1])/rxy[1])**2.)
            
            circle_mask = np.transpose(dist<1.)
            avg_in = np.average(pdCopy[circle_mask])
            pdCopy[~circle_mask] = avg_in
            
            # set mask to slightly smaller to avoid border problems
            rxy = [ mw/2.-self.artdet_bkscale, mh/2.-self.artdet_bkscale]
            dist = np.sqrt( ((X-cxy[0])/rxy[0])**2. + ((Y-cxy[1])/rxy[1])**2.)
            circle_mask = np.transpose(dist<1.)

        if self.artdet_UseStructure: # use structure detector
            # Apart from spots, we also see lines appearing as large artefacts.
            # To detect the later, the slower "structure" is needed
            have_big_mem = True
            if have_big_mem:
                cs.art_image = StructureDetector(pdCopy,bksigma=self.artdet_bkscale,uiobject=uiobject) # remove background trend
            else:
                cs.art_image = SplitMergeStructureDetector(pdCopy,bksigma=self.artdet_bkscale,uiobject=uiobject) # remove background trend

        else:
            # We see only spots, so we do not need the slower "structure"
            cs.art_image = LocalSNR(pdCopy,sigma=self.artdet_fgscale,bksigma=self.artdet_bkscale,uiobject=uiobject)# remove background trend

        """
        find clusters of pixels BELOW a certain threshold
        """
        try:
            cs.ArtefactExcludeItems() # For L50, first exclude contrast plugs:
        except AttributeError:
            pass


        if uiobject:
            uiobject.pbar.startProgress(2,"checking histogram...")
        cwid = np.shape(cs.art_image)[0]
        chei = np.shape(cs.art_image)[1]
        clusters = []
        npix = 0
        # check histogram
        cs.art_clusters = []
        cs.art_borderpx = borderpx
        cwid4 = int(cwid/4)
        chei4 = int(cwid/4)
        smallIm = cs.art_image[cwid4:3*cwid4, chei4:3*chei4]
        smMean = np.mean(smallIm)
        smStDev = np.std(smallIm)
        #print('art',smMean, smStDev, (self.artdet_threshold-smMean)/smStDev)
        if self.artdet_UseStructure:
            self.artdet_threshold = smMean+10*smStDev
            #print('art', self.artdet_threshold)
            if not np.abs(cs.art_image).max()>self.artdet_threshold: # > above for Structure!
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs.art_image[::3, ::3], bins=[min(cs.art_image.min(),-self.artdet_threshold),-self.artdet_threshold,self.artdet_threshold,max(cs.art_image.max(),self.artdet_threshold)],density=False)
            frac = 1.*(hist[0]+hist[-1])/( (cwid/3)*(chei/3) ) # above for Structure!
        else:
            self.artdet_threshold = smMean+5*smStDev
            #print('art', self.artdet_threshold)
            if not np.abs(cs.art_image).max()>self.artdet_threshold: #<
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs.art_image[::3, ::3], bins=[min(cs.art_image.min(),-self.artdet_threshold),-self.artdet_threshold,self.artdet_threshold,max(cs.art_image.max(),self.artdet_threshold)],density=False)
            frac = 1.*(hist[0]+hist[-1])/( (cwid/3)*(chei/3) )

        cs.art_threshold = self.artdet_threshold 

        if frac > .1:
            if uiobject:
                uiobject.pbar.endProgress()
                uiobject.pbar.label.setText("Corrupt image!")
            print("[SNRArtefacts] ERROR! Corrupt image!",hist,self.artdet_threshold)
            if uiobject:
                uiobject.pbar.label.setText("Finished.")
            return error

        if uiobject:
            uiobject.pbar.doProgress("clustering artefacts...")
        cca = wadwrapper_lib.connectedComponents()

        # extra mapout
        if self.artdet_UseStructure:
            if border_is_circle: # exclude outside max circle if so defined
                cs.art_image[~circle_mask] = self.artdet_threshold-1e-3
            else:
                if borderpx[0]>0:
                    cs.art_image[:(borderpx[0]+1),:] = self.artdet_threshold-1e-3
                if borderpx[1]>0:
                    cs.art_image[-borderpx[1]:,:]    = self.artdet_threshold-1e-3
                if borderpx[2]>0:
                    cs.art_image[:,:(borderpx[2]+1)] = self.artdet_threshold-1e-3
                if borderpx[3]>0:
                    cs.art_image[:,-borderpx[3]:]    = self.artdet_threshold-1e-3
            cca_in = np.abs(cs.art_image)>self.artdet_threshold
        else:
            if border_is_circle: # exclude outside max circle if so defined
                cs.art_image[~circle_mask] = -self.artdet_threshold+1e-3
            else:
                if borderpx[0]>0:
                    cs.art_image[:(borderpx[0]+1),:] = -self.artdet_threshold+1e-3
                if borderpx[1]>0:
                    cs.art_image[-borderpx[1]:,:]    = -self.artdet_threshold+1e-3
                if borderpx[2]>0:
                    cs.art_image[:,:(borderpx[2]+1)] = -self.artdet_threshold+1e-3
                if borderpx[3]>0:
                    cs.art_image[:,-borderpx[3]:]    = -self.artdet_threshold+1e-3
            cca_in = np.abs(cs.art_image)>self.artdet_threshold

        cc_art_image,nclusters = cca.run(cca_in) 
        for a in range(1,nclusters+1):
            clusters.append(cca.indicesOfCluster(a))

        nc = len(clusters)
        if cs.verbose:
            print("[Artefacts] found",npix ," artefact pixels in",nc,"clusters")

        if uiobject:
            uiobject.pbar.endProgress()
        if cs.verbose:
            print("[Artefacts] Found",len(clusters),"clusters")

        cs.art_rois = []
        for c in clusters:
            npclus = np.array(c)
            xvals = npclus[:,0]
            yvals = npclus[:,1]
            minx = np.min(xvals)
            maxx = np.max(xvals)
            miny = np.min(yvals)
            maxy = np.max(yvals)
            rad = max(maxy-miny,maxx-minx,1.)/2.
            cs.art_rois.append([(maxx+minx+1)/2.,(maxy+miny+1)/2., rad])

            if cs.verbose:
                print("[Artefacts]...size",len(c))

        if cs is not None:
            cs.art_clusters = copy.deepcopy(clusters)
            cs.art_borderpx = borderpx
            cs.art_threshold = self.artdet_threshold

        if uiobject:
            uiobject.pbar.label.setText("Finished.")
        error = False
        return error

#----------------------------------------------------------------------
    def drawThickCircle(self,draw,x,y,rad,color,thick):
        for t in range(-int((thick-1)/2),int((thick+1)/2)):
            r1 = rad+t
            draw.ellipse((x-r1,y-r1,x+r1,y+r1), outline=color)

    def saveAnnotatedArtefactImage(self,cs,fname):
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        im = scipy.misc.toimage(cs.art_image.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

        # now draw all rois in reserved color
        if len(cs.art_rois) > 0:
            draw = ImageDraw.Draw(im)
            for xyr in cs.art_rois:
                self.drawThickCircle(draw, xyr[0],xyr[1], rad=xyr[2],color=0,thick=5)
            del draw

        # convert to RGB for JPG, cause JPG doesn't do PALETTE and PNG is much larger
        im = im.convert("RGB")

        imsi = im.size
        if max(imsi)>2048:
            ratio = 2048./max(imsi)
            im = im.resize( (int(imsi[0]*ratio+.5), int(imsi[1]*ratio+.5)),Image.ANTIALIAS)
        im.save(fname)

    #----- NEW-----
    def RestrictROINormi13(self, cs, borderpx=[0,0,0,0]):
        """
        Check if we need to look for high or low or both
        """
        error = False

        small_offset_px = 0# int(3./cs.pixDim(0))# 0 # set > 0 to make found region smaller with this offset from all edges

        nparts = 5
        widthpx  = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        xmin0_px = borderpx[0]
        xmax0_px = widthpx-borderpx[1]-1
        ymin0_px = borderpx[2]# cs.border_offset_px
        ymax0_px = heightpx-borderpx[3] -1 #cs.border_offset_px -1

        stepx = int(widthpx/nparts)
        stepy = int(heightpx/nparts)
        sx2 = int(stepx/2)
        sy2 = int(stepy/2)
        x0 = int(widthpx/2)
        y0 = int(heightpx/2)

        # compare to center
        smallIm = cs.pixeldataIn[(x0-sx2):(x0+sx2),(y0-sy2):(y0+sy2)]
        mean  = np.mean(smallIm)

        # find edge_val
        num = 5
        edge_max = np.max([np.max(cs.pixeldataIn[:num,:]),
                           np.max(cs.pixeldataIn[-num:,:]),
                           np.max(cs.pixeldataIn[:,:num]),
                           np.max(cs.pixeldataIn[:,-num:])])
        edge_min = np.min([np.min(cs.pixeldataIn[:num,:]),
                           np.min(cs.pixeldataIn[-num,:]),
                           np.min(cs.pixeldataIn[:,:num]),
                           np.min(cs.pixeldataIn[:,-num:])])

        thresh_hi = thresh_lo = None
        if edge_max > mean: # edge is higher val than object
            thresh_hi = edge_max - .2*(edge_max-mean)

        if edge_min < mean: # edge is lower val than object
            thresh_lo = edge_min + .2*(mean-edge_min)

        if thresh_hi is None: thresh_hi=cs.max_pixel_value
        if thresh_lo is None: thresh_lo=0

        edge_max1 = edge_max
        edge_min1 = edge_min
        while (thresh_lo>0 and edge_min1 == edge_min) or (thresh_hi<cs.max_pixel_value and edge_max1 == edge_max):
            roipts = self._RestrictROINormi13(cs, range_x=[xmin0_px,xmax0_px], range_y=[ymin0_px,ymax0_px], thresh_hi=thresh_hi, thresh_lo=thresh_lo, small_offset_px=small_offset_px)
            xmin1_px, xmax1_px, ymin1_px, ymax1_px =  roipts[0][0], roipts[2][0], roipts[0][1], roipts[2][1]
            cropim = cs.pixeldataIn[xmin1_px:xmax1_px, ymin1_px:ymax1_px]
            edge_max1 = np.max([np.max(cropim[:num,:]),
                                np.max(cropim[-num:,:]),
                                np.max(cropim[:,:num]),
                                np.max(cropim[:,-num:])])
            edge_min1 = np.min([np.min(cropim[:num,:]),
                                np.min(cropim[-num,:]),
                                np.min(cropim[:,:num]),
                                np.min(cropim[:,-num:])])
            xmin0_px = xmin1_px
            xmax0_px = xmax1_px
            ymin0_px = ymin1_px
            ymax0_px = ymax1_px

        cs.unif_crop = [roipts[0][0], roipts[2][0], roipts[0][1], roipts[2][1]]
        # here check for valid values
        return error

    def _RestrictROINormi13(self, cs, range_x, range_y, thresh_hi=None, thresh_lo=None, small_offset_px=0):
        # attempt to crop to a box without over-exposed edges
        # define 3 horizontal rois: x= 20-40%, 40-60%, 60-80%; y=0-20%
        # calculate for each roi mean, sd; position (mean-n.sd)<=val<=(mean+n.sd); cut-off is max pos
        # check from center towards edge because edge might be min surrounded by max.

        roipts = []
        nparts = 5
        xmin0_px, xmax0_px = range_x
        ymin0_px, ymax0_px = range_y

        widthpx  = xmax0_px-xmin0_px
        heightpx = ymax0_px-ymin0_px

        stepx = int(widthpx/nparts)
        stepy = int(heightpx/nparts)
        sx2 = 5#int(stepx/2)
        sy2 = 5#int(stepy/2)
        x0 = int((xmin0_px+xmax0_px)/2)
        y0 = int((ymin0_px+ymax0_px)/2)
        valmax = cs.max_pixel_value

        # y first
        ymin_px = ymin0_px
        ymax_px = ymax0_px

        # min y pos
        for iy in range(ymin0_px, ymax0_px):
            # only include values not max and not min
            vals = []
            for ix in range((x0-sx2),(x0+sx2)):
                v = cs.pixeldataIn[ix, iy]
                if v >1 and v<valmax:
                    vals.append(v)
            if len(vals) == 0:
                continue
            val = np.mean(vals)
            if val>thresh_lo and val<thresh_hi:
                ymin_px = max(ymin_px, iy+1+small_offset_px)
                break

        # max y pos
        for iy in reversed(range(ymin0_px, ymax0_px)):
            vals = []
            for ix in range((x0-sx2),(x0+sx2)):
                v = cs.pixeldataIn[ix, iy]
                if v >1 and v<valmax:
                    vals.append(v)
            if len(vals) == 0:
                continue
            val = np.mean(vals)
            if val>thresh_lo and val<thresh_hi:
                ymax_px = min(ymax_px, iy-1-small_offset_px)
                break
        print('[restrictnormi13_y]', ymin0_px,ymax0_px, thresh_lo,thresh_hi, ymin_px,ymax_px)

        xmin_px = xmin0_px
        xmax_px = xmax0_px
        # min x pos
        for ix in range(xmin0_px, xmax0_px):
            vals = []
            for iy in range((y0-sy2),(y0+sy2)):
                v = cs.pixeldataIn[ix, iy]
                if v >1 and v<valmax:
                    vals.append(v)
            if len(vals) == 0:
                continue
            val = np.mean(vals)
            if val>thresh_lo and val<thresh_hi:
                xmin_px = max(xmin_px, ix+1+small_offset_px)
                break

        # max x pos
        for ix in reversed(range(xmin0_px, xmax0_px)):
            vals = []
            for iy in range((y0-sy2),(y0+sy2)):
                v = cs.pixeldataIn[ix, iy]
                if v >1 and v<valmax:
                    vals.append(v)
            if len(vals) == 0:
                continue
            val = np.mean(vals)
            if val>thresh_lo and val<thresh_hi:
                xmax_px = min(xmax_px, ix-1-small_offset_px)
                break
        print('[restrictnormi13_x]', xmin0_px, xmax0_px, thresh_lo,thresh_hi, xmin_px,xmax_px)

        roipts = [
            [xmin_px,ymin_px],
            [xmax_px,ymin_px],
            [xmax_px,ymax_px],
            [xmin_px,ymax_px]]
        area_in = 1.*(xmax_px-xmin_px)*(ymax_px-ymin_px)
        area_all = np.prod(np.shape(cs.pixeldataIn))#1.*(xmax0_px-xmin0_px)/(ymax0_px-ymin0_px)
        cs.unif_crop_frac = area_in/area_all

        area_out = area_all-area_in
        mean_in = np.mean(cs.pixeldataIn[xmin_px:xmax_px:3, ymin_px:ymax_px:3])
        mean_all = np.mean(cs.pixeldataIn[xmin0_px:xmax0_px:3,ymin0_px:ymax0_px:3])
        mean_out = (area_all*mean_all-area_in*mean_in)/area_out
        cs.unif_crop_inoutoverin = (mean_in-mean_out)/mean_in
        return roipts

    def NeedsCropping2(self, cs, mode='normi13'):
        """
        Simple test: If along any of the edges the value 2^15 -1 is reached or, a value higher than 1.4* center than yes
        Also check the fraction of off pixels and the absolute value of the pixel
        """

        wid  = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        hei = np.shape(cs.pixeldataIn)[1]

        pixel_spacing = cs.pixDim(0)## spacing in mm

        wid_roi_mm = 30.0 #width of roi in mm
        wid_roi = int(wid_roi_mm/pixel_spacing)  #width of roi in px
        hei_roi_mm = 30.0 #height of roi in mm
        hei_roi = int(hei_roi_mm/pixel_spacing)  #height of roi in px

        x0 = int((wid-wid_roi)/2.)
        y0 = int((hei-hei_roi)/2.)

        # first the central value
        cmean = np.mean(cs.pixeldataIn[x0-int(wid_roi/2):x0+int(wid_roi/2),y0-int(hei_roi/2):y0+int(hei_roi/2)])

        hist,bins = np.histogram(cs.pixeldataIn[::3,::3], bins=500,density=False)
        thresh = self.otsu(hist,bins) # distinguish between background and foreground
        total = np.sum(hist)
        hist[bins[:-1]>=thresh] = 0
        low = np.sum(hist)
        frc = 1.*low/total
        if frc>0.97: # if almost all is background, likely background and foreground are swapped!
            frc = 1.-0.97

        num = 5
        edge_max = np.max([np.max(cs.pixeldataIn[:num,:]),
                           np.max(cs.pixeldataIn[-num:,:]),
                           np.max(cs.pixeldataIn[:,:num]),
                           np.max(cs.pixeldataIn[:,-num:])])
        edge_min = np.min([np.min(cs.pixeldataIn[:num,:]),
                           np.min(cs.pixeldataIn[-num,:]),
                           np.min(cs.pixeldataIn[:,:num]),
                           np.min(cs.pixeldataIn[:,-num:])])

        if mode == 'normi13':
            minfraction = 0.25
            if edge_min < max(0.05*cs.max_pixel_value+1.,0.6*cmean): 
                if frc>minfraction:
                    print('[NeedsCropping2] Needs -cropping (%f).'%(frc), edge_min,edge_max)
                    return True 

        elif mode == 'uniformity':
            minfraction = .025
            if edge_min < max(0.05*cs.max_pixel_value+1.,0.6*cmean) or edge_min <2 or edge_max>cs.max_pixel_value-1: 
                if frc>minfraction or edge_min <2 or edge_max>cs.max_pixel_value-1:
                    print('[NeedsCropping2] Needs -cropping (%f).'%(frc), edge_min,edge_max)
                    return True 

        print('[NeedsCropping2] Does not need cropping (%f).'%(1.*low/total))
        return False

    def otsu(self, hist, bins):
        """
        Otsu thresholding to distinguish between foreground and background
        """
        currentMax = 0
        threshold = 0
        sumTotal, sumForeground, sumBackground = 0., 0., 0.
        weightBackground, weightForeground = 0, 0
        totalPixels = np.sum(hist)

        # Calculate the total of the data
        for b,h in zip(bins,hist):
            sumTotal += b * h

        for b,h in zip(bins,hist):
            # Calculate the weight of the background
            weightBackground += h
            if( weightBackground == 0 ):
                continue

            # Calculate the weight of the foreground
            weightForeground = totalPixels - weightBackground
            if ( weightForeground == 0 ):
                break

            sumBackground += b*h

            # Calculate the mean of the classes
            meanB = 1.*sumBackground / weightBackground
            meanF = 1.*(sumTotal - sumBackground) / weightForeground

            # Calculate variance between classes
            varBetween = 1.*weightBackground*weightForeground
            varBetween *= (meanB-meanF)*(meanB-meanF)

            # Check if the variance between classes is greater than the current best
            if(varBetween > currentMax):
                currentMax = varBetween
                threshold = b

        return threshold

    #----------------------------------------------------------------------
def SplitMergeStructureDetector(inImage, bksigma = None,uiobject=None):
    """
    L50 images cause a huge memory hog.
    To save some memory, the idea is to split the input image into parts,
    process each part, and reassemble before return.
    Workflow:
    1. start with an empty result LIST
    2. define a part of the inImage, taking into account overlapping for blurring
    3. make a copy of that subpart of type float (inImage is explicitely expected to be INT)
    4. process that sub part without bkremoval
    5, append the processed image (with exception of overlap) to the result
    6. repeat 2-5 for all other parts
    7. convert result to np.array, do bkRemoval and return
    """
    rows = np.shape(inImage)[0] ## width/height in pixels
    cols = np.shape(inImage)[1]

    #print("[SplitMergeStructureDetector] inImage first/last/rows",0,rows-1,rows)

    # 1. start with an empty result LIST
    result = np.empty_like(inImage, dtype=float)
    overlap = 4*(5+2)+1 #4*(iscale+dscale)+1; filter is 4 sigma wide
    xparts = [ # first last keepfirst keeplast
               [0,             int(rows/2) +overlap,  0,     int(rows/2)-1],
               [int(rows/2)-overlap,rows,             int(rows/2),rows-1]
               ]

    for xpart in xparts:
        pdCopy = (inImage[xpart[0]:xpart[1]:]).astype(float)
        #print("[SplitMergeStructureDetector] pdCopy first/last/rows",xpart[0],xpart[1]-1,np.shape(pdCopy)[0])
        pdProc = StructureDetector(pdCopy, None,uiobject)

        firstkeep = xpart[2]-xpart[0]
        lastkeep  = xpart[3]-xpart[0]
        offset    = xpart[2]-firstkeep
        #print("[SplitMergeStructureDetector] Restore first/last",firstkeep+offset,lastkeep+offset)
        for i,row in enumerate(pdProc):
            if i<firstkeep or i>lastkeep:
                continue
            result[i+offset] = row

    if bksigma is not None:
        blurRes = scind.gaussian_filter(result,bksigma,order=[0,0])
        result -= blurRes
        if uiobject:
            print("[SplitMergeStructureDetector] 5/4 trend removed")

    return result

def StructureDetector(inImage, bksigma = None,uiobject=None):
    """
    Code is helaas slecht leesbaar geworden, maar gebruikte teveel geheugen.
    """
    #import resource
    #print("[Structure] resMB",1,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) # max usage so far)

    """
    Return numpy array with largest eigenvalues of Structure Tensor
    (integrated over some area)

    Structure tensor: S = [Lx^2 LxLy]
                          [LyLx Ly^2]
    integrate = blur each element over gauss
    Warning: X and Y direction changed in pyqtgraph
    WATCH IT: should be float, else overflow errors later (not mentioned!)
    """
    dscale = 2.0
    iscale = 5.0

    #Step 1: make sure inImage is of type float
    try:
        uiobject.pbar.startProgress(7,"...converting to float")
    except:
        pass
    pdCopy = inImage
    if type(inImage[0,0]) != type(1.0):
        pdCopy = inImage.astype(float)

    #Step 2: Elements of Structure tensor
    try:
        uiobject.pbar.doProgress("...calculating Lx**2")
    except:
        #print("Cannot do progress")
        pass

    Lx = (scind.gaussian_filter(pdCopy,dscale,order=[1,0],mode='mirror'))
    Lxx = Lx**2
    try:
        uiobject.pbar.doProgress("...calculating Ly**2")
    except:
        pass
    Ly = (scind.gaussian_filter(pdCopy,dscale,order=[0,1],mode='mirror'))
    Lyy = Ly**2
    try:
        uiobject.pbar.doProgress("...calculating Lx*Ly")
    except:
        pass

    # optimize copies in mem: 4 in use
    Lx *= Ly # actually LxLy Lxy = Lx*Ly
    Ly = None # mark as free: 3 in use

    #Step 3 Integrate over iscale
    try:
        uiobject.pbar.doProgress("...smoothing Lx**2")
    except:
        pass

    Lxx = scind.gaussian_filter(Lxx,iscale,order=[0,0],mode='mirror')
    try:
        uiobject.pbar.doProgress("...smoothing Ly**2")
    except:
        pass

    Lyy = scind.gaussian_filter(Lyy,iscale,order=[0,0],mode='mirror')
    try:
        uiobject.pbar.doProgress("...smoothing Lx*Ly")
    except:
        pass

    Lx = scind.gaussian_filter(Lx,iscale,order=[0,0],mode='mirror') # Lxy = scind.gaussian_filter(Lxy,iscale,order=[0,0]); 3 in use

    # Step 2: calculate the eigenvalues of the Hessian matrix:
    #  Lpp = 0.5*(Lyy+Lxx- sqrt((Lyy+Lxx)**2-4*(LyyLxx-Lyx**2)))
    #  Lqq = 0.5*(Lyy+Lxx+ sqrt((Lyy+Lxx)**2-4*(LyyLxx-Lyx**2)))
    try:
        uiobject.pbar.doProgress("...calculating ev")
    except:
        pass

    Lxx += Lyy                # Ma= Lyy+Lxx ; 3 in use
    Lyy = (Lxx-Lyy)*Lyy # Lyy*Lxx-Lxy**2: 3 in use
    Lyy -= Lx**2 # Lyy*Lxx-Lxy**2: 3 in use
    Lx = None # Mark as free: 2 in use
    Lyy = np.sqrt(Lxx**2-4*Lyy) # Mb= np.sqrt(Ma**2-4*(Lyy*Lxx-Lxy**2)) # 2 in use

    Lxx -= Lyy        # Lpp = 0.5*(Ma-Mb): 2 in use
    Lyy = Lxx + 2*Lyy # Lqq = 0.5*(Ma+Mb): 2 in use

    Lxx = np.maximum(Lxx,Lyy)/2. # Lout = np.maximum(Lpp,Lqq) : 2 in use
    Lyy = None # Mark as free:1 in use

    if bksigma is not None:
        Lyy = scind.gaussian_filter(Lxx,bksigma,order=[0,0],mode='mirror')
        Lxx -= Lyy
        if uiobject:
            print("[StructureDetector] 5/4 trend removed")

    try:
        uiobject.pbar.endProgress()
    except:
        pass

    return Lxx

def LocalSNR(pSrc, sigma,bksigma = None, uiobject=None):
    """
    Local SNR: [ mean(x,y) ]/stdev(x,y)
    Can be approximated as [ Gauss{I,sigma}] / sqrt[ Gauss{I-Gauss{I,sigma},sigma}^2]
    """
    if uiobject:
        uiobject.pbar.startProgress(3,"Calculating LocSNR")
    blurIm = scind.gaussian_filter(pSrc,sigma,order=[0,0])
    if uiobject:
        print("[LocalSNR] 1/4 blur'd")
    devIm = pSrc-blurIm
    if uiobject:
        uiobject.pbar.doProgress("dev'd")
        print("[LocalSNR] 2/4 dev'd")
    sdIm  = np.sqrt(scind.gaussian_filter(devIm**2,sigma,order=[0,0]))
    sdIm[sdIm<1.e-6]=1. # prevent div by zero
    if uiobject:
        uiobject.pbar.doProgress("std'd")
        print("[LocalSNR] 3/4 std'd")
    locnormIm = blurIm/sdIm
    if uiobject:
        print("[LocalSNR] 4/4 snr'd")
    if bksigma is not None:
        blurIm = scind.gaussian_filter(locnormIm,bksigma,order=[0,0])
        locnormIm -= blurIm
        if uiobject:
            print("[LocalSNR] 5/4 trend removed")
    if uiobject:
        uiobject.pbar.endProgress()
    #blurIm = scind.gaussian_filter(locnormIm,sigma,order=[2,0])**2+scind.gaussian_filter(locnormIm,sigma,order=[0,2])**2
    #return blurIm
    return locnormIm

def LocalNorm(pSrc, sigma,bksigma = None, uiobject=None):
    """
    Local Norm: [ I - mean(x,y) ]/stdev(x,y)
    Can be approximated as [ I-Gauss{I,sigma}] / sqrt[ Gauss{I-Gauss{I,sigma},sigma}^2]
    """
    if uiobject:
        uiobject.pbar.startProgress(4,"Calculating LocalNorm")
    blurIm = scind.gaussian_filter(pSrc,sigma,order=[0,0])
    if uiobject:
        print("[LocalNorm] 1/4 blur'd")
    devIm = pSrc-blurIm
    if uiobject:
        uiobject.pbar.doProgress("dev'd")
        print("[LocalNorm] 2/4 dev'd")
    sdIm  = np.sqrt(scind.gaussian_filter(devIm**2,sigma,order=[0,0]))
    sdIm[sdIm<1.e-6]=1. # prevent div by zero
    if uiobject:
        uiobject.pbar.doProgress("std'd")
        print("[LocalNorm] 3/4 std'd")
    locnormIm = devIm/sdIm
    if uiobject:
        print("[LocalNorm] 4/4 normed'd")
    if uiobject:
        uiobject.pbar.endProgress()

    return locnormIm

