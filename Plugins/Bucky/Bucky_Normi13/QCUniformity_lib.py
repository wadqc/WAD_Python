# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 08:50:51 2014

@author: aschilha
TODO:

Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

Changelog:
    20160205: Fix integer overflow on 32bit
    20160202: Finished uniformity
    20160201: Split Uniformity/Artefact detection off from QCMammo to enable recycling; starting from v20150522
"""
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
if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
    raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")
try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib

__qcversion__=20160205

class UnifStruct:
    verbose = False

    # input image
    dcmInfile = None
    pixeldataIn = None
    mustbeinverted = False
    
    # for matlib plotting
    hasmadeplots = False

    # uniformity
    means = []
    stdevs = []
    unif_pct = -1
    snr_hol = -1
    unif_rois = [] # xy roi definitions # format: x0,wid, y0,hei
    
    # artefacts
    art_clusters = []
    art_image = None
    art_borderpx = []
    art_threshold = -1
    art_rois = [] # x0,y0,rad

    # working with a smaller FOV
    expertmode = False
    expert_roipts = []
    expert_frac = -1
    expert_inoutoverin = -1

    def __init__ (self,dcmInfile,pixeldataIn):
        self.verbose = False
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.hasmadeplots = False
        self.means = []
        self.stdevs = []
        self.unif_pct = -1
        self.snr_hol = -1
        self.unif_rois = []
        self.art_clusters = []
        self.art_image = None
        self.art_borderpx = []
        self.art_threshold = -1
        self.art_rois = []
        self.expertmode = False
        self.expert_roipts = []
        self.expert_frac = -1
        self.expert_inoutoverin = -1

    def pixDim(self,axis=0):
        if self.dcmInfile is None:
            return None
        try:
            pix_mm = self.dcmInfile.PixelSpacing[axis]
        except:
            pix_mm = self.dcmInfile.ImagerPixelSpacing[axis]
        return pix_mm

class Uniformity_QC:
    qcversion = __qcversion__

    def __init__(self):
        self.artefactDetectorParameters() # default for non-L50 mammo

    def artefactDetectorParameters(self,UseStructure=False,
                                   bkscale=25, # dscale
                                   fgscale=5.0, # iscale
                                   threshold=15,
                                   ):
        self.artdet_UseStructure = UseStructure # use structure detector or LocalSNR
        # for Hologic Selenia there is a part of the image (12 px) that should be ignored
        self.artdet_bkscale = bkscale # if given, then subtract Gaussian(bkscale) from filtered image
        self.artdet_fgscale = fgscale # Gaussian scale of features
        self.artdet_threshold = threshold # threshold for artefact or not
        
    def otsu(self,hist,bins):
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

    def RestrictROI(self,cs,otsu=True,borderpx=[0,0,0,0]):
        roipts = []
        cs.expert_roipts = []
        try:
            small_offset_px = 0 # set > 0 to make found region smaller with this offset from all edges

            widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
            heightpx = np.shape(cs.pixeldataIn)[1]
            # check histogram
            ypos = heightpx/2
            xmin0_px = borderpx[0]
            xmax0_px = widthpx-borderpx[1]-1
            
            hist,bins = np.histogram(cs.pixeldataIn[xmin0_px:xmax0_px:3, ypos], bins=500,density=False)
            thresh = self.otsu(hist,bins) # distinguish between background and foreground
            print '[otsu]x',thresh
            if not otsu:
                if cs.mustbeinverted:
                    hist[bins[:-1]>=thresh] = 0 # ignore background peak
                    thresh = bins[hist.argmax()]#1.05*bins[hist.argmax()] # find foreground peak
                else:
                    hist[bins[:-1]<=thresh] = 0 # ignore background peak
                    thresh = bins[hist.argmax()]#.95*bins[hist.argmax()] # find foreground peak
                    
            if cs.verbose:
                print "[RestrictRoi] htresh=",thresh
            xmin_px = xmin0_px
            xmax_px = xmax0_px
            for ix in range(xmin0_px,xmax0_px+1):
                if cs.mustbeinverted:
                    if(cs.pixeldataIn[ix,ypos]>0.1*thresh and cs.pixeldataIn[ix,ypos]<thresh):
                        xmin_px = ix+small_offset_px
                        break
                else:
                    if(cs.pixeldataIn[ix,ypos]<1.9*thresh and cs.pixeldataIn[ix,ypos]>thresh):
                        xmin_px = ix+small_offset_px
                        break
            for ix in reversed(range(xmin0_px,xmax0_px+1)):
                if cs.mustbeinverted:
                    if(cs.pixeldataIn[ix,ypos]>0.1*thresh and cs.pixeldataIn[ix,ypos]<thresh):
                        xmax_px = ix-small_offset_px
                        break
                else:
                    if(cs.pixeldataIn[ix,ypos]<1.9*thresh and cs.pixeldataIn[ix,ypos]>thresh):
                        xmax_px = ix-small_offset_px
                        break
            print '[restrict_x]',thresh,xmin_px,xmax_px

            xpos = widthpx/2
            ymin0_px = borderpx[2]# cs.border_offset_px
            ymax0_px = heightpx-borderpx[3] -1 #cs.border_offset_px -1
            hist,bins = np.histogram(cs.pixeldataIn[xpos,ymin0_px:ymax0_px:3], bins=500,density=False)
            thresh = self.otsu(hist,bins) #.95*hist[1][hist[0].argmax()]
            print '[otsu]y',thresh
            if not otsu:
                if cs.mustbeinverted:
                    hist[bins[:-1]>=thresh] = 0
                    thresh = bins[hist.argmax()]#1.05*bins[hist.argmax()]
                else:
                    hist[bins[:-1]<=thresh] = 0
                    thresh = bins[hist.argmax()]#.95*bins[hist.argmax()]
            if cs.verbose:
                print "[RestrictRoi] vtresh=",thresh
            ymin_px = ymin0_px
            ymax_px = ymax0_px
            for iy in range(ymin0_px,ymax0_px+1):
                if cs.mustbeinverted:
                    if(cs.pixeldataIn[xpos,iy]>0.1*thresh and cs.pixeldataIn[xpos,iy]<thresh):
                        ymin_px = iy+small_offset_px
                        break
                else:
                    if(cs.pixeldataIn[xpos,iy]<1.9*thresh and cs.pixeldataIn[xpos,iy]>thresh):
                        ymin_px = iy+small_offset_px
                        break
            for iy in reversed(range(ymin0_px,ymax0_px+1)):
                if cs.mustbeinverted:
                    if(cs.pixeldataIn[xpos,iy]>0.1*thresh and cs.pixeldataIn[xpos,iy]<thresh):
                        ymax_px = iy-small_offset_px
                        break
                else:
                    if(cs.pixeldataIn[xpos,iy]<1.9*thresh and cs.pixeldataIn[xpos,iy]>thresh):
                        ymax_px = iy-small_offset_px
                        break
            print '[restrict_y]',thresh,ymin_px,ymax_px
            roipts = [
                [xmin_px,ymin_px],
                [xmax_px,ymin_px],
                [xmax_px,ymax_px],
                [xmin_px,ymax_px]]
            cs.expert_roipts = [xmin_px,xmax_px, ymin_px,ymax_px]
            area_in = 1.*(xmax_px-xmin_px)*(ymax_px-ymin_px)
            area_all = 1.*(xmax0_px-xmin0_px)/(ymax0_px-ymin0_px)
            cs.expert_frac = area_in/area_all

            area_out = area_all-area_in
            mean_in = np.mean(cs.pixeldataIn[xmin_px:xmax_px:3, ymin_px:ymax_px:3])
            mean_all = np.mean(cs.pixeldataIn[xmin0_px:xmax0_px:3,ymin0_px:ymax0_px:3])
            mean_out = (area_all*mean_all-area_in*mean_in)/area_out
            cs.expert_inoutoverin = (mean_in-mean_out)/mean_in
        except:
            print "EXCEPTION!"
            pass
        return roipts

    def invertmaxval(self,cs):
        # Not stored properly in self.dcmfileIn.BitsStored
        dicomfields = [ ["0028,0101",  "Bits Stored"]]
        key = dicomfields[0][0]
        dicvalue = cs.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
        if dicvalue != "":
            return dicvalue
        else:
            return 0

    def NeedsCropping2(self,cs,mode='normi13'):
        """
        Simple test: If along any of the edges the value 2^15 -1 is reached or, a value higher than 1.4* center than yes
        32767 
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
        cmean = np.mean(cs.pixeldataIn[x0-wid_roi/2:x0+wid_roi/2,y0-hei_roi/2:y0+hei_roi/2])

        hist,bins = np.histogram(cs.pixeldataIn[::3,::3], bins=500,density=False)
        thresh = self.otsu(hist,bins) # distinguish between background and foreground
        total = np.sum(hist)
        hist[bins[:-1]>=thresh] = 0
        low = np.sum(hist)

        if mode == 'normi13':
            minfraction = 0.25
        elif mode == 'uniformity':
            minfraction = .025
        num = 5
        if cs.mustbeinverted:
            invertmax = 2**(self.invertmaxval(cs))-1
            edge_value = np.max([np.max(cs.pixeldataIn[:num,:]),
                                 np.max(cs.pixeldataIn[-num:,:]),
                                 np.max(cs.pixeldataIn[:,:num]),
                                 np.max(cs.pixeldataIn[:,-num:])])
            if edge_value > min(.95*invertmax -1.,1.4*cmean): 
                if 1.*(total-low)/total>minfraction:
                    print '[NeedsCropping2] Needs +cropping (%f).'%(1.*(total-low)/total)
                    return True 
        else:
            edge_value = np.min([np.min(cs.pixeldataIn[:num,:]),
                                 np.min(cs.pixeldataIn[-num,:]),
                                 np.min(cs.pixeldataIn[:,:num]),
                                 np.min(cs.pixeldataIn[:,-num:])])
            if edge_value < max(0.05*invertmax+1.,0.6*cmean): 
                if 1.*low/total>minfraction:
                    print '[NeedsCropping2] Needs -cropping (%f).'%(1.*low/total)
                    return True 

        print '[NeedsCropping2] Does not need cropping (%f).'%(1.*low/total)
        return False


    def NeedsCropping(self,cs):
        """
        Simple test: If the area of the automatic ROI is too small, the wrong phantom size is used,
        so cropping must be applied. Note that this means that QC is in fact invalid, since only part
        of the detector area can be checked.
        """
        cs_test = copy.deepcopy(cs)
        self.RestrictROI(cs_test)
        #print "[NeedsCropping] frac=",cs_test.expert_frac
        cs.expert_frac        = cs_test.expert_frac
        cs.expert_inoutoverin = cs_test.expert_inoutoverin
        #if cs_test.expert_frac < .75:
        #  return True
        #return False
        #print "[NeedsCropping] ratio=",cs_test.expert_inoutoverin
        if cs_test.expert_inoutoverin > .5:
            print "[NeedsCropping] frac,ratio=",cs_test.expert_frac,cs_test.expert_inoutoverin
            return True
        return False


    def Uniformity(self,cs,borderpx=[0,0,0,0]):
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
        if(cs.expertmode == True):
            print '[Uniformity] using expertmode'
            field_min_x_px = cs.expert_roipts[0]
            field_max_x_px = cs.expert_roipts[1]
            field_min_y_px = cs.expert_roipts[2]
            field_max_y_px = cs.expert_roipts[3]
        ## According to EUREF_DMWG_Protocol_2003 the areas need to be 4cm^2

        ## ROI definitions:

        ## 0    1
        ##
        ##   4  5
        ##
        ## 2    3

        width_roi_mm = 20.0 #width of roi in mm
        width_roi_px = int(width_roi_mm/pixel_spacing)  #width of roi in px
        height_roi_mm = 20.0 #height of roi in mm
        height_roi_px = int(height_roi_mm/pixel_spacing)  #height of roi in px

        cs.unif_rois = [] # format: x0,wid, yo,hei
        x0 = field_min_x_px
        y0 = field_min_y_px
        x1 = field_max_x_px-width_roi_px
        y1 = field_max_y_px-height_roi_px
        
        cs.unif_rois.append([x0,width_roi_px, y0,height_roi_px])
        cs.unif_rois.append([x1,width_roi_px, y0,height_roi_px])
        cs.unif_rois.append([x0,width_roi_px, y1,height_roi_px])
        cs.unif_rois.append([x1,width_roi_px, y1,height_roi_px])
        x0 = field_min_x_px+(field_max_x_px-field_min_x_px)/2
        y0 = field_min_y_px+(field_max_y_px-field_min_y_px)/2

        cs.unif_rois.append([x0-width_roi_px/2,width_roi_px, y0-height_roi_px/2,height_roi_px]) # ? or 0?
        cs.unif_rois.append([x1,width_roi_px, y0-height_roi_px/2,height_roi_px])

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
    def SplitMergeStructureDetector(self,inImage, bksigma = None,uiobject=None):
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

        #print "[SplitMergeStructureDetector] inImage first/last/rows",0,rows-1,rows

        # 1. start with an empty result LIST
        result = np.empty_like(inImage, dtype=float)
        overlap = 4*(5+2)+1 #4*(iscale+dscale)+1; filter is 4 sigma wide
        xparts = [ # first last keepfirst keeplast
                   [0,             rows/2 +overlap,  0,     rows/2-1],
                   [rows/2-overlap,rows,             rows/2,rows-1]
                   ]

        for xpart in xparts:
            pdCopy = (inImage[xpart[0]:xpart[1]:]).astype(float)
            #print "[SplitMergeStructureDetector] pdCopy first/last/rows",xpart[0],xpart[1]-1,np.shape(pdCopy)[0]
            pdProc = self.StructureDetector(pdCopy, None,uiobject)

            firstkeep = xpart[2]-xpart[0]
            lastkeep  = xpart[3]-xpart[0]
            offset    = xpart[2]-firstkeep
            #print "[SplitMergeStructureDetector] Restore first/last",firstkeep+offset,lastkeep+offset
            for i,row in enumerate(pdProc):
                if i<firstkeep or i>lastkeep:
                    continue
                result[i+offset] = row

        if bksigma is not None:
            blurRes = scind.gaussian_filter(result,bksigma,order=[0,0])
            result -= blurRes
            if uiobject:
                print "[SplitMergeStructureDetector] 5/4 trend removed"

        return result

    def StructureDetector(self,inImage, bksigma = None,uiobject=None):
        """
        Code is helaas slecht leesbaar geworden, maar gebruikte teveel geheugen.
        """
        #import resource
        #print "[Structure] resMB",1,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) # max usage so far

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
            #print "Cannot do progress"
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
                print "[StructureDetector] 5/4 trend removed"

        try:
            uiobject.pbar.endProgress()
        except:
            pass

        return Lxx

    def LocalSNR(self,pSrc, sigma,bksigma = None, uiobject=None):
        """
        Local SNR: [ mean(x,y) ]/stdev(x,y)
        Can be approximated as [ Gauss{I,sigma}] / sqrt[ Gauss{I-Gauss{I,sigma},sigma}^2]
        """
        if uiobject:
            uiobject.pbar.startProgress(3,"Calculating LocSNR")
        blurIm = scind.gaussian_filter(pSrc,sigma,order=[0,0])
        if uiobject:
            print "[LocalSNR] 1/4 blur'd"
        devIm = pSrc-blurIm
        if uiobject:
            uiobject.pbar.doProgress("dev'd")
            print "[LocalSNR] 2/4 dev'd"
        sdIm  = np.sqrt(scind.gaussian_filter(devIm**2,sigma,order=[0,0]))
        sdIm[sdIm<1.e-6]=1. # prevent div by zero
        if uiobject:
            uiobject.pbar.doProgress("std'd")
            print "[LocalSNR] 3/4 std'd"
        locnormIm = blurIm/sdIm
        if uiobject:
            print "[LocalSNR] 4/4 snr'd"
        if bksigma is not None:
            blurIm = scind.gaussian_filter(locnormIm,bksigma,order=[0,0])
            locnormIm -= blurIm
            if uiobject:
                print "[LocalSNR] 5/4 trend removed"
        if uiobject:
            uiobject.pbar.endProgress()
#        blurIm = scind.gaussian_filter(locnormIm,sigma,order=[2,0])**2+scind.gaussian_filter(locnormIm,sigma,order=[0,2])**2
#        return blurIm
        return locnormIm

#----------------------------------------------------------------------
    def Artefacts(self,cs,borderpx=[0,0,0,0],uiobject=None):
        error = True
        doUseStructureAlways = False
        print '[Artefacts]',borderpx
        # copy image, and cutoff black borders
        field_max_x_px = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        field_max_y_px = np.shape(cs.pixeldataIn)[1]
        field_min_x_px = 0
        field_min_y_px = 0
        if cs.expertmode:
            field_min_x_px = cs.expert_roipts[0]
            field_max_x_px = cs.expert_roipts[1]
            field_min_y_px = cs.expert_roipts[2]
            field_max_y_px = cs.expert_roipts[3]

        field_max_x_px -= borderpx[1] 
        field_max_y_px -= borderpx[3]
        field_min_x_px += borderpx[0]
        field_min_y_px += borderpx[2]

        pdCopy =  cs.pixeldataIn[field_min_x_px:field_max_x_px,field_min_y_px:field_max_y_px]

        if self.artdet_UseStructure: # use structure detector
            # Apart from spots, we also see lines appearing as large artefacts.
            # To detect the later, the slower "structure" is needed
            have_big_mem = True
            if have_big_mem:
                cs.art_image = self.StructureDetector(pdCopy,bksigma=self.artdet_bkscale,uiobject=uiobject) # remove background trend
            else:
                cs.art_image = self.SplitMergeStructureDetector(pdCopy,bksigma=self.artdet_bkscale,uiobject=uiobject) # remove background trend
                
        else:
            # We see only spots, so we do not need the slower "structure"
            cs.art_image = self.LocalSNR(pdCopy,sigma=self.artdet_fgscale,bksigma=self.artdet_bkscale,uiobject=uiobject)# remove background trend

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
        cs.art_threshold = self.artdet_threshold 
        if self.artdet_UseStructure:
            if not np.abs(cs.art_image).max()>self.artdet_threshold: # > above for Structure!
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs.art_image[::3, ::3], bins=[min(cs.art_image.min(),-self.artdet_threshold),-self.artdet_threshold,self.artdet_threshold,max(cs.art_image.max(),self.artdet_threshold)],density=False)
            frac = 1.*(hist[0]+hist[-1])/( (cwid/3)*(chei/3) ) # above for Structure!
        else:
            if not np.abs(cs.art_image).max()>self.artdet_threshold: #<
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs.art_image[::3, ::3], bins=[min(cs.art_image.min(),-self.artdet_threshold),-self.artdet_threshold,self.artdet_threshold,max(cs.art_image.max(),self.artdet_threshold)],density=False)
            frac = 1.*(hist[0]+hist[-1])/( (cwid/3)*(chei/3) )

        if frac > .1:
            if uiobject:
                uiobject.pbar.endProgress()
                uiobject.pbar.label.setText("Corrupt image!")
            print "[SNRArtefacts] ERROR! Corrupt image!",hist,self.artdet_threshold
            if uiobject:
                uiobject.pbar.label.setText("Finished.")
            return error

        if uiobject:
            uiobject.pbar.doProgress("clustering artefacts...")
        cca = wadwrapper_lib.connectedComponents()
        if self.artdet_UseStructure:
            if borderpx[0]>0:
                cs.art_image[:(borderpx[0]+1),:] = self.artdet_threshold-1
            if borderpx[1]>0:
                cs.art_image[-borderpx[1]:,:]    = self.artdet_threshold-1
            if borderpx[2]>0:
                cs.art_image[:,:(borderpx[2]+1)] = self.artdet_threshold-1
            if borderpx[3]>0:
                cs.art_image[:,-borderpx[3]:]    = self.artdet_threshold-1
            cca_in = np.abs(cs.art_image)>self.artdet_threshold
        else:
            if borderpx[0]>0:
                cs.art_image[:(borderpx[0]+1),:] = -self.artdet_threshold+1
            if borderpx[1]>0:
                cs.art_image[-borderpx[1]:,:]    = -self.artdet_threshold+1
            if borderpx[2]>0:
                cs.art_image[:,:(borderpx[2]+1)] = -self.artdet_threshold+1
            if borderpx[3]>0:
                cs.art_image[:,-borderpx[3]:]    = -self.artdet_threshold+1
            cca_in = np.abs(cs.art_image)>self.artdet_threshold
        cc_art_image,nclusters = cca.run(cca_in) 
        for a in range(1,nclusters+1):
            clusters.append(cca.indicesOfCluster(a))

        nc = len(clusters)
        if cs.verbose:
            print "[Artefacts] found",npix ," artefact pixels in",nc,"clusters"

        if uiobject:
            uiobject.pbar.endProgress()
        if cs.verbose:
            print "[Artefacts] Found",len(clusters),"clusters"

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
                print "[Artefacts]...size",len(c)

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
        for t in range(-(thick-1)/2,(thick+1)/2):
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
