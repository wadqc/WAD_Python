# -*- coding: utf-8 -*-
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

"""
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

Changelog:
    20161220: removed class variables; removed testing stuff
    20160902: sync with wad2.0
    20150522: Force scipy version at least 0.10.1 to avoid problems with PIL/Pillow mixing
    20150430: allow crash if no PIL.ImageDraw available
    20150428: prefer ImageDraw from PIL to fix error
    20150420: replaced connectedComponents with scipy's version
    20150213: fixed missing scannername; no matplots if not in gui
    20150205: DICOMInfo according to general interface
    20141024: DoseRatio for Hologic Dimensions
    20141006: all DICOM read through wadwrapper_lib
    20140501: Annotated image downscaled since it is a jpeg anyway
    20140430: SplitMerge in two parts for L50 StructureDetector to save mem
    20140429: Always call L50 contrast; save jpg of annotated contrast; fixed roisize=0
    20140428: Also BK removal for Structure Artefact Detector; will keep both detectors Structure for L50 and SNR for Selenia; fix L50 contrast
    20140425: L50 contrast added
    20140424: New clustering, reverted artefact detection, gui visibility
    20140416: More sensitive doseratio
    20140402: Fix div zero in LocalSNR; added protocol name to dicominfo
    20140328: Compatibility with WAD; crop detection; new artefact detection
    20131011: Finished linearity
    20131010: FFU calc of rad10 and rad20 by Euclidan distance transform
    20131009: Finished SNR; finished ArtLevel; finish FloodField Uniformity
"""
__version__ = '20161220'
__author__ = 'aschilham'

import numpy as np
from math import pi
import copy
import scipy.ndimage as scind

try:
    # wad2.0 runs each module stand alone
    import QCMammo_constants as lit
except ImportError:
    from . import QCMammo_constants as lit
    
import matplotlib.pyplot as plt
from PIL import Image # image from pillow is needed
from PIL import ImageDraw # imagedraw from pillow is needed, not pil
import scipy.misc
# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
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
        # wad1.0 solutions failed, try wad2.0 from system package wad_qc
        from wad_qc.modulelibs import wadwrapper_lib

class MammoStruct:
    def __init__ (self, dcmInfile, pixeldataIn):
        self.verbose = False

        # input image
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn

        # for matlib plotting
        self.hasmadeplots = False

        # uniformity
        self.means = []
        self.stdevs = []
        self.unif_pct = -1
        self.snr_hol = -1
        self.unif_rois = [] # xy roi definitions # format: x0,wid, y0,hei

        # dose ratio
        self.doseratio = -1

        # artefacts
        self.art_clusters = []
        self.art_image = None
        self.art_borderpx_lrtb = [0, 0, 0, 0] # for selenia there is a part of the image that should be ignored
        self.art_threshold = -1
        self.art_rois = [] # x0,y0,rad

        # working with a smaller FOV
        self.expertmode = False
        self.expert_roipts = []
        self.expert_frac = -1
        self.expert_inoutoverin = -1

        # contrast L50
        self.contrast_rois = [] # x0,y0,rad
        self.contrast_snr = []
        self.contrast_mean = []
        self.contrast_sd = []

        # identifiers
        self.filtername  = lit.stUnknown
        self.scannername = lit.stUnknown
        self.DetermineScannerID()

    def DetermineScannerID(self):
        """
        Actually only needed to determine if L50 contrast module should run
        """
        error = True
        self.scannername = lit.stUnknown
        dicomvalue = wadwrapper_lib.readDICOMtag("0008,1090", self.dcmInfile, imslice=0) #Manufacturer's Model Name
        dicomvalue = str(dicomvalue).lower()
        if dicomvalue.find("l50")>-1:
            self.scannername = lit.stL50
            error = False
        elif dicomvalue.find("lorad selenia")>-1:
            self.scannername = lit.stSelenia
            error = False
        elif dicomvalue.find("dimensions")>-1:
            self.scannername = lit.stDimensions
            error = False
        elif dicomvalue.find("affirm prone")>-1:
            self.scannername = lit.stAffirmProne
            error = False
        return error


class Mammo_QC:
    def __init__(self,guimode=False):
        self.qcversion = __version__
        self.guimode = guimode
        pass

    def readDICOMtag(self,cs_mam,key,imslice=0): # slice=2 is image 3
        value = wadwrapper_lib.readDICOMtag(key,cs_mam.dcmInfile,imslice)
        return value

    def pixDim(self,dcmInfile,axis=0):
        try:
            pix_mm = dcmInfile.PixelSpacing[axis]
        except:
            pix_mm = dcmInfile.ImagerPixelSpacing[axis]
        return pix_mm

    def otsu(self,hist,bins):
        currentMax = 0
        threshold = 0
        sumTotal, sumForeground, sumBackground = 0, 0, 0
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
            meanB = sumBackground / weightBackground
            meanF = (sumTotal - sumBackground) / weightForeground

            # Calculate variance between classes
            varBetween = weightBackground*weightForeground
            varBetween *= (meanB-meanF)*(meanB-meanF)

            # Check if the variance between classes is greater than the current best
            if(varBetween > currentMax):
                currentMax = varBetween
                threshold = b

        return threshold

    def RestrictROI(self,cs_mam):
        roipts = []
        cs_mam.expert_roipts = []
        try:
            pixel_spacing  = self.pixDim(cs_mam.dcmInfile) #  self.dcmInfile.PixelSpacing[0] ## spacing in mm
            border_offset_mm = 10.0 #this is the
            border_offset_px = int(border_offset_mm/pixel_spacing)
            small_offset_px = 10

            widthpx = np.shape(cs_mam.pixeldataIn)[0] ## width/height in pixels
            heightpx = np.shape(cs_mam.pixeldataIn)[1]
            # check histogram
            ypos = int(heightpx/2)
            xmin_px = 0
            xmax_px = widthpx-1
            hist,bins = np.histogram(cs_mam.pixeldataIn[::3, ypos], bins=500,density=False)
            thresh = self.otsu(hist,bins) # distinguish between background and forground
            hist[bins[:-1]<=thresh] = 0 # ignore background peak
            thresh = .95*bins[hist.argmax()] # find foreground peak
            if cs_mam.verbose:
                print("[RestrictRoi] htresh=",thresh)
            for ix in range(0,widthpx):
                if(cs_mam.pixeldataIn[ix,ypos]>thresh):
                    xmin_px = ix+small_offset_px
                    break
            for ix in reversed(range(0,widthpx)):
                if(cs_mam.pixeldataIn[ix,ypos]>thresh):
                    xmax_px = ix-small_offset_px
                    break
            xpos = int(widthpx/2)
            ymin_px = border_offset_px
            ymax_px = heightpx-border_offset_px -1
            hist,bins = np.histogram(cs_mam.pixeldataIn[xpos,border_offset_px:heightpx-border_offset_px:3], bins=500,density=False)
            thresh = self.otsu(hist,bins) #.95*hist[1][hist[0].argmax()]
            hist[bins[:-1]<=thresh] = 0
            thresh = .95*bins[hist.argmax()]
            if cs_mam.verbose:
                print("[RestrictRoi] vtresh=",thresh)
            for iy in range(border_offset_px,heightpx-border_offset_px):
                if(cs_mam.pixeldataIn[xpos,iy]>thresh):
                    ymin_px = iy+small_offset_px
                    break
            for iy in reversed(range(border_offset_px,heightpx-border_offset_px)):
                if(cs_mam.pixeldataIn[xpos,iy]>thresh):
                    ymax_px = iy-small_offset_px
                    break
            roipts = [
                [xmin_px,ymin_px],
                [xmax_px,ymin_px],
                [xmax_px,ymax_px],
                [xmin_px,ymax_px ]]
            cs_mam.expert_roipts = [xmin_px,xmax_px, ymin_px,ymax_px]
            cs_mam.expert_frac = 1.*(xmax_px-xmin_px)*(ymax_px-ymin_px)/(widthpx*(heightpx-2*border_offset_px))
            area_in = 1.*(xmax_px-xmin_px)*(ymax_px-ymin_px)
            area_all = 1.*(widthpx*(heightpx-2*border_offset_px))
            area_out = area_all-area_in
            mean_in = np.mean(cs_mam.pixeldataIn[xmin_px:xmax_px:3, ymin_px:ymax_px:3])
            mean_all = np.mean(cs_mam.pixeldataIn[::3,border_offset_px:heightpx-border_offset_px:3])
            mean_out = (area_all*mean_all-area_in*mean_in)/area_out
            cs_mam.expert_inoutoverin = (mean_in-mean_out)/mean_in
        except:
            print("EXCEPTION!")
            pass
        return roipts

    def NeedsCropping(self,cs_mam):
        """
        Simple test: If the area of the automatic ROI is too small, the wrong phantom size is used,
        so cropping must be applied. Note that this means that QC is in fact invalid, since only part
        of the detector area can be checked.
        """
        cs_test = copy.deepcopy(cs_mam)
        self.RestrictROI(cs_test)
#        print("[NeedsCropping] frac=",cs_test.expert_frac)
        cs_mam.expert_frac        = cs_test.expert_frac
        cs_mam.expert_inoutoverin = cs_test.expert_inoutoverin
#        if cs_test.expert_frac < .75:
#            return True
#        return False
#        print("[NeedsCropping] ratio=",cs_test.expert_inoutoverin)
        if cs_test.expert_inoutoverin > .5:
            print("[NeedsCropping] frac,ratio=",cs_test.expert_frac,cs_test.expert_inoutoverin)
            return True
        return False


    def L50Contrast(self,cs_mam):
        """
        Use the six contrast inserts of the Microdose phantom to calculate contrast
        """
        if cs_mam.scannername != lit.stL50:
            return False

        error = True
        pixel_spacing_x = self.pixDim(cs_mam.dcmInfile,0)## spacing in mm
        pixel_spacing_y = self.pixDim(cs_mam.dcmInfile,1)##cs_mam.dcmInfile.PixelSpacing[1]

        if pixel_spacing_x == pixel_spacing_y :
            pixel_spacing = pixel_spacing_x

        hei = np.shape(cs_mam.pixeldataIn)[1]
        ## ROI definitions:
        rad_px = int(8./pixel_spacing)
        pos_x = int(110./pixel_spacing)
        pos_y = int(125./pixel_spacing)
        dy = int(54./pixel_spacing)
        dx = int(45./pixel_spacing)

        yoff = 0
        if self.NeedsCropping(cs_mam):
            roipts = self.RestrictROI(cs_mam)
            yoff = hei
            for rp in roipts:
                yoff = min(yoff,rp[1])
            if yoff == hei:
                yoff = 0
            else:
                yoff += 420 # empicallly if not placed properly

        tgxy = []
        bgxy = []
        for j in range(2):
            for i in range(3):
                x0 = pos_x+j*dx
                y0 = pos_y+i*dy+yoff
                ry = min(y0+2*rad_px,hei-1)-y0
                if ry> rad_px: # not placed properly, so lower part is missing
                    partimage = (cs_mam.pixeldataIn[x0-2*rad_px:x0+2*rad_px,y0-ry:y0+ry]).astype(float)
                    error,shiftxypx = self.FindPinShift2D(cs_mam,partimage)
                    tgxy.append([x0+shiftxypx[0],y0+shiftxypx[1]])
                    if j == 0:
                        bgxy.append([x0+shiftxypx[0]-3*rad_px,y0+shiftxypx[1]])
                    else:
                        bgxy.append([x0+shiftxypx[0]+3*rad_px,y0+shiftxypx[1]])

        # first avg and sd
        tg_mean = []
        tg_sd = []
        np.seterr(over='raise')
        for rxy in tgxy:
            sumv = np.float32(0.)
            sumv2 = np.float32(0.)
            count = 0
            for ky in range(-rad_px,rad_px):
                for kx in range(-rad_px,rad_px):
                    if kx**2+ky**2<rad_px**2 :
                        val = np.float32(cs_mam.pixeldataIn[rxy[0]+kx,rxy[1]+ky])
                        sumv += val
                        try:
                            sumv2 += val*val
                        except:
                            print("[L50Contrast] ERROR!",type(val),type(sumv),type(sumv2),val,sumv,sumv2,count)

                        count += 1
            if count>0:
                tg_mean.append(sumv/count)
                tg_sd.append( np.sqrt( (count*sumv2-sumv*sumv)/count/(count-1.) ) )
            else:
                tg_mean.append(0)
                tg_sd.append(0)

        bg_mean = []
        bg_sd = []
        np.seterr(over='raise')
        for rxy in bgxy:
            sumv = np.float32(0.)
            sumv2 = np.float32(0.)
            count = 0
            for ky in range(rad_px):
                for kx in range(0,rad_px):
                    if (kx-rad_px)**2+(ky-rad_px)**2<rad_px**2 :
                        val = np.float32(cs_mam.pixeldataIn[rxy[0]+kx,rxy[1]+ky])
                        sumv += val
                        try:
                            sumv2 += val*val
                        except:
                            print("[L50Contrast] ERROR!",type(val),type(sumv),type(sumv2),val,sumv,sumv2,count)

                        count += 1
            if count>0:
                bg_mean.append(sumv/count)
                bg_sd.append( np.sqrt( (count*sumv2-sumv*sumv)/count/(count-1.) ) )
            else:
                bg_mean.append(0)
                bg_sd.append(0)

        cs_mam.contrast_snr = []
        cs_mam.contrast_mean = []
        cs_mam.contrast_sd = []
        for k in range(len(tg_mean)):
            sd = np.sqrt(0.5*(tg_sd[k]**2.+bg_sd[k]**2.))
            if sd <1e-6:
                sd = 1.
            cs_mam.contrast_snr.append(np.abs(tg_mean[k]-bg_mean[k])/sd)
            cs_mam.contrast_mean.append(tg_mean[k])
            cs_mam.contrast_sd.append(tg_sd[k])

        cs_mam.contrast_rois = []
        for tr in tgxy:
            cs_mam.contrast_rois.append([tr[0],tr[1], rad_px])

        return error

    def FindPinShift2D(self,cs,pixeldata):
        error = True

        shiftxypx = [0,0]
        """
        Needs 2d input pixeldata. If 3D, then input pixeldata[slice]
        Concept:
            1. blur with object dependend scale
            2. Find bgval=average in border
            3. Find fgval=average center

        """
        if(np.shape(np.shape(pixeldata))[0]!=2):
            print("[FindCenterShift2D] Error, called with non-2D data!")
            return error

        wid = np.shape(pixeldata)[0] ## width/height in pixels
        hei = np.shape(pixeldata)[1]
        borderpx = int(max(2,min(wid,hei)/20))

        dscale = 25.0

        # 1. blur object
        blurIm = scind.gaussian_filter(pixeldata, sigma=dscale)

        # 2. Find bk value and fg value:
        bkval = np.average((np.average(pixeldata[0:borderpx,:]),
                            np.average(pixeldata[-borderpx:,:]),
                            np.average(pixeldata[:,0:borderpx,]),
                            np.average(pixeldata[:,-borderpx:])
                            ))
        fgval = np.average(pixeldata[int(wid/2)-borderpx:int(wid/2)+borderpx,int(hei/2)-borderpx:int(hei/2)+borderpx])

        thresh = (fgval*2.+bkval)/3.
        if cs.verbose:
            print("[FindCenterShift] bpix/bkval/fgval/thresh",borderpx,bkval,fgval,thresh)
        if bkval>fgval:
            lowmode = True
        else:
            lowmode = False

        # 2. only look for voxels below/above a threshold (eg. look for air)
        arrayLow = []
        for iy in range(0,hei):
            arrayLow.append(0)
            for ix in range(0,wid):
                if(lowmode):
                    if(blurIm[ix,iy]<thresh):
                        arrayLow[iy] += 1
                else:
                    if(blurIm[ix,iy]>thresh):
                        arrayLow[iy] += 1

        if(cs.verbose and self.guimode):
            plt.figure()
            plt.plot(arrayLow)
            plt.title("vertical ")
            cs.hasmadeplots = True

        # 2.1 find left first pos without voxels below threshLow
        minLowId = 0
        minLow = arrayLow[minLowId]
        for iy in range(0,hei):
            if(arrayLow[iy]>minLow):
                minLowId = iy
                minLow = arrayLow[minLowId]
        shiftxypx[1] = minLowId
        if(cs.verbose):
            print("[FindCenterShift] vertical ",minLowId)
        # 2.2 find right first pos without voxels below threshLow
        minLowId = hei-1
        minLow = arrayLow[minLowId]
        for iy in reversed(range(0,hei)):
            if(arrayLow[iy]>minLow):
                minLowId = iy
                minLow = arrayLow[minLowId]
        # 2.3 mid is halfway left and right pos
        shiftxypx[1] =(shiftxypx[1]+minLowId-(hei-1))/2
        if(cs.verbose):
            print("[FindCenterShift] vertical ",minLowId,(hei-1)/2,shiftxypx[1])
        # repeat for horizontal
        arrayLow = []
        for ix in range(0,wid):
            arrayLow.append(0)
            for iy in range(0,hei):
                if(lowmode):
                    if(blurIm[ix,iy]<thresh):
                        arrayLow[ix] += 1
                else:
                    if(blurIm[ix,iy]>thresh):
                        arrayLow[ix] += 1

        if(cs.verbose and self.guimode):
            plt.figure()
            plt.plot(arrayLow)
            plt.title("horizontal ")
            cs.hasmadeplots = True

        # 2.1 find left first pos without voxels below threshLow
        minLowId = 0
        minLow = arrayLow[minLowId]
        for iy in range(0,wid):
            if(arrayLow[iy]>minLow):
                minLowId = iy
                minLow = arrayLow[minLowId]
        shiftxypx[0] = minLowId
        if(cs.verbose):
            print("[FindCenterShift] horizontal ",minLowId)

        # 2.2 find right first pos without voxels below threshLow
        minLowId = wid-1
        minLow = arrayLow[minLowId]
        for iy in reversed(range(0,wid)):
            if(arrayLow[iy]>minLow):
                minLowId = iy
                minLow = arrayLow[minLowId]
        # 2.3 mid is halfway left and right pos
        shiftxypx[0] =(shiftxypx[0]+minLowId-(wid-1))/2
        if(cs.verbose):
            print("[FindCenterShift] horizontal ",minLowId,(wid-1)/2,shiftxypx[0])

        error = False
        return error,shiftxypx


    def Uniformity(self,cs_mam):
        """
        Calculation of non-uniformity

        Defined in EUREF_DMWG_Protocol_2003 as max deviation between avg pix value in five ROIs
        """
        error = True

        width_array = np.shape(cs_mam.pixeldataIn)[0] ## width/height in pixels
        height_array = np.shape(cs_mam.pixeldataIn)[1]

        pixel_spacing_x = self.pixDim(cs_mam.dcmInfile,0)## spacing in mm
        pixel_spacing_y = self.pixDim(cs_mam.dcmInfile,1)##cs_mam.dcmInfile.PixelSpacing[1]

        if pixel_spacing_x == pixel_spacing_y :
            pixel_spacing = pixel_spacing_x

        field_max_x_px = width_array ## total field width/height in px
        field_max_y_px = height_array
        field_min_x_px = 0
        field_min_y_px = 0
        if(cs_mam.expertmode == True):
            field_min_x_px = cs_mam.expert_roipts[0]
            field_max_x_px = cs_mam.expert_roipts[1]
            field_min_y_px = cs_mam.expert_roipts[2]
            field_max_y_px = cs_mam.expert_roipts[3]

        border_offset_mm = 10.0 #this is the
        border_offset_px = int(border_offset_mm/pixel_spacing)
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

        cs_mam.unif_rois = [] # format: x0,wid, yo,hei
        x0 = field_min_x_px+border_offset_px
        y0 = field_min_y_px+border_offset_px
        x1 = field_max_x_px-width_roi_px-border_offset_px
        y1 = field_max_y_px-height_roi_px-border_offset_px
        cs_mam.unif_rois.append([x0,width_roi_px, y0,height_roi_px])
        cs_mam.unif_rois.append([x1,width_roi_px, y0,height_roi_px])
        cs_mam.unif_rois.append([x0,width_roi_px, y1,height_roi_px])
        cs_mam.unif_rois.append([x1,width_roi_px, y1,height_roi_px])
        x0 = int(field_min_x_px+(field_max_x_px-field_min_x_px)/2)
        y0 = int(field_min_y_px+(field_max_y_px-field_min_y_px)/2)
        if cs_mam.scannername == lit.stL50:
            cs_mam.unif_rois.append([x0,width_roi_px, y0-int(3.5*border_offset_px),height_roi_px])
        else:
            cs_mam.unif_rois.append([x0,width_roi_px, y0,height_roi_px])
        cs_mam.unif_rois.append([x1,width_roi_px, y0,height_roi_px])

        cs_mam.means  = []
        cs_mam.stdevs = []
        for r in cs_mam.unif_rois:
            arr = cs_mam.pixeldataIn[r[0]:(r[0]+r[1]),r[2]:(r[2]+r[3])]
            cs_mam.means.append(np.mean(arr))
            cs_mam.stdevs.append(np.std(arr))

        maxdev = abs(cs_mam.means[0]-cs_mam.means[4])
        for i in range(1,4):
            maxdev = max(maxdev,abs(cs_mam.means[i]-cs_mam.means[4]))

        cs_mam.unif_pct = 100.*maxdev/cs_mam.means[4]
        cs_mam.snr_hol = cs_mam.means[5]/cs_mam.stdevs[5]
        if cs_mam.verbose:
            print("[Uniformity] maxdev="+str(maxdev)+" unif="+str(cs_mam.unif_pct)+" snr="+str(cs_mam.snr_hol))
        error = False
        return error

    def DoseRatio(self,cs_mam):
        """
        New dose ratio. The idea is to produce a number that says something about the dose settings (calibration?)
        without actually measuring the dose (with an external probe), which is not influenced by unimportant
        common mistakes in producing the images (like an improper compression).
        Assume dose at detector (distance R) is constant D(R). Than the dose D(R-t) at entrance point
        of breast of compressed thickness t is D(R-t) = D(R) . [R/(R-t)]^2 . exp[mu.t]
        Next, the Dose will be scaled by uAs. This all results in something which is not entirely correct, but gives
        a good predictor of the dose settings.
        """
        error = True
        # add EntranceDose info
        # Determine MOLYBDENUM or RHODIUM: different offset
        dicomfields = [
            ["0018,1110",  "DistanceSourceToDetector"],
            ["0018,11A0",  "BodyPartThickness"],
            ["0040,8302",  "EntranceDose_(mGy)"],
            ["0018,1153",  "muAs"],
            ["0018,7050",  "FilterMaterialLT"],
            ["0008,1090",  "ModelName"]
        ]

        # calculate distance = L0-Thick
        key = dicomfields[0][0]
        L0 = self.readDICOMtag(cs_mam,key)
        key = dicomfields[1][0]
        thick = self.readDICOMtag(cs_mam,key)
        key = dicomfields[2][0]
        entdose = self.readDICOMtag(cs_mam,key)
        key = dicomfields[3][0]
        umAs = self.readDICOMtag(cs_mam,key)

        relR2 = (L0/(L0-thick))**2.
        key = dicomfields[5][0]
        model = self.readDICOMtag(cs_mam,key)

        # Determine MOLYBDENUM or RHODIUM: different offset
        key = dicomfields[4][0]
        filt = self.readDICOMtag(cs_mam,key)
        if filt == "RHODIUM":
            cs_mam.filtername = "RH"
            if "Lorad Selenia" in model:
                # for Lorad Selenia
                slope  = 0.02177294034021
                offset = -8.80671341620691
            elif "Selenia Dimensions" in model:
                # for Hologic Dimensions
                slope  = 0.0276049841785717
                offset = -10.225707588108031
            else:
                print("[DoseRatio] Error! Unknown device model "+model)
                cs_mam.doseratio = -1
                return False
        elif filt == "MOLYBDENUM":
            # for Lorad Selenia
            slope  = 0.02076894376829
            offset = -8.37788116074648
            cs_mam.filtername = "MO"
        elif filt == "ALUMINUM":
            # for Philips L30/L50
            slope  = 0.03225282432099
            offset = -8.91017235713734
            key = ["0019,10c2",  "EnergyComponent"]
            cs_mam.filtername =  self.readDICOMtag(cs_mam,key[0])
            key = ["0008,0068",  "PresentationIntent"]
            val = self.readDICOMtag(cs_mam,key[0])
            if "PROCESSING" in val:
                cs_mam.filtername += "PROC"
            elif "PRESENTATION" in val:
                cs_mam.filtername += "PRES"
        elif filt == "SILVER":
            # for Hologic Dimensions
            slope  = 0.0272018880924507
            offset = -9.738774100294515
            cs_mam.filtername = "AG"
        else:
            print("[DoseRatio] Error! Unknown filtermaterial "+(filt))
            cs_mam.doseratio = -1
            cs_mam.filtername = filt
            return False

        pred = np.exp( (slope*thick+offset)*relR2 )*umAs

        # Divide EntranceDose by Disc Integral
        cs_mam.doseratio = entdose/pred
        error = False
#        print(model,cs_mam.filtername,L0,thick,relR2,umAs,pred,entdose)
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

#        print("[SplitMergeStructureDetector] inImage first/last/rows",0,rows-1,rows)

        # 1. start with an empty result LIST
        result = np.empty_like(inImage, dtype=float)
        overlap = 4*(5+2)+1 #4*(iscale+dscale)+1; filter is 4 sigma wide
        rows2 = int(rows/2)
        xparts = [ # first last keepfirst keeplast
                   [0,             rows2 +overlap,  0,     rows2-1],
                   [rows2-overlap,rows,             rows2,rows-1]
                ]

        for xpart in xparts:
            pdCopy = (inImage[xpart[0]:xpart[1]:]).astype(float)
#            print("[SplitMergeStructureDetector] pdCopy first/last/rows",xpart[0],xpart[1]-1,np.shape(pdCopy)[0])
            pdProc = self.StructureDetector(pdCopy, None,uiobject)

            firstkeep = xpart[2]-xpart[0]
            lastkeep  = xpart[3]-xpart[0]
            offset    = xpart[2]-firstkeep
#            print("[SplitMergeStructureDetector] Restore first/last",firstkeep+offset,lastkeep+offset)
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

    def StructureDetector(self,inImage, bksigma = None,uiobject=None):
        """
        Code is helaas slecht leesbaar geworden, maar gebruikte teveel geheugen.
        """
#        import resource
#        print("[Structure] resMB",1,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) # max usage so far)

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
        #try:
        try:
            uiobject.pbar.doProgress("...calculating Lx**2")
        except:
            pass
        #except:
        #  print("Cannot do progress")
        #  pass
        Lx = (scind.gaussian_filter(pdCopy,dscale,order=[1,0]))
        Lxx = Lx**2
        try:
            uiobject.pbar.doProgress("...calculating Ly**2")
        except:
            pass
        Ly = (scind.gaussian_filter(pdCopy,dscale,order=[0,1]))
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
        Lxx = scind.gaussian_filter(Lxx,iscale,order=[0,0])
        try:
            uiobject.pbar.doProgress("...smoothing Ly**2")
        except:
            pass
        Lyy = scind.gaussian_filter(Lyy,iscale,order=[0,0])
        try:
            uiobject.pbar.doProgress("...smoothing Lx*Ly")
        except:
            pass
        Lx = scind.gaussian_filter(Lx,iscale,order=[0,0]) # Lxy = scind.gaussian_filter(Lxy,iscale,order=[0,0]); 3 in use

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
            Lyy = scind.gaussian_filter(Lxx,bksigma,order=[0,0])
            Lxx -= Lyy
            if uiobject:
                print("[StructureDetector] 5/4 trend removed")

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
#        blurIm = scind.gaussian_filter(locnormIm,sigma,order=[2,0])**2+scind.gaussian_filter(locnormIm,sigma,order=[0,2])**2
#        return blurIm
        return locnormIm

#----------------------------------------------------------------------
    def Artefacts(self,cs_mam,uiobject=None):
        error = True
        doUseStructureAlways = False

        # copy image, and cutoff black borders
        field_max_x_px = np.shape(cs_mam.pixeldataIn)[0] - cs_mam.art_borderpx_lrtb[1] ## width/height in pixels
        field_max_y_px = np.shape(cs_mam.pixeldataIn)[1] - cs_mam.art_borderpx_lrtb[3]
        field_min_x_px = cs_mam.art_borderpx_lrtb[0]
        field_min_y_px = cs_mam.art_borderpx_lrtb[2]
        if cs_mam.expertmode:
            field_min_x_px = cs_mam.expert_roipts[0]
            field_max_x_px = cs_mam.expert_roipts[1]
            field_min_y_px = cs_mam.expert_roipts[2]
            field_max_y_px = cs_mam.expert_roipts[3]

        pdCopy =  cs_mam.pixeldataIn[field_min_x_px:field_max_x_px,field_min_y_px:field_max_y_px]

        if cs_mam.scannername == lit.stL50:
            # Apart from spots, we also see lines appearing as large artefacts.
            # To detect the later, the slower "structure" is needed
            dscale = 25.0
###            import resource
###            print("[Structure] resMB",1,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) # max usage so far)
#            cs_mam.art_image = self.StructureDetector(pdCopy,bksigma=dscale,uiobject=uiobject) # remove background trend
            cs_mam.art_image = self.SplitMergeStructureDetector(pdCopy,bksigma=dscale,uiobject=uiobject) # remove background trend
###            print("[Structure] resMB",1,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) # max usage so far)
###            return
            thresh = 3000 # with bk removal, lower value required
##            cs_mam.art_image = self.StructureDetector(pdCopy,uiobject=uiobject) # do not remove background trend
##            thresh = 4000
        else:
            if doUseStructureAlways:
                #cs_mam.art_image = self.StructureDetector(pdCopy,uiobject=uiobject)
                #thresh = 5
                dscale = 25.0
                cs_mam.art_image = self.StructureDetector(pdCopy,bksigma=dscale,uiobject=uiobject) # remove background trend
                thresh = 10 # less sensitive 5
            else:
                # We see only spots, so we do not need the slower "structure"
                iscale = 5.0
                dscale = 25.0
                cs_mam.art_image = self.LocalSNR(pdCopy,sigma=iscale,bksigma=dscale,uiobject=uiobject)# remove background trend
                #cs_mam.art_image = self.LocalSNR(pdCopy,sigma=iscale,uiobject=uiobject)# do not remove background trend
                #return error
                thresh = -15. # more sensitive -20.

        """
        find clusters of pixels BELOW a certain threshold
        """
        # For L50, first exclude contrast plugs:
        if cs_mam.scannername == lit.stL50:
            self.L50Contrast(cs_mam)
            pixel_spacing = self.pixDim(cs_mam.dcmInfile,0)## spacing in mm
            rad_px = int(8.5/pixel_spacing)
            wid = np.shape(cs_mam.art_image)[0] ## width/height in pixels
            hei = np.shape(cs_mam.art_image)[1] ## width/height in pixels

            yoff = 0
            if self.NeedsCropping(cs_mam):
                roipts = self.RestrictROI(cs_mam)
                yoff = hei
                for rp in roipts:
                    yoff = min(yoff,rp[1])
                if yoff == hei:
                    yoff = 0

            rad_in = rad_px-5*5 # iscale = 5
            rad_out = rad_px + 5*5
            rin2 = rad_in**2
            rout2 = rad_out**2
            for rxy in cs_mam.contrast_rois:
                x0 = rxy[0]
                y0 = rxy[1]-yoff
                for ky in range(-min(rad_out,y0),min(rad_out,hei-y0-1)):
                    for kx in range(-min(rad_out,x0),min(rad_out,wid-x0-1)):
                        d2 = kx**2+ky**2
                        if rin2 < d2 < rout2:
                            cs_mam.art_image[x0+kx,y0+ky] = 0

        if uiobject:
            uiobject.pbar.startProgress(2,"checking histogram...")
        cwid = np.shape(cs_mam.art_image)[0]
        chei = np.shape(cs_mam.art_image)[1]
        clusters = []
        npix = 0
        # check histogram
        cs_mam.art_clusters = []
        cs_mam.art_threshold = thresh
        if cs_mam.scannername == lit.stL50 or doUseStructureAlways:
            if not cs_mam.art_image.max()>thresh: # above for Structure!
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs_mam.art_image[::3, ::3], bins=[thresh,cs_mam.art_image.max()],density=False)
            frac = 1.*hist[0]/( (cwid/3)*(chei/3) ) # above for Structure!
        else:
            if not cs_mam.art_image.min()<thresh:
                error = False
                if uiobject:
                    uiobject.pbar.label.setText("Finished.")
                return error
            hist,bins = np.histogram(cs_mam.art_image[::3, ::3], bins=[cs_mam.art_image.min(),thresh],density=False)
            frac = 1.*hist[0]/( (cwid/3)*(chei/3) )

        if frac > .1:
            if uiobject:
                uiobject.pbar.endProgress()
                uiobject.pbar.label.setText("Corrupt image!")
            print("[SNRArtefacts] ERROR! Corrupt image!",hist)
            if uiobject:
                uiobject.pbar.label.setText("Finished.")
            return error

        if uiobject:
            uiobject.pbar.doProgress("clustering artefacts...")
        cca = wadwrapper_lib.connectedComponents()
        if cs_mam.scannername == lit.stL50 or doUseStructureAlways:
            cca_in = cs_mam.art_image>thresh
        else:
            cca_in = cs_mam.art_image<thresh
        cc_art_image,nclusters = cca.run(cca_in) 
        for a in range(1,nclusters+1):
            clusters.append(cca.indicesOfCluster(a))

        nc = len(clusters)
        if cs_mam.verbose:
            print("[Artefacts] found",npix ," artefact pixels in",nc,"clusters")

        if uiobject:
            uiobject.pbar.endProgress()
        if cs_mam.verbose:
            print("[Artefacts] Found",len(clusters),"clusters")

        cs_mam.art_rois = []
        for c in clusters:
            npclus = np.array(c)
            xvals = npclus[:,0]
            yvals = npclus[:,1]
            minx = np.min(xvals)
            maxx = np.max(xvals)
            miny = np.min(yvals)
            maxy = np.max(yvals)
            rad = max(maxy-miny,maxx-minx,1.)/2.
            cs_mam.art_rois.append([(maxx+minx+1)/2.,(maxy+miny+1)/2., rad])

            if cs_mam.verbose:
                print("[Artefacts]...size",len(c))

        if cs_mam is not None:
            cs_mam.art_clusters = copy.deepcopy(clusters)
            cs_mam.art_threshold = thresh

        if uiobject:
            uiobject.pbar.label.setText("Finished.")
        error = False
        return error

#----------------------------------------------------------------------
    def DICOMInfo(self,cs,info='dicom'):
        # Different from ImageJ version; tags "0008","0104" and "0054","0220"
        #  appear to be part of sequences. This gives problems (cannot be found
        #  or returning whole sequence blocks)
        # Possibly this can be solved by using if(type(value) == type(dicom.sequence.Sequence()))
        #  but I don't see the relevance of these tags anymore, so set them to NO

        if(info == "dicom"):
            dicomfields = [ ["0008,0021",  "Series Date"],
                            ["0008,0031",  "Series Time"],
                            ["0008,0070",  "Manufacturer"],
                            ["0008,0080",  "InstitutionName"],
                            ["0008,1010",  "StationName"],
                            ["0008,1030",  "StudyDescription"],
                            ["0008,103E",  "SeriesDescription"],
                            ["0008,1070",  "Operator"],
                            ["0010,0020",  "PatientID"],
                            ["0018,0060",  "kVp"],
                            ["0018,1000",  "DeviceSerialNumber"],
                            ["0018,1020",  "SoftwareVersions"],
                            ["0018,1110",  "DistanceSourceToDetector"],
                            ["0018,1111",  "DistanceSourceToPatient"],
                            ["0018,1150",  "ExposureTime_(ms)"],
                            ["0018,1151",  "TubeCurrent_(mA)"],
                            ["0018,1153",  "muAs"],
                            ["0018,1164",  "ImagerPixelSpacing"],
                            ["0018,1166",  "Grid"],
                            ["0018,1190",  "FocalSpot"],
                            ["0018,1191",  "AnodeTargetMaterial"],
                            ["0018,11A0",  "BodyPartThickness"],
                            ["0018,11A2",  "CompressionForce"],
                            ["0018,1405",  "RelativeXRayExposure"],
                            ["0018,700A",  "DetectorID"],
                            ["0018,700C",  "DateOfLastDetectorCalibration"],
                            ["0018,7050",  "FilterMaterialLT"],
                            ["0019,1029",  "---"],
                            ["0028,0101",  "BitsStored"],
                            ["0028,1040",  "PixelensityRelationship"],
                            ["0028,1052",  "Rescaleercept"],
                            ["0028,1053",  "RescaleSlope"],
                            ["0028,1054",  "RescaleType"],
                            ["0040,0302",  "EntranceDose"],
                            ["0040,0314",  "HalfValueLayer_(mm)"],
                            ["0040,0316",  "OrganDose"],
                            ["0040,8302",  "EntranceDose_(mGy)"],
                            ["0000,0000",  "NOViewCodeSequence"],
                            ["0000,0000",  "NOViewCodeMeaning"]]
            if cs.scannername == lit.stL50:
                dicomfields[0] = ["0008,0022",  "Acquisition Date"]
                dicomfields[1] = ["0008,0030",  "Acquisition Time"]

        elif(info == "qc"):
            dicomfields = [["0008,0021",  "Series Date"],
                           ["0008,1010",  "StationName"],
                           ["0008,1070",  "Operator"],
                           ["0018,0060",  "kVp"],
                           ["0018,1020",  "SoftwareVersions"],
                           ["0018,1030",  "ProtocolName"],
                           ["0018,1110",  "DistanceSourceToDetector"],
                           ["0018,1111",  "DistanceSourceToPatient"],
                           ["0018,1153",  "muAs"],
                           ["0018,1166",  "Grid"],
                           ["0018,1190",  "FocalSpot"],
                           ["0018,1191",  "AnodeTargetMaterial"],
                           ["0018,11A0",  "BodyPartThickness"],
                           ["0018,11A2",  "CompressionForce"],
                           ["0018,700A",  "DetectorID"],
                           ["0018,700C",  "DateOfLastDetectorCalibration"],
                           ["0018,7050",  "FilterMaterialLT"],
                           ["0040,0314",  "HalfValueLayer_(mm)"],
                           ["0040,0316",  "OrganDose"],
                           ["0040,8302",  "EntranceDose_(mGy)"]]
            if cs.scannername == lit.stL50:
                dicomfields[0] = ["0008,0022",  "Acquisition Date"]
                dicomfields.append( ["0019,10c1",  "EnergyComponent"] )
        elif info == "id":
            dicomfields = [
                ["0008,0022",  "AcquisitionDate"],
                ["0008,0032",  "AcquisitionTime"],
                ["0018,7050",  "FilterMaterialLT"]
            ]

        results = []
        for df in dicomfields:
            key = df[0]
            value = ""
            try:
#                value =  cs.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
                value = self.readDICOMtag(cs,key)
            except:
                value = ""

            results.append( (df[1],value) )

        return results

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
