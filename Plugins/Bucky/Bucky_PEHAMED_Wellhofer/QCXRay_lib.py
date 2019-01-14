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
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED! And make sure rescaling is corrected!

TODO:
Changelog:
    20180318: fix thresh undefined
    20171116: fix scipy version 1.0
    20161220: Removed class variables; removed testing stuff
    20161216: allow override of auto determine maybeinvert; try other field is pixelspacing missing
    20160902: sync wad2.0 with pywad1.0
    20151112: Bugfix in calculation of lowcontrast
    20151029: Changed output levels and descriptions etc.
    20151027: Added sign to XRayDev
    20150817: Removed hardcoded roomdefinitions. All in config files now!
    20150618: Tried DDL version of fftorientation, but that performs worse here. Left interface intact
    20150610: Bugfix: uniformity was calculated over copper only, instead of non-copper!; moved some math to wadwrapper_lib
    20150522: Force scipy version at least 0.10.1 to avoid problems with PIL/Pillow mixing
    20150521: Again, fix mixing PIL and not PIL
    20150430: allow crash if no PIL.ImageDraw available
    20150428: ImageDraw and Draw both from PIL or both NOT from PIL
    20150331: Fix found positions outside image dimensions; fixed DiDi distance detectot to phantom now 65 (t2w 1630)
    20150330: Made x-ray edge finding bit more robust
    20150205: DICOMInfo according to general interface
    20141017: Room CALHOS1 added; determine 90,180,270 rotation from std now instead of avg
    20141013: Fix WKZ threshold TableOrWall (string was compared to int)
    20140815: Fix xRayedges where found threshhigh<threslow
    20140717: Fix expertbox QC, adjusted WKZ1 S Threshold
    20140626: Moved MaybeInvert and DetermineDeviceID into init of XRayStruct
    20140625: Auto unrotate Wellhofer; sorting MTF; simplified most mentions of mustinvert;
              fix for po_box on annotation; fix for invert LowContrast
    20140623: First attempt to rewrite into WAD module; speedup and bugfix of Uniformity()
"""
__version__ = '20180318'
__author__ = 'aschilham'

try:
    import pydicom as dicom
except ImportError:
    import dicom
import numpy as np
import scipy.ndimage as scind
    
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

import operator
import matplotlib.pyplot as plt
import copy
from numpy import fft as fftpack
from scipy import stats
import numpy.ma as ma
from PIL import Image # image from pillow is needed
from PIL import ImageDraw # imagedraw from pillow is needed, not pil
import scipy.misc
# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
if scipy_version[0] == 0:
    if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
        raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")

try:
    # wad2.0 runs each module stand alone
    import QCXRay_constants as lit
    import QCXRay_math as mymath
except ImportError:
    from . import QCXRay_constants as lit
    from . import QCXRay_math as mymath

class Room :
    def __init__ (self,_name, outvalue=-1, tablesid=-1, wallsid=-1, tablepid=-1, wallpid=-1,
                  phantom=lit.stPehamed,sdthresh=-1,sens_threshold = [],
                  linepairmarkers = {}, mustbeinverted=None):
        self.name = _name # identifier of room
        
        self.outvalue = outvalue # value of pixels outside x-ray field
        
        self.pidtablemm = tablepid # distance between mid phantom and detector in mm
        self.pidwallmm  = wallpid
        self.sidtablemm = tablesid
        self.sidwallmm = wallsid
        self.phantom = phantom

        self.sdthresh = sdthresh             # threshold on stdev for determin table or wall, now only for CALHOS, maybe WKZ?
        self.sens_threshold = sens_threshold # list of (date,sensitivity) describing from dates and max threshold on sensitivity for table
        self.mustbeinverted = mustbeinverted # allow force

        self.skipFFT = False # only for wellhofer
        if(phantom == lit.stWellhofer):
            self.skipFFT = True

        if len(linepairmarkers)>0:
            self.xy06mm = linepairmarkers['xymm0.6'] 
            self.xy14mm = linepairmarkers['xymm1.4']
            self.xy18mm = linepairmarkers['xymm1.8']
            self.xy46mm = linepairmarkers['xymm4.6']
        elif phantom == lit.stWellhofer:
            # wellhofer WKZ
            self.xy18mm = [53.7, 27.1] # x,y position in mm of decimal dot in 1.8 lp/mm 
            self.xy06mm = [80.6,-04.7] # x,y position in mm of decimal dot in 0.6 lp/mm 
            self.xy14mm = [59.9,-24.9] # x,y position in mm of decimal dot in 1.4 lp/mm 
            self.xy46mm = [28.4, 01.8] # x,y position in mm of decimal dot in 4.6 lp/mm 
        else: # pehamed phantom nr 100876
            self.xy18mm = [53.1, 27.1] # x,y position in mm of decimal dot in 1.8 lp/mm 
            self.xy06mm = [80.1,-04.5] # x,y position in mm of decimal dot in 0.6 lp/mm 
            self.xy14mm = [59.4,-24.7] # x,y position in mm of decimal dot in 1.4 lp/mm 
            self.xy46mm = [28.0, 02.0] # x,y position in mm of decimal dot in 4.6 lp/mm 
        
    def setPIDs(self,_pidtable, _pidwall):
        self.pidtablemm = _pidtable
        self.pidwallmm = _pidwall
    def setSIDS(self,_sidtable, _sidwall):
        self.sidtablemm = _sidtable
        self.sidwallmm = _sidwall

class XRayStruct:
    roomUnknown = Room(lit.stUnknown)

    class UnifStruct :
        def __init__ (self):
            self.ROIuniformity = self.LRuniformity = 0. # fraction
            self.BKmean = self.BKsdev = 0. # gridless bk
            self.peakValue = 0. # max valule over uniformity line
            self.posval = []
            self.intens = []
            self.trendMax = self.trendMin = 0.
            self.trend = []
            self.roi = []
            self.verticalROI = False

    class CuStruct :
        def __init__ (self):
            self.roi_snr = []
            self.roi_cnr = []
            self.roi_mmcu = []
            self.roi_mean = []
            self.roi_sdev = []
            self.dynamicRange = -1
            self.guesskVp     = -1
            self.slope        = 0.
            self.roi = []
            self.step_rois = []
            self.wedge_confidence = -1.

    class LoCoStruct :
        def __init__ (self):
            self.low_cnr = []
            self.mean_sg = []
            self.mean_bk = []
            self.sdev_sg = []
            self.sdev_bk = []
            self.lo_rois = []

    class MTFStruct:
        def __init__ (self):
            self.dotxys = []
            self.contrast_response = []
            self.contrast_high = []
            self.contrast_low = []
            self.ctfmtf = []
            self.contrast_freqs = []
            self.contrast_tops = []
            self.contrast_bots = []
            self.calc_freqs = []  # frequencies in phantom units as measured from extremes
            self.mtf_aapm = -1.
            self.freq_confidence = -1.
            self.pos_confidence = -1.
            self.roi = []

    def maybeInvert(self):
        if not self.mustbeinverted is None:
            print("Must be Inverted (use)",self.mustbeinverted)
            return
        
        if self.dcmInfile is None:
            return
        self.mustbeinverted = False
        if self.dcmInfile.PhotometricInterpretation == "MONOCHROME2":
            self.mustbeinverted = True
        print("Must be Inverted",self.mustbeinverted)


    def __init__ (self,dcmInfile,pixeldataIn,room):
        self.verbose = False

        # input image
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.knownTableOrWall = None

        self.expertOverridepixToGridScaleCm = None
        self.forceRoom = room
        self.lastimage = None # GUI feedback

        # phantom orientation
        self.po_center_roi = [] # xmid,ymid,rad
        self.po_roi = [] # list [ [x,y] ]
        self.po_rot = 0
        self.po_fftroi = []
        self.test_rois = []
        self.bbox_confidence = 0.

        # xray field
        self.xrayNSWEmm = []
        self.xr_roi = []

        # horz. uniformity
        self.unif = self.UnifStruct()

        # Cu Wedge
        self.cuwedge = self.CuStruct()

        # MTF
        self.mtf = self.MTFStruct()

        # Low Contrast
        self.loco = self.LoCoStruct()

        self.mustbeinverted = None # if pixval(Cu) = high and pixval(air) = low, then mustbeinverted
        if room.mustbeinverted is None:
            self.maybeInvert()
        else:
            self.mustbeinverted = room.mustbeinverted

        # for matlib plotting
        self.hasmadeplots = False

class XRayQC:
    def __init__(self):
        self.qcversion = __version__
        self.boxradmm   = 110  # choose 11 cm or 8 cm for clean surroundings
        self.adjustmtfangledeg = 0. # if consistency check fails, add a little angle
        self.bShowMTFDetail = False
        self.bShowCTF = False
        self.bIgnoreMTFError = False # in rare cases MTF is "ok" but looks wrong
    
        self.crLimit = 0.1 # contrast limit for MTF
        self.sigma_ext = 1.5 # gaussian blur factor for extreme finder

    def pixToGridScaleCm(self,cs):
        if not cs.expertOverridepixToGridScaleCm is None:
            return cs.expertOverridepixToGridScaleCm

        try:
            pixel_spacing_x = cs.dcmInfile.PixelSpacing[0]
        except AttributeError:
            pixel_spacing_x = cs.dcmInfile.ImagerPixelSpacing[0]
            
        stand = self.TableOrWall(cs)
        # calculate distance = L0-Thick
        if not 'DistanceSourceToDetector' in cs.dcmInfile:
            if stand == lit.stTable:
                sidmm = cs.forceRoom.sidtablemm
                pidmm = cs.forceRoom.pidtablemm
            else:
                sidmm = cs.forceRoom.sidwallmm
                pidmm = cs.forceRoom.pidwallmm
        else:
            distance = ""
            try:
                distance = cs.dcmInfile.DistanceSourceToDetector
            except:
                distance = ""

            sidmm = -1.
            pidmm = 100  # distance phantom to detector
            if(distance != ""):
                sidmm = distance
                if(stand == lit.stTable):
                    pidmm = cs.forceRoom.pidtablemm # DiDi: 8cm +half phantom thickness
                else:
                    pidmm = cs.forceRoom.pidwallmm # DiDi: 8cm +half phantom thickness
        return  pixel_spacing_x*(sidmm-pidmm)/sidmm

    def pix2phantomm(self, cs, pix):
        pix2phantommm = self.pixToGridScaleCm(cs)
        return pix*pix2phantommm

    def phantommm2pix(self, cs,mm):
        pix2phantommm = self.pixToGridScaleCm(cs)
        return mm/pix2phantommm

    def diamondNESW(self,roipts_orig):
        # turn into diamond
        roipts = []
        #North
        ix1 = 0
        ix2 = 3
        roipts.append([(int)(0.5+0.5*(roipts_orig[ix1][0]+roipts_orig[ix2][0])),(int)(0.5+0.5*(roipts_orig[ix1][1]+roipts_orig[ix2][1]))])
        #East
        ix1 = 2
        ix2 = 3
        roipts.append([(int)(0.5+0.5*(roipts_orig[ix1][0]+roipts_orig[ix2][0])),(int)(0.5+0.5*(roipts_orig[ix1][1]+roipts_orig[ix2][1]))])
        #South
        ix1 = 1
        ix2 = 2
        roipts.append([(int)(0.5+0.5*(roipts_orig[ix1][0]+roipts_orig[ix2][0])),(int)(0.5+0.5*(roipts_orig[ix1][1]+roipts_orig[ix2][1]))])
        #West
        ix1 = 0
        ix2 = 1
        roipts.append([(int)(0.5+0.5*(roipts_orig[ix1][0]+roipts_orig[ix2][0])),(int)(0.5+0.5*(roipts_orig[ix1][1]+roipts_orig[ix2][1]))])
        return roipts

    def phantomposmm2pix(self, roipts_orig, xmm, ymm):
        """
        Convert position as indicated on phantom to image pixels

        input: x and y coordinates in phantom grid mm
        output: x and y in pix
        """

        dirNESW = self.diamondNESW(roipts_orig)

        xmid = 0.
        ymid = 0.
        for rp in roipts_orig:
            xmid += rp[0]
            ymid += rp[1]
        xmid /= len(roipts_orig)
        ymid /= len(roipts_orig)

        xpx = xmid + ymm/self.boxradmm*(dirNESW[0][0]-xmid) +xmm/self.boxradmm*(dirNESW[1][0]-xmid)
        ypx = ymid +ymm/self.boxradmm*(dirNESW[0][1]-ymid)+xmm/self.boxradmm*(dirNESW[1][1]-ymid)

        return xpx,ypx

#----------------------------------------------------------------------
    def invertmaxval(self,cs):
        # Not stored properly in self.dcmfileIn.BitsStored
        dicomfields = [ ["0028,0101",  "Bits Stored"]]
        key = dicomfields[0][0]
        dicvalue = cs.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
        if dicvalue != "":
            return dicvalue
        else:
            return 0

    def TableOrWallFromSTDEV(self,cs):
        """
        Voor CALHOS (fosfor) kan Wand/Tafel onderscheid gemaakt worden met std van groot gedeelte image
        Misschien ook voor andere fosfor (WKZ)?
        :param cs:
        :return:
        """
        # 1. Maak box van x = 100/350 tot 250/350, en dan vierkant voor y
        # 2. Bereken avg en stdev; voor tafel typisch 50, voor wand typisch 30->grens op 40

        # 1.1 find approx center of phantom (screw)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        rad = min(widthpx/2,heightpx/2)
        immidx = int(.5*(widthpx-1)+.5)
        immidy = int(.5*(heightpx-1)+.5)
        xlo = int(immidx-rad/2+.5)
        xhi = int(immidx+rad/2+.5)
        ylo = int(immidy-rad/2+.5)
        yhi = int(immidy+rad/2+.5)
        roipts = [
            [xlo,ylo],
            [xlo,yhi],
            [xhi,yhi],
            [xhi,ylo]
        ]
        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1
        if cs.mustbeinverted:
            smallimage = invertmax-cs.pixeldataIn[xlo:xhi,ylo:yhi]
        else:
            smallimage = cs.pixeldataIn[xlo:xhi,ylo:yhi]

        avg = np.mean(smallimage)
        std = np.std(smallimage)
        if std<cs.forceRoom.sdthresh:
            return lit.stWall
        else:
            return lit.stTable

    def TableOrWall(self, cs):
        if cs.knownTableOrWall is not None:
            return cs.knownTableOrWall

        result = ""
        if cs.forceRoom.sidtablemm < 0. and cs.forceRoom.sidwallmm > 0.:
            result =  lit.stWall
        elif cs.forceRoom.sidtablemm > 0. and cs.forceRoom.sidwallmm < 0.:
            result =  lit.stTable
        elif cs.forceRoom.sdthresh>0.:
            result =  self.TableOrWallFromSTDEV(cs)
        elif len(cs.forceRoom.sens_threshold) > 0:
            # check on sensitivity
            sens = cs.dcmInfile.Sensitivity
            cdate = int(cs.dcmInfile.SeriesDate)
            thresh = cs.forceRoom.sens_threshold[-1][1]
            for maxdate,maxthresh in cs.forceRoom.sens_threshold:
                if cdate < maxdate:
                    thresh = maxthresh
                    break

            if sens<thresh:
                result = lit.stTable
            else:
                result = lit.stWall
        else:
            try:
                sid = cs.dcmInfile.DistanceSourceToDetector
            except:
                return lit.stUnknown

            if sid>1600.:
                result = lit.stWall
            else:
                result = lit.stTable
#        print("[TableOrWall] %s %d %d (%d): %s" %(cs.forceRoom.name,int(cs.dcmInfile.SeriesDate),sens,thresh,result))
        if result != lit.stUnknown:
            cs.knownTableOrWall = result
        else:
            print("ERROR! Cannot determine if Wall or Table is used!")
        return result

#----------------------------------------------------------------------
    def findPhantomOrientation(self,cs,ebbox=None):
        error = True

        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        # 1.1 find approx center of phantom (screw)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        sqsize = min(widthpx,heightpx)
        midx = int(.5*(widthpx-1)+.5)
        midy = int(.5*(heightpx-1)+.5)
        sqsize3 = int(sqsize/3)
        if cs.mustbeinverted:
            smallimage = invertmax - wadwrapper_lib.extract(cs.pixeldataIn, [sqsize3,sqsize3],[midx,midy])
        else:
            smallimage = wadwrapper_lib.extract(cs.pixeldataIn, [sqsize3,sqsize3],[midx,midy])

        smallimage = scind.gaussian_filter(smallimage, sigma=5)
        cs.lastimage = smallimage

        x0,y0 = np.unravel_index(smallimage.argmin(), smallimage.shape)
        immidx = midx-smallimage.shape[0]/2+x0
        immidy = midy-smallimage.shape[1]/2+y0
        rad=20
        cs.po_center_roi = [immidx,immidy,rad] # feedback for GUI

        # 1.2 position 10x10phantomcm box and manipulate corner positions
        pix2phantommm = self.pixToGridScaleCm(cs)
        # choose 11cm each direction or 8 cm
        if cs.forceRoom.phantom == lit.stWellhofer: 
            self.boxradmm = 80
        rad = 2*self.boxradmm/pix2phantommm
        roipts = [
            [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
            [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
            [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
            [int(immidx+rad/2+.5),int(immidy-rad/2+.5)]
        ]
        if ebbox is not None:
            roipts = ebbox
            radpts = np.sqrt( (roipts[1][0]-roipts[0][0])**2. + (roipts[1][1]-roipts[0][1])**2. )
            mindist = radpts
            for i in [110,80,90,70,60]:
                radtest = 2*i/pix2phantommm
                dist = np.abs(radtest-radpts)
                if dist < mindist:
                    mindist = dist
                    self.boxradmm = i
            print("found mindist=",mindist,"rad=",self.boxradmm)
        # try to find rotation of phantom
        if not cs.forceRoom.skipFFT:
            error,roipts,rotangledeg = self.FieldRotationFFT(cs,roipts)
            if error: # probably too small angle, so try again
                error,roipts,rotangledeg = self.FieldRotationFFT(cs,roipts,10.)

        else:
            rotangledeg = 0.
            error = False

        cs.bbox_confidence = 0.
        if not error:
            last_attempt= False
            radmmtests = [110,80,90,70,60]
            if cs.forceRoom.phantom == lit.stWellhofer: 
                radmmtests = [80,90,70,60,110]
            for radmm in radmmtests:
                self.boxradmm = radmm
                rad = 2*self.boxradmm/pix2phantommm
                roipts = [
                    [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
                    [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
                    [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
                    [int(immidx+rad/2+.5),int(immidy-rad/2+.5)] ]
                roipts = self.RotateBoundingBox(roipts,rotangledeg)
                error,cs.bbox_confidence = self.AlignROI(cs,roipts,"BoundingBox",last_attempt)
                for (xa,ya) in roipts:
                    if xa<0 or ya<0:
                        error = True
                    elif xa>=widthpx or ya>=heightpx:
                        error = True
                if error == True and np.abs(rotangledeg)<2.:
                    roipts = [
                        [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
                        [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
                        [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
                        [int(immidx+rad/2+.5),int(immidy-rad/2+.5)] ]
                    if radmm<61:
                        last_attempt = True
                    error,cs.bbox_confidence = self.AlignROI(cs,roipts,"BoundingBox",last_attempt)
                if not error:
                    break

        print("using rad,conf:",self.boxradmm,cs.bbox_confidence)
        cs.po_roi = roipts
        return error

    def RotateBoundingBox(self, roipts_orig,rotangledeg):
        xmid = 0.
        ymid = 0.
        for rp in roipts_orig:
            xmid += rp[0]
            ymid += rp[1]
        xmid /= len(roipts_orig)
        ymid /= len(roipts_orig)
        xmid  = int(xmid+.5)
        ymid  = int(ymid+.5)

        anglerad = (rotangledeg/180.*np.pi)
        costerm = np.cos(-anglerad)
        sinterm = np.sin(-anglerad)
        roipts = []
        for rp in roipts_orig:
            xp = xmid +(rp[0]-xmid)*costerm-(rp[1]-ymid)*sinterm
            yp = ymid +(rp[0]-xmid)*sinterm+(rp[1]-ymid)*costerm
            roipts.append([xp,yp])
        return roipts

    def _fieldRotationFFT(self,cs,smallimage,initangle=None):
        # also tried runmode==2, which is the DDL version (with slight mods), but this old one works far better
        runmode = 1
        # Start FFT
        F1 = fftpack.fft2(smallimage-np.mean(smallimage))
        # Now shift the quadrants around so that low spatial frequencies are in
        # the center of the 2D fourier transformed image.
        F2 = fftpack.fftshift( F1 )
        if runmode == 1:
            # Calculate a 2D power spectrum
            psd2D = np.abs( F2 )**2.
            fwid = psd2D.shape[0]
            fhei = psd2D.shape[1]
            cs.lastimage = psd2D
    
            fwid2 = int(fwid/2)
            fhei2 = int(fhei/2)
            kwart = np.zeros((fwid2,fhei2),dtype=float)
            for x in range(0,fwid2):
                for y in range (0,fhei2):
                    kwart[x,y] += psd2D[fwid2+x,fhei2+y]
                    kwart[x,y] += psd2D[fhei2-y,fwid2+x]
    
            kwartmax = np.max(kwart)
            # Find local maxima
            posxmax = []
            posymax = []
            kwart2 = np.zeros((fwid2,fhei2),dtype=float)
            while len(posxmax)<20 and np.max(kwart)>0.05:
                xx,yy = np.unravel_index(kwart.argmax(), kwart.shape)
                if( not(xx<3 and yy<3) ):
                    posxmax.append(xx)
                    posymax.append(yy)
                    kwart2[xx,yy] = kwart[xx,yy]
                kwart[xx,yy] = 0
    
            # a bit of sorting:
            index = list(range(len(posxmax)))
            index.sort(key = posxmax.__getitem__)
            posxmax[:] = [posxmax[i] for i in index]
            posymax[:] = [posymax[i] for i in index]
            #for x,y in zip(posxmax,posymax):
            #    print(x,y)
    
            maxoff = len(posxmax)-5
            r2angleoff = []
            r2intslope = []
            for off in range(0,maxoff):
                # Following kludge is needed to prevent division by zero in linregress for small angle
                nonidentical = True
                ssxm, ssxym, ssyxm, ssym = np.cov(posxmax[off:len(posxmax)],posymax[off:len(posxmax)], bias=1).flat
                if ssxm == 0:
                    nonidentical = False
                if nonidentical:
                    slope, intercept, r1_value, p_value, std_err = stats.linregress(posxmax[off:len(posxmax)],posymax[off:len(posxmax)])
                    anglerad = np.arctan2(1.,slope) # needs inverse angle
                    r2 = r1_value**2
                    r2angleoff.append( (r2,anglerad,off) )
                    r2intslope.append( (r2,intercept,slope,off) )
                #print(r2,anglerad,off)
                nonidentical = True
                ssxm, ssxym, ssyxm, ssym = np.cov(posxmax[0:len(posxmax)-off],posymax[0:len(posxmax)-off], bias=1).flat
                if ssxm == 0:
                    nonidentical = False
                    anglerad = np.pi*2.
                    r1_value = 0.
                    intercept = 0.
                    slope = 0.
                if nonidentical:
                    slope, intercept, r1_value, p_value, std_err = stats.linregress(posxmax[0:len(posxmax)-off],posymax[0:len(posxmax)-off])
                    anglerad = np.arctan2(1.,slope) # needs inverse angle
                r2 = r1_value**2
                r2angleoff.append( (r2,anglerad,-off) )
                r2intslope.append( (r2,intercept,slope,-off) )
                    
                #print(r2,anglerad,-off)
    
            r2angleoff = sorted(r2angleoff)
            r2,anglerad,off = r2angleoff[-1]
    
            if cs.verbose:
                dummy1,intercept,slope ,dummy2 = sorted(r2intslope)[-1]
                plt.figure()
                plt.plot(posxmax,posymax,'bo')
                dafit = [intercept+slope*x for x in posxmax]
                plt.plot(posxmax,dafit,'b-')
                plt.title("Orientation")
                print('orientation fit: %f + %f*x; offset=%d'%(intercept,slope,dummy2))
                print('angledeg=%f,initangledeg=%f;'%(anglerad/np.pi*180.,0 if initangle is None else initangle))
                cs.hasmadeplots = True
        else: # 'new'
            # must yield anglerad and r2
            #Let's log-scale the image so we're dealing with something in uint8 range.
            # Calculate a 2D power spectrum
            abs_data = np.abs( F2 )**2.
            cs.lastimage = abs_data

            # Threshold the lower 25% of the peak
            max_peak = np.max(abs_data)
            abs_data[abs_data < (max_peak * 0.25)] = 0
            
            # Log-scale the data
            abs_data += 1
            c = 255.0 / np.log(1 + max_peak)
            log_data = c * np.log(abs_data)
            
            # fit data through points within 90% of the max peak of the scaled image
            # run for both angles (horz grid and vert grid)
            mask = np.zeros(np.shape(log_data),dtype=bool)
            wid,hei = np.shape(log_data)
            wid2 = int(wid/2)
            hei2 = int(hei/2)
            mask[0:wid2,0:hei2] = True
            mask[wid2:wid,hei2:hei] = True
            mask = ~mask
            
            # we will sort on number of points, as it is easier to fit to less points, but more points should be prefered
            r2anglenum = []
            for i in [0,1]:
                log_data1 = log_data.copy()
                log_data1[mask] = 0
                max_scaled_peak = np.max(log_data1)

                # Determine the angle of two high-peak points in the image
                rows, cols = np.where(log_data1 > (max_scaled_peak * 0.90))
                # Following kludge is needed to prevent division by zero in linregress for small angle
                nonidentical = True
                ssxm, ssxym, ssyxm, ssym = np.cov(rows,cols, bias=1).flat
                if ssxm == 0:
                    nonidentical = False
                if nonidentical:
                    slope, intercept, r1_value, p_value, std_err = stats.linregress(rows,cols)
                    anglerad = np.arctan2(1.,slope) -i*np.pi/2.# needs inverse angle
                else:
                    anglerad = np.pi*2.
                    r1_value = 0.
                    intercept = 0.
                    slope = 0.
                print('run %d angle: %f (%f), r2=%f, n=%d'%(i,anglerad,anglerad/np.pi*180.,r1_value**2.,len(rows)))
                r2anglenum.append( (r1_value**2.,anglerad,len(rows)))
                mask = ~mask
                if cs.verbose:
                    plt.figure()
                    plt.plot(rows,cols,'bo')
                    dafit = [intercept+slope*x for x in rows]
                    plt.plot(rows,dafit,'b-')
                    plt.title("Orientation")
                    print('orientation fit: %f + %f*x'%(intercept,slope))
                    print('angledeg=%f,initangledeg=%f;'%(anglerad/np.pi*180.,0 if initangle is None else initangle))
                    #plt.show()
                    cs.hasmadeplots = True
                
            avg_angle = (r2anglenum[0][1]+r2anglenum[1][1])/2.
            wavg_angle = (r2anglenum[0][1]*r2anglenum[0][2]+r2anglenum[1][1]*r2anglenum[1][2])/(r2anglenum[0][2]+r2anglenum[1][2])
            avg_r2 = np.sqrt( 0.5*(r2anglenum[0][0]**2. +r2anglenum[1][0]**2.) )
            print("avgangle: %f"%(avg_angle/np.pi*180.+(0 if initangle is None else initangle)))
            r2anglenum.sort(key=operator.itemgetter(2)) # 0 = maxr2, 2 = maxnum
            #r2anglenum = sorted(r2anglenum)
            r2,anglerad,num = r2anglenum[-1]
            off = 0
            #r2,anglerad = avg_r2,avg_angle
            #r2,anglerad = avg_r2,wavg_angle
        return anglerad,r2,off # must return angle in rad, confidence measure and offset if any


    def FieldRotationFFT(self,cs,roipts_orig,initangle=None):
        """
        Cut out a part of the phantom that should contain only grid
        Make a FFT, select one quadrant and find the locations of the first N maximums
        Fit a line through the maximums
        """
        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        error = True
        xmid = 0.
        ymid = 0.
        for rp in roipts_orig:
            xmid += rp[0]
            ymid += rp[1]
        xmid /= len(roipts_orig)
        ymid /= len(roipts_orig)
        xmid  = int(xmid+.5)
        ymid  = int(ymid+.5)

        boxradpx = int(self.phantommm2pix(cs,20.))
        yoffpx = int(self.phantommm2pix(cs,10.))
        ll = [xmid-boxradpx,ymid-yoffpx]
        ur = [xmid+boxradpx,ymid-yoffpx-2*boxradpx]
        cs.po_fftroi = [ ll, [ur[0],ll[1]], ur, [ll[0],ur[1]] ] # x0,wid y0,height

        xlo = min(ll[0],ur[0])
        ylo = min(ll[1],ur[1])
        xhi = max(ll[0],ur[0])
        yhi = max(ll[1],ur[1])
        if not (initangle is None): # allow for extra artificial rotation of phantom
            smallimage = scind.interpolation.rotate(cs.pixeldataIn,initangle, axes=(1, 0), reshape=True, output=None, order=3, mode='constant', cval=0.0, prefilter=True)[xlo:xhi+1,ylo:yhi+1]
        else:
            smallimage = cs.pixeldataIn[xlo:xhi+1,ylo:yhi+1]
        if cs.mustbeinverted:
            smallimage = invertmax - smallimage

        cs.lastimage = smallimage
        #return True, roipts_orig,-360

        anglerad,confidence,off = self._fieldRotationFFT(cs, smallimage,initangle)

        if initangle != None:
            anglerad += initangle/180*np.pi
        if anglerad > np.pi/4.:
            anglerad = anglerad - np.pi/2.
        rotangledeg = (anglerad/np.pi*180.)

        if confidence<0.85:
            label = "first try"
            if initangle!=None:
                label = "Error!"
            print("FieldRotationFFT:",label,"confidence too low:",confidence,off)
            #print(offanglerad)
            return error,roipts_orig,rotangledeg
        error = False

        print("rotangledegFFT:",rotangledeg,confidence,off)

        roipts = self.RotateBoundingBox(roipts_orig,rotangledeg)
        return error,roipts,rotangledeg

    def AlignROI(self, cs, roipts, what,blast_attempt=True):
        error = True
        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        pix2phantommm = self.pixToGridScaleCm(cs)
        # choose 11cm each direction or 8 cm
        searchrad = int(2.5/pix2phantommm+.5) # 2 mm around location

        workimage = cs.pixeldataIn
        if what == "MTF": #
            searchrad = int(0.5/pix2phantommm+.5)
#            workimage = self.templateMatchDisc(cs)
#            cs.lastimage = workimage

        searchrad = max(1,searchrad)
        print("%s searchrad="%what,searchrad)
        widthpx = np.shape(workimage)[0] ## width/height in pixels
        heightpx = np.shape(workimage)[1]

        conf_pts = []
#        conf_pts.append((self.ROIConfidence(cs,roipts,what),copy.deepcopy(roipts)))
        conf_pts.append((0,copy.deepcopy(roipts))) # never not select default pointset
        roipts0 = copy.deepcopy(roipts)
        for kk in range(0,6):
            for i in range(0, len(roipts)):
                rp = roipts[i]
                x0 = rp[0]
                minx = int(max(0,x0-searchrad))
                maxx = int(min(widthpx-2,x0+searchrad))
                y0 = rp[1]
                miny = int(max(0,y0-searchrad))
                maxy = int(min(heightpx-2,y0+searchrad))
                if cs.mustbeinverted:
                    cropped = workimage[minx:maxx+1,miny:maxy+1]
                else:
                    cropped = invertmax - workimage[minx:maxx+1,miny:maxy+1]
#                plt.figure()
#                plt.imshow(cropped)
#                plt.show()
                # smoothing to get rid of noise
                sigma = max(1.5,cropped.shape[0]/8.)
                cropped = scind.gaussian_filter(cropped, sigma)
                # find maximum location in projections
                xflat = np.mean(cropped,axis=0)
                yflat = np.mean(cropped,axis=1)
                x1 = minx+np.unravel_index(yflat.argmax(), yflat.shape)[0]
                y1 = miny+np.unravel_index(xflat.argmax(), xflat.shape)[0]
                if cs.verbose:
                    if what == "MTF" and i==2:
                        plt.figure()
                        plt.imshow(cropped)
                        plt.title('Align '+str(kk)+ ' '+what+str(i))
                        print(sigma,"shift ",i," of point ",kk," =(",x1-x0,",",y1-y0,")")

                rp[0] = x1
                rp[1] = y1
            roipts = self.ConsistencyAlign(cs,roipts0,roipts,what)
            conf_pts.append((self.ROIConfidence(cs,roipts,what),copy.deepcopy(roipts)))

        # Just take the last one; should be best!
#        conf_pts = sorted(conf_pts) #sorting orientation/distance box does worsen results
        if what == "MTF": # sorting MTF seems to help
#            print("-------------")
#            for i,copts in enumerate(conf_pts):
#                print(i,copts[0])
#            print("-------------")
            conf_pts = sorted(conf_pts)

        for i in range(0,len(roipts)):
            roipts[i][0] = conf_pts[-1][1][i][0]
            roipts[i][1] = conf_pts[-1][1][i][1]
        confidence = conf_pts[-1][0]
        if confidence> 0.8:
            error = False
        else:
            label = ""
            if blast_attempt==True:
                label = "Error!"
            print("AlignRoi (",what,"):",label,", confidence too low:",confidence)

        return error,confidence

    def ConsistencyAlign(self, cs,roiptsold, roiptsnew,what):
        """
        4 types of motion:
            1. Shift (all points move in same direction)
            2. Scale (all distance to center change by same amount)
            3. Rotation (all points rotate along same angle around center)
            4. Shear (some points rotate more than others. Do not allow)
        Some heuristics to correct for one corner point being very obviously wrong
        """
        midx0 = 0.
        midy0 = 0.
        for rp in roiptsold:
            midx0 += rp[0]
            midy0 += rp[1]
        midx0 /= len(roiptsold)
        midy0 /= len(roiptsold)

#        shifts = []
        shiftdist = []
#        scales = []
#        angles = []
        for rp0,rp1 in zip(roiptsold, roiptsnew):
#            shifts.append([rp1[0]-rp0[0],rp1[1]-rp0[1]])
            shiftdist.append(np.sqrt((rp1[0]-rp0[0])**2+(rp1[1]-rp0[1])**2))
#            dist1 = np.sqrt((rp1[0]-midx0)**2+(rp1[1]-midy0)**2)
#            scales.append(dist1/d0)
#            angles.append(np.arctan2(rp1[1]-rp0[1],rp1[0]-rp0[0]))

        if(what == "BoundingBox"):
            """
            [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
             [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
              [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
               [int(immidx+rad/2+.5),int(immidy-rad/2+.5)] ]
            """
            idcmp = [ [[0,1],[3,2]], [[3,2],[0,1]], [[0,3],[1,2]], [[1,2],[0,3]] ]
            for abc in idcmp:
                i0 = abc[0][0]
                i1 = abc[0][1]
                i2 = abc[1][0]
                i3 = abc[1][1]
                l1 = np.sqrt( (roiptsnew[i1][0]-roiptsnew[i0][0])**2+ (roiptsnew[i1][1]-roiptsnew[i0][1])**2)
                if(np.abs(self.pix2phantomm(cs,l1)-self.boxradmm*2)>5.):
                    if(shiftdist[i0]>shiftdist[i1]):
                        roiptsnew[i0][0] = roiptsnew[i1][0]-(roiptsnew[i3][0]-roiptsnew[i2][0])
                        roiptsnew[i0][1] = roiptsnew[i1][1]-(roiptsnew[i3][1]-roiptsnew[i2][1])
                    else:
                        roiptsnew[i1][0] = roiptsnew[i0][0]+(roiptsnew[i3][0]-roiptsnew[i2][0])
                        roiptsnew[i1][1] = roiptsnew[i0][1]+(roiptsnew[i3][1]-roiptsnew[i2][1])

        if(what == "MTF"):
            """
            roipts = [ [x18px,y18px],[x06px,y06px],[x14px,y14px],[x46px,y46px] ]
            """
            self.adjustmtfangledeg = 0.
            id18 = 0
            id06 = 1
            id14 = 2
            id46 = 3
            idcmp = [ [id18,id06], [id06,id14], [id14,id46], [id46,id18] ]
            dd18_46 = [roiptsnew[id46][0]-roiptsnew[id18][0],roiptsnew[id46][1]-roiptsnew[id18][1]]
            dd06_14 = [roiptsnew[id14][0]-roiptsnew[id06][0],roiptsnew[id14][1]-roiptsnew[id06][1]]
            lengthmms = [41.4298705358224,29.1192821547155,41.1481478200461,35.5992147241798 ]
            len18_46 = lengthmms[3]
            len06_14 = lengthmms[1]
            for abc,l0mm in zip(idcmp,lengthmms):
                i0 = abc[0]
                i1 = abc[1]
                l1 = np.sqrt( (roiptsnew[i1][0]-roiptsnew[i0][0])**2+ (roiptsnew[i1][1]-roiptsnew[i0][1])**2)
                if(np.abs(self.pix2phantomm(cs,l1)-l0mm)>.45):
                    alter = i1
                    if(shiftdist[i0]>shiftdist[i1]):
                        alter = i0
                    if (alter == id18): #18:
                        roiptsnew[alter][0] = int(0.5+roiptsnew[id46][0]-len18_46/len06_14*dd06_14[0])
                        roiptsnew[alter][1] = int(0.5+roiptsnew[id46][1]-len18_46/len06_14*dd06_14[1])
                        self.adjustmtfangledeg = 0.5
                    elif(alter == id06): #06
                        roiptsnew[alter][0] = int(0.5+roiptsnew[id14][0]-len06_14/len18_46*dd18_46[0])
                        roiptsnew[alter][1] = int(0.5+roiptsnew[id14][1]-len06_14/len18_46*dd18_46[1])
                    elif(alter == id14): #14
                        roiptsnew[alter][0] = int(0.5+roiptsnew[id06][0]+len06_14/len18_46*dd18_46[0])
                        roiptsnew[alter][1] = int(0.5+roiptsnew[id06][1]+len06_14/len18_46*dd18_46[1])
                    elif(alter == id46): #46
                        roiptsnew[alter][0] = int(0.5+roiptsnew[id18][0]+len18_46/len06_14*dd06_14[0])
                        roiptsnew[alter][1] = int(0.5+roiptsnew[id18][1]+len18_46/len06_14*dd06_14[1])
                        self.adjustmtfangledeg = 0.5
            for abc,l0mm in zip(idcmp,lengthmms):
                i0 = abc[0]
                i1 = abc[1]
                l1 = np.sqrt( (roiptsnew[i1][0]-roiptsnew[i0][0])**2+ (roiptsnew[i1][1]-roiptsnew[i0][1])**2)

        #sanity check
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        for (x0,y0) in roiptsnew:
            if x0 <0 or y0<0:
                return roiptsold
            if x0>=widthpx or y0>=heightpx:
                return roiptsold
        return roiptsnew

    def ROIConfidence(self,cs,roipts,what):
        confidence = 0.

        # First the lengths of just the sides of the box
        pix2phantommm = self.pixToGridScaleCm(cs)
        redlength = self.boxradmm/pix2phantommm
        lengths = [ ]
        for i in range(0,3):
            dx = roipts[i][0]-roipts[i+1][0]
            dy = roipts[i][1]-roipts[i+1][1]
            lengths.append(np.sqrt(dx*dx+dy*dy))
        dx = roipts[3][0]-roipts[0][0]
        dy = roipts[3][1]-roipts[0][1]
        lengths.append(np.sqrt(dx*dx+dy*dy))
        # now add diagonals
        dx = roipts[2][0]-roipts[0][0]
        dy = roipts[2][1]-roipts[0][1]
        lengths.append(np.sqrt(dx*dx+dy*dy))
        dx = roipts[1][0]-roipts[3][0]
        dy = roipts[1][1]-roipts[3][1]
        lengths.append(np.sqrt(dx*dx+dy*dy))

        if(what=="BoundingBox"):
            invertmax = 2**(self.invertmaxval(cs))-1
            """
            Calculate lengths and
		    // need to check only is NS^2 == len0^2+len1^2 and EW^2=len2^2+len3^2
		    // 1cm dev must lead to fail. So lengths (12cm) minus 11cm. Scaling is from (sid+10)/sid.
		    // find magnification
		    // we want a cm dev to drop confidence to 0.5, so all lengths (emperically) need to be reduced by magnified 11cm
            """
            for (x0,y0) in roipts:
                if cs.pixeldataIn[int(x0)][int(y0)] == 0 or cs.pixeldataIn[int(x0)][int(y0)] ==invertmax: # on annotation
                    return 0

            confidence = 1.
                # ns via nes east+north
            length1 = lengths[0]*lengths[0]+lengths[1]*lengths[1]
            length2 = lengths[4]*lengths[4]
            if(length2>length1):
                confidence *= length1/length2
            else:
                confidence *=length2/length1

                # ns via west+south
            length1 = lengths[2]*lengths[2]+lengths[3]*lengths[3]
            if(length2>length1):
                confidence *= length1/length2
            else:
                confidence *= length2/length1

            # ew via north+west
            length1 = lengths[1]*lengths[1]+lengths[2]*lengths[2]
            length2 = lengths[5]*lengths[5]
            if(length2>length1):
                confidence *= length1/length2
            else:
                confidence *= length2/length1

            # ew via east+south
            length1 = lengths[3]*lengths[3]+lengths[0]*lengths[0]
            if(length2>length1):
                confidence *= length1/length2
            else:
                confidence *= length2/length1

            # punish for not being equal to input boxsize
            for p in range(0,4):
                confidence *= min(2*redlength,lengths[p])/max(2*redlength,lengths[p])
            confidence = confidence **6.

        if(what == "MTF"):
#            colimit = .1 # just a scaling for confidence
            confidence = 1.
            # order lengths and see if the ratio's are within proper limits
            lengths = sorted(lengths)
            for p in reversed(range(0,4)):
                lengths[p] /= lengths[0]
            #	    grlen = [1.000, 1.2281, 1.4177, 1.4302] # pehamed
            #    grlen = [1.000, 1.2344, 1.4305, 1.4400] # wellhofer
            grlen = [1.000, 1.2312, 1.4241, 1.4351] # avg pehamed and wellhofer
            for p in range(0,4):
                confidence *= (1.0 - np.abs(lengths[p]-grlen[p])/grlen[p] )

        print(what+"Confidence = ", (confidence*100.),"%")
        return confidence

#----------------------------------------------------------------------
    def XRayField(self,cs,roipts=None):
        """
        Use either min across line between box and edge, or use corner values
        """
        error = False
        if roipts is None:
            roipts = cs.po_roi

        cs.xrayNSWEmm = []
        # north- and southside
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'N'))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'S'))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'W'))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'E'))
        if min(cs.xrayNSWEmm)<1.:
            error = True
        else:
            error = False
        print('Edge [N/S/W/E] cm = %.1f %.1f %.1f %.1f' % (cs.xrayNSWEmm[0]/10., cs.xrayNSWEmm[1]/10., cs.xrayNSWEmm[2]/10., cs.xrayNSWEmm[3]/10. ))
        cs.xr_roi = []
        xco,yco = self.phantomposmm2pix(roipts,-cs.xrayNSWEmm[2],cs.xrayNSWEmm[0])
        cs.xr_roi.append([xco,yco])
        xco,yco = self.phantomposmm2pix(roipts,cs.xrayNSWEmm[3],cs.xrayNSWEmm[0])
        cs.xr_roi.append([xco,yco])
        xco,yco = self.phantomposmm2pix(roipts,cs.xrayNSWEmm[3],-cs.xrayNSWEmm[1])
        cs.xr_roi.append([xco,yco])
        xco,yco = self.phantomposmm2pix(roipts,-cs.xrayNSWEmm[2],-cs.xrayNSWEmm[1])
        cs.xr_roi.append([xco,yco])

        return error

    def FindXRayEdge(self,cs,roipts,side):
        invertmax = 2**(self.invertmaxval(cs))-1

        pix2phantommm = self.pixToGridScaleCm(cs)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]

        outvalue = cs.forceRoom.outvalue
        # for DiDi, just take the minimal corner value
        if outvalue<0:
            outvalue = min(cs.pixeldataIn[0][0], cs.pixeldataIn[-1][0],cs.pixeldataIn[0][-1],cs.pixeldataIn[-1][-1])
        if cs.mustbeinverted:
            outvalue = invertmax-outvalue
        """
        0 ll [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
        1 ul [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
        2 ur [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
        3 lr [int(immidx+rad/2+.5),int(immidy-rad/2+.5)] ]
        """
        # north- and southside
        baseidmax = 1
        baseidmin = 0
        horizontal = True
        startmax = False
        if side == 'N':
            baseids = [ [1,0],[2,3] ]
            horizontal = False
            startmax = False
        elif side == 'S':
            baseids = [ [1,0],[2,3] ]
            horizontal = False
            startmax = True
        elif side == 'E':
            baseids = [ [3,0],[2,1] ]
            horizontal = True
            startmax = True
        else: #(side == 'W'):
            baseids = [ [3,0],[2,1] ]
            horizontal = True
            startmax = False

        found = False
        edgemm = []
        for ba in baseids:
            baseidmax = ba[0]
            baseidmin = ba[1]
            if horizontal:
                dypx = 1.*(roipts[baseidmax][1]-roipts[baseidmin][1])/(roipts[baseidmax][0]-roipts[baseidmin][0]) # diff in x if 1 px to y
            else:
                dxpy = 1.*(roipts[baseidmax][0]-roipts[baseidmin][0])/(roipts[baseidmax][1]-roipts[baseidmin][1]) # diff in x if 1 px to y
            posvec = []
            valvec = []
            id = 0
            inrange = True
            while inrange:
                if horizontal:
                    if startmax:
                        xpos = id
                    else:
                        xpos = -id
                    ypos = dypx*xpos
                else:
                    if startmax:
                        ypos = id
                    else:
                        ypos = -id
                    xpos = dxpy*ypos

                pos = np.sqrt(xpos*xpos+ypos*ypos)

                # start from maxpoint, and increase with double linear interpolation
                if startmax:
                    x0 = roipts[baseidmax][0] +int(xpos)
                    y0 = roipts[baseidmax][1] +int(ypos)
                else:
                    x0 = roipts[baseidmin][0] +int(xpos)
                    y0 = roipts[baseidmin][1] +int(ypos)
                if xpos<0:
                    x0 -= 1
                x1 = x0+1
                if x0<0 or x1>(widthpx-1) or y0<0 or y0>(heightpx-1):
                    inrange = False
                    if len(valvec)==0:
                        xa = 0 if x0<0 else widthpx-1 
                        ya = 0 if y0<0 else heightpx-1
                        posvec.append(pos)
                        valvec.append(int(cs.pixeldataIn[xa,ya]))
                    break

                val00 = int(cs.pixeldataIn[int(x0),int(y0)])
                val10 = int(cs.pixeldataIn[int(x1),int(y0)])
                val05 = 1.*val00+(xpos-(int)(xpos))*(val10-val00)
                if cs.mustbeinverted:
                    val05 = invertmax-val05
                posvec.append(pos)
                valvec.append(val05)
                if val00 == outvalue:
                    break
                id += 1

            minval = min(valvec)
            meanval = np.mean(valvec)
            threshLow = (9.*minval+meanval)/10.
            threshHigh = outvalue
            if cs.verbose:
                plt.figure()
                plt.plot(posvec,valvec)
                plt.title(side+" "+str(threshLow)+" "+str(threshHigh))
                cs.hasmadeplots = True

            found = False
            for p,v in zip(posvec,valvec):
                if v<threshLow or (threshHigh>threshLow and v>=threshHigh):
                    found = True
                    edgemm.append( p*pix2phantommm+self.boxradmm )
                    if cs.verbose:
                        plt.plot(p,v,'bo')
                    break
            if not found:
                edgemm.append(0)
        
        return max(edgemm)

#----------------------------------------------------------------------
    def HorizontalUniformity(self,cs,roipts_orig=None):
        error = True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        invertmax = 2**(self.invertmaxval(cs))-1
        # Make box 1 cm from edges of x-ray E and W, from 7cm to 9cm mark on N

        # box runs from (xmm,ymm) = (edgeW+10,70) to (edgeE-10,90)
        xmidpx = 0.
        ymidpx = 0.
        for rp in roipts_orig:
            xmidpx += rp[0]
            ymidpx += rp[1]
        xmidpx /= len(roipts_orig)
        ymidpx /= len(roipts_orig)
        xmidpx = int(xmidpx)
        ymidpx = int(ymidpx)
        yhipx = ymidpx-int(self.phantommm2pix(cs,90.))
        ylopx = ymidpx-int(self.phantommm2pix(cs,70.))
        seppx = int(self.phantommm2pix(cs,10.)+.5)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels

        outvalue = cs.forceRoom.outvalue
    # for DiDi, just take the minimal corner value
        if outvalue<0:
            outvalue = min(cs.pixeldataIn[0][0], cs.pixeldataIn[-1][0],cs.pixeldataIn[0][-1],cs.pixeldataIn[-1][-1])
        if cs.mustbeinverted:
            outvalue = invertmax-outvalue

        # Try lefthand side first
        hiyarr = []
        loyarr = []
        for x in range(0,xmidpx):
            if cs.mustbeinverted:
                hiyarr.append(invertmax-cs.pixeldataIn[x,yhipx])
                loyarr.append(invertmax-cs.pixeldataIn[x,ylopx])
            else:
                hiyarr.append(cs.pixeldataIn[x,yhipx])
                loyarr.append(cs.pixeldataIn[x,ylopx])
        minval = min(hiyarr)
        meanval = np.mean(hiyarr)
        threshLow = (9.*minval+meanval)/10.
        threshHigh = outvalue
        for p,v in zip(reversed(range(0,xmidpx)),reversed(hiyarr)):
            if v<=threshLow or (threshHigh>threshLow and v>=threshHigh):
                x1px = p
                break
        minval = min(loyarr)
        meanval = np.mean(loyarr)
        threshLow = (9.*minval+meanval)/10.
        threshHigh = outvalue
        for p,v in zip(reversed(range(0,xmidpx)),reversed(loyarr)):
            if v<=threshLow or (threshHigh>threshLow and v>=threshHigh):
                x2px = p
                break
        xlopx = max(x1px,x2px)+seppx

        # Now righthand side
        hiyarr = []
        loyarr = []
        for x in range(xmidpx,widthpx):
            if cs.mustbeinverted:
                hiyarr.append(invertmax-cs.pixeldataIn[x,yhipx])
                loyarr.append(invertmax-cs.pixeldataIn[x,ylopx])
            else:
                hiyarr.append(cs.pixeldataIn[x,yhipx])
                loyarr.append(cs.pixeldataIn[x,ylopx])
        minval = min(hiyarr)
        meanval = np.mean(hiyarr)
        threshLow = (9.*minval+meanval)/10.
        threshHigh = outvalue
        for p,v in zip(range(xmidpx,widthpx),hiyarr):
            if v<=threshLow or (threshHigh>threshLow and v>=threshHigh):
                x1px = p
                break
        minval = min(loyarr)
        meanval = np.mean(loyarr)
        threshLow = (9.*minval+meanval)/10.
        threshHigh = outvalue
        for p,v in zip(range(xmidpx,widthpx),loyarr):
            if v<=threshLow or (threshHigh>threshLow and v>=threshHigh):
                x2px = p
                break
        xhipx = min(x1px,x2px)-seppx

        """
	    step 2 run Uniformity check and report only roiuniformity
        """
        cs.unif.roi = [ [xlopx,yhipx],[xhipx,yhipx],[xhipx,ylopx],[xlopx,ylopx] ]
        error = self.Uniformity(cs,cs.unif.roi)

        return error

    def Uniformity(self,cs,roipts=None,what="",bshowplot=False):
        """ Concept:
        1. Cut out rectangular roi
        2. mask out the copper grid
        3. calculate line profile
        4. Check for difference between left and right
        """
        error = True
        if roipts is None:
            roipts = cs.unif.roi

        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        # 1. cut out roi
        xmin=roipts[0][0]
        xmax=roipts[1][0]
        ymin=roipts[1][1]
        ymax=roipts[2][1]
        for rp in roipts:
            xmin = min(xmin,rp[0])
            xmax = max(xmax,rp[0])
            ymin = min(ymin,rp[1])
            ymax = max(ymax,rp[1])

        if cs.mustbeinverted:
            smallimage = invertmax-cs.pixeldataIn[xmin:xmax+1,ymin:ymax+1].astype(float)
        else:
            smallimage = cs.pixeldataIn[xmin:xmax+1,ymin:ymax+1].astype(float)
        wid = smallimage.shape[0]
        hei = smallimage.shape[1]
        # 2. mask out copper grid (afterwards copper = 1)
        mask = wadwrapper_lib.threshold_adaptive(smallimage, int(self.phantommm2pix(cs,50)), method='gaussian', offset=10, mode='reflect', param=None)
        enhance = 3
        if what == "lowcontrast":
            enhance = 1

        for i in range (0,enhance):
            mask = scind.binary_erosion(mask)

        if cs.verbose:
            plt.figure()
            plt.imshow(smallimage)
            plt.title('Uniformity image')
            plt.figure()
            plt.imshow(mask)
            plt.title('Uniformity Mask')

        # 3. calculate line profile
        posval = []
        intens = []
        mincount = min(9,hei/3)

        mask = ~mask # masked = true = do NOT include
        mx = ma.masked_array(smallimage,mask)
        mx2 = ma.masked_array((smallimage.astype(float))**2,mask)
        BKcount = mx.count()
        BKmean  = mx.mean()*BKcount
        BKsumsq = mx2.mean()*BKcount
        counts = mx.count(axis=1)
        sums   = mx.sum(axis=1).data
        for ix in range(0,wid):
            if counts[ix]>mincount:
                value = sums[ix]/counts[ix]
                intens.append(value)
                posval.append(self.pix2phantomm(cs,ix))

        if len(posval)==0:
            print("ERROR: Uniformity: no valid pixels found.")
            return error

        if BKcount>0:
            BKmean /= BKcount
            BKsdev = np.sqrt(BKsumsq/BKcount-BKmean*BKmean)
        if what == "lowcontrast":
            cs.unif.BKmean = BKmean
            cs.unif.BKsdev = BKsdev
            return

        cufraction = 1.*(wid*hei-BKcount)/(wid*hei)
        if cufraction<.1 or cufraction>.9:
            print("ERROR: Uniformity: invalid Cu fraction ",cufraction)
            return error

        print(cufraction)
        if cs.verbose or bshowplot==True:
            plt.figure()
            plt.plot(posval,intens)
            plt.title('Uniformity')
            plt.xlabel('positon on phantom [mm]')
            plt.ylabel('average intensity')
            cs.hasmadeplots = True

        ## line uniformity
        intemin = np.min(intens)
        intemax = np.max(intens)
        inteavg = np.mean(intens)
        overlengthuniformity = np.max([intemax-inteavg,inteavg-intemin])/inteavg

        ## left/right uniformity
        inteavgL = 0.
        count = 0
        cutpos = self.pix2phantomm(cs,hei)
        for pos,val in zip(posval,intens):
            if pos<cutpos:
                inteavgL += val
                count += 1
        if count>0:
            inteavgL /= count

        inteavgR = 0.
        count = 0
        cutpos = self.pix2phantomm(cs,wid-hei)
        for pos,val in zip(posval,intens):
            if(pos>cutpos):
                inteavgR += val
                count += 1
        if count>0:
            inteavgR /= count
        LRuniformity = 1
        if inteavgR+inteavgL>0:
            LRuniformity = np.abs(inteavgR-inteavgL)/(inteavgR+inteavgL)
        ## ROI uniformity
        inteavgROI = 0.
        count = 0
        cutpos1 = self.pix2phantomm(cs,(wid-hei)/2.)
        cutpos2 = self.pix2phantomm(cs,(wid+hei)/2.)
        for pos,val in zip(posval,intens):
            if pos>cutpos1 and pos<cutpos2:
                inteavgROI += val
                count += 1
        inteavgROI /= count
        ROIuniformity = np.max([np.abs(inteavgR-inteavgROI),np.abs(inteavgROI-inteavgL)])/inteavgROI
        print("LineUniformity%=",100.*overlengthuniformity)
        print("LRuniformity%=",100.*LRuniformity)
        print("ROIuniformity%=",100.*ROIuniformity)
        print("AAPMROIlimit%=",10)

        cs.unif.ROIuniformity = ROIuniformity
        cs.unif.LRuniformity = LRuniformity
        cs.unif.BKmean = BKmean
        cs.unif.BKsdev = BKsdev
        cs.unif.posval = copy.deepcopy(posval)
        cs.unif.intens = copy.deepcopy(intens)

        """
		// Report
		/* note AAPM_39 states:
		   For film output (hard-copy), optical densities are measured in the center of each quadrant
		   of the film and in the center position to determine absolute density and spatial uniformity.
		   Central film density is acceptable if within 0.10 OD of the programmed OD value (usually
		   1.20). Spatial uniformity is acceptable when all measured OD values are within ~10% of the
		   average OD. For soft-copy evaluation of the images on a workstation, the average digital value
		   of each ROI should be within 10% of the global average. Standard deviation should also be sim-
		   ilar in each of the five ROIs.
		   SD/Av < 5%
		*/
        """
        error = False
        return error
#----------------------------------------------------------------------
    def CuWedge(self,cs, roipts_orig=None):
        # 1. Make box around wedge
        # 2. Find all steps and do some statistcs
        error = True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)
        bump = 2.
        xmmll = -55.+bump
        ymmur = -55.+bump
        xmmur =  55.-bump
        ymmll = -40.-bump
        if cs.forceRoom.phantom == lit.stWellhofer:
            xmmll = -55.+4.
            ymmur = -60.+3.
            xmmur =  55.-4.
            ymmll = -43.-2.

        xpxll,ypxll = self.phantomposmm2pix(roipts_orig,xmmll,ymmll) # lowerlef
        xpxur,ypxur = self.phantomposmm2pix(roipts_orig,xmmur,ymmur) # upperright
        xpxul,ypxul = self.phantomposmm2pix(roipts_orig,xmmll,ymmur) # upperleft
        xpxlr,ypxlr = self.phantomposmm2pix(roipts_orig,xmmur,ymmll) # lowerright

        # box needs to be horizontal
        xlo = int(.5+ max(xpxll,xpxul))
        xhi = int(.5+ min(xpxlr,xpxur))
        ylo = int(.5+ max(ypxlr,ypxll))
        yhi = int(.5+ min(ypxur,ypxul))
        if ylo>yhi:
            print("[CuWedge]: Error, phantom angle too large, cannot sample wedge")
            return error

        # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)
        cs.cuwedge.roi = [ [xlo,yhi],[xhi,yhi],[xhi,ylo],[xlo,ylo] ]
        error = self.AnalyseWedge(cs,cs.cuwedge.roi)
        return error

    def AnalyseWedge(self,cs,roipts_orig=None):
        """
        Concept:
            1. Cut out rectangular roi
            2. For each cu step, measure mean and sd
            4. Calculate SNR and CNR
            5. Calculate kV

        NB excel data fit shows perfectly linear pehamed Cu intensity with {0.55, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00} (F10 after recalib 20090720)
        to make the 0.00 Cu lie on the line, all thicknesses other than 0.00 should be increased by 0.125mm
        the wellhofer wkz table1 data perfectly agrees with {0.57, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00} including zero point

        Assuming that we have a linear read out (as it should be for Qc) meaning pixel value ~ -mu x
        By fitting a linear line through the mean pixelvalues of each step, we obtain I = a x + b
        If we compare the found coefficient a with the reference value which was obtained for a calibrated kV, we can
        calculate the current kV by linear interpolation from the NIST X-RAY Data for Cu; find the kV for which the ratio a/a_ref is the
        same as mu/mu_ref perfect predict if use 0.6 instead of 0.55

        Contrast from Cu Mean/SD?
        """
        error = True
        if roipts_orig is None:
            roipts_orig = cs.cuwedge.roi

        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        iec_wellhofer = [0.6, 0.5,  0.4,  0.3,  0.2,  0.1,  0.]   # mm Cu; wellhofer
#        din_pehamed   = [2.3, 1.85, 1.40, 1.00, 0.65, 0.30, 0.00] # mm Cu; IBA, Digi13, pehamed?
        #    1. Cut out rectangular roi
        xmin = roipts_orig[0][0]
        xmax = roipts_orig[1][0]
        ymin = roipts_orig[2][1]
        ymax = roipts_orig[1][1]
        if ymin>ymax:
            print("[AnalyseWedge]: Error, phantom angle too large, cannot sample wedge")
            return error

        if cs.mustbeinverted:
            smallimage = invertmax-cs.pixeldataIn[xmin:xmax+1,ymin:ymax+1].astype(float)
        else:
            smallimage = cs.pixeldataIn[xmin:xmax+1,ymin:ymax+1].astype(float)

        wid = smallimage.shape[0]
        hei = smallimage.shape[1]
        if cs.verbose:
            plt.figure()
            plt.imshow(smallimage)
            plt.title('Cu Wedge')

        # 2.1 make a profile
        profile = np.zeros(wid,dtype=float)
        for ix in range(0,wid):
            for iy in range(0,hei):
                profile[ix] += smallimage[ix,iy]
            profile[ix]/=hei

        profile = scind.gaussian_filter1d(profile,sigma=3.,order=1)
        if cs.verbose:
            plt.figure()
            plt.plot(profile)
            plt.title('Cu Wedge derivative')

        # 2.2 Find the edges between the steps
        n_edges = 6
        posedges = []
        flatpix = int( self.phantommm2pix(cs,2.5)+.5)
        for ix in range(0,n_edges):
            profminid = np.unravel_index(profile.argmax(), profile.shape)[0]
            posedges.append(profminid)
            miniy = max(0,posedges[-1]-flatpix)
            maxiy = min(wid-1,posedges[-1]+flatpix)
            flatval = .5*(profile[maxiy]+profile[miniy])
            for iy in range(miniy,maxiy+1):
                profile[iy] = flatval
        posedges = sorted(posedges)

        cs.cuwedge.wedge_confidence = 1.
        cs.cuwedge.wedge_confidence *= (1.*min(n_edges,len(posedges))/max(n_edges,len(posedges)))**3. # must be 6! 
        avg_dist = 1.*(posedges[-1]-posedges[0])/(len(posedges)-1)
        for ix in range(1,len(posedges)):
            dist = 1.*(posedges[ix]-posedges[ix-1])
            cs.cuwedge.wedge_confidence *= min(avg_dist,dist)/max(avg_dist,dist)

        if cs.verbose:
            print("Edge 0 at ",posedges[0])
            for ix in range(1,n_edges):
                print("Edge ",ix," at ",posedges[ix]," sep= ",posedges[ix]-posedges[ix-1])

        # 2.3 Calculate statistics for each step
        cs.cuwedge.roi_mean = []
        cs.cuwedge.roi_sdev = []

        xlo = 0
        ylo = 0   # flatpix
        yhi = hei # - flatpix
        if ylo>yhi:
            print("[AnalyseWedge]: Error, phantom angle too large, cannot sample wedge")
            return error

        cs.cuwedge.step_rois = []
        for p in posedges:
            xlo += flatpix
            xhi = p-flatpix
            value = np.mean(smallimage[xlo:xhi,ylo:yhi])
            cs.cuwedge.roi_mean.append( value )
            cs.cuwedge.roi_sdev.append( np.std(smallimage[xlo:xhi,ylo:yhi]) )
            roipts = [ [xlo+xmin,yhi+ymin],[xhi+xmin,yhi+ymin],[xhi+xmin,ylo+ymin],[xlo+xmin,ylo+ymin] ]
            cs.cuwedge.step_rois.append(roipts)
            xlo = p
        xlo += flatpix
        xhi = wid-flatpix
        value = np.mean(smallimage[xlo:xhi,ylo:yhi])

        cs.cuwedge.roi_mean.append( value )
        cs.cuwedge.roi_sdev.append( np.std(smallimage[xlo:xhi,ylo:yhi]) )
        roipts = [ [xlo+xmin,yhi+ymin],[xhi+xmin,yhi+ymin],[xhi+xmin,ylo+ymin],[xlo+xmin,ylo+ymin] ]
        cs.cuwedge.step_rois.append(roipts)

        cs.cuwedge.roi_mmcu = []
        cs.cuwedge.roi_snr = []
        cs.cuwedge.roi_cnr = []
        for ix in range(0,n_edges+1):
            cs.cuwedge.roi_snr.append( cs.cuwedge.roi_mean[ix]/cs.cuwedge.roi_sdev[ix] )
            #roi_cnr.append( np.sqrt(2.*(roi_mean[ix]-roi_mean[n_edges])**2/(roi_sdev[ix]**2+roi_sdev[n_edges]**2) )
            if ix<n_edges:
                cs.cuwedge.roi_cnr.append( np.abs(cs.cuwedge.roi_mean[ix]-cs.cuwedge.roi_mean[ix+1])/np.sqrt(0.5*(cs.cuwedge.roi_sdev[ix]**2+cs.cuwedge.roi_sdev[ix+1]**2)) )
            else:
                cs.cuwedge.roi_cnr.append(0.)
            cs.cuwedge.roi_mmcu.append( iec_wellhofer[ix] )

        cs.cuwedge.dynamicRange = max(cs.cuwedge.roi_mean[n_edges]/cs.cuwedge.roi_mean[0],cs.cuwedge.roi_mean[0]/cs.cuwedge.roi_mean[n_edges])

        if cs.verbose:
            print("mmCu","SNR","CNR")
            for m,s,c in zip(cs.cuwedge.roi_mmcu,cs.cuwedge.roi_snr,cs.cuwedge.roi_cnr):
                print(m,s,c)

            #cu_cs.guesskVp = guesskVp
        """
	    guesskVp = EstimatekV(roi_mean, roi_mmcu,give_info,cs);

        """
        error = False

        return error
#----------------------------------------------------------------------
    def MTF(self,cs,roipts_orig=None):
        error = True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        x18px,y18px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy18mm[0],cs.forceRoom.xy18mm[1])
        x06px,y06px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy06mm[0],cs.forceRoom.xy06mm[1])
        x46px,y46px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy46mm[0],cs.forceRoom.xy46mm[1])
        x14px,y14px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy14mm[0],cs.forceRoom.xy14mm[1])

        roipts = [ [x18px,y18px],[x06px,y06px],[x14px,y14px],[x46px,y46px] ]
#        print(roipts)
#        error = False
#        confid = 1.
#        cs.lastimage = self.templateMatchDisc(cs)
        error, confid = self.AlignROI(cs,roipts,"MTF")

        cs.mtf.roi = roipts
        cs.mtf.pos_confidence = confid

        if not error:
            error = self.AnalyseMTF(cs,roipts)
        return error

    def AnalyseMTF(self,cs,roipts_orig=None):
        """
        Find rotation angle of linepairs insert. Cut out the back rotated insert.
        For each linepair, do an analysis
        """
        error = True
        invertmax = 2**(self.invertmaxval(cs))-1

        if roipts_orig is None:
            roipts_orig = cs.mtf.roi

        extend18 = self.phantommm2pix(cs,2.8)  # extend bbox beyond dot in '1.8' [mm]
        extend46 = self.phantommm2pix(cs,3.2)  # extend beyond dot in '4.6' [mm]
        print("2.8",extend18)
        print("3.2",extend46)
        # First cut out rotated roi
        id18 = 0
        id06 = 1
        id14 = 2
        id46 = 3
        len1846 = np.sqrt((roipts_orig[id18][0]-roipts_orig[id46][0])**2+(roipts_orig[id18][1]-roipts_orig[id46][1])**2)
        print("1846=",len1846)
        extend18 = 0.0786*len1846
        extend46 = 0.0898*len1846
        copyimage = cs.pixeldataIn.astype(float)
        rotanglerad = 3.*np.pi/2.-.5*(np.arctan2((roipts_orig[id18][1]-roipts_orig[id46][1]),(roipts_orig[id18][0]-roipts_orig[id46][0]))+np.arctan2((roipts_orig[id06][1]-roipts_orig[id14][1]),(roipts_orig[id06][0]-roipts_orig[id14][0])))
        rotanglerad += self.adjustmtfangledeg/180*np.pi
        rotangledeg = (rotanglerad/np.pi*180.)
        print("MTF at",rotangledeg, "degrees")
        rotimage = scind.interpolation.rotate(copyimage, rotangledeg, axes=(1, 0), reshape=False, output=None, order=3, mode='constant', cval=0.0, prefilter=True)

        costerm = np.cos(rotanglerad)
        sinterm = np.sin(rotanglerad)
        xc = cs.pixeldataIn.shape[0]/2.
        yc = cs.pixeldataIn.shape[1]/2.
        roipts = []
        for rp in roipts_orig:
            xp = xc +(rp[0]-xc)*costerm-(rp[1]-yc)*sinterm
            yp = yc +(rp[0]-xc)*sinterm+(rp[1]-yc)*costerm
            roipts.append([xp,yp])
        cs.mtf.dotxys = copy.deepcopy(roipts)

        roipts[id18][1] -= extend18
        roipts[id46][1] += extend46
        minxco = roipts[0][0]
        minyco = roipts[0][1]
        maxxco = roipts[0][0]
        maxyco = roipts[0][1]
        for rp in roipts:
            minxco = min(minxco,rp[0])
            maxxco = max(maxxco,rp[0])
            minyco = min(minyco,rp[1])
            maxyco = max(maxyco,rp[1])

        if cs.mustbeinverted:
            smallimage = invertmax-rotimage[int(minxco):int(maxxco)+1,int(minyco):int(maxyco)+1]
        else:
            smallimage = rotimage[int(minxco):int(maxxco)+1,int(minyco):int(maxyco)+1]

        for rp in roipts:
            rp[0] -= minxco
            rp[1] -= minyco


        if self.bShowMTFDetail or cs.verbose:
            plt.figure()
            plt.imshow(smallimage)
            plt.title("MTF Detail")
            cs.hasmadeplots = True

        # Calculate contrast for each line pair
        contrast_freqs = [0.0,
                          0.6, 0.7, 0.8, 0.9,
                          1.0, 1.2, 1.4, 1.6,
                          1.8, 2.0, 2.2, 2.5,
                          2.8, 3.1, 3.4, 3.7,
                          4.0, 4.3, 4.6, 5.0 ]

        num_freq = len(contrast_freqs)
        contrast_response = [0.]*num_freq
        contrast_high = [0.]*num_freq
        contrast_low  = [0.]*num_freq
        contrast_tops = [0]*num_freq
        contrast_bots = [0]*num_freq
        calc_freq     = [0]*num_freq

        startpos = []
        endpos = []
        startpos,endpos = self.FillMTFBarDetails(cs,smallimage)

        # Find high contrast response of line patterns
        for vpi in range(1,num_freq):
            contrast_response[vpi] = 0.
            contrast_high[vpi]     = 0.
            contrast_low[vpi]      = 0.
            contrast_tops[vpi]     = 0
            contrast_bots[vpi]     = 0
            calc_freq[vpi]         = 0.
            contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = self.AnalyseMTF_Part(cs,smallimage, startpos[vpi-1],endpos[vpi-1], vpi)

        if contrast_tops[0] == 0:
            contrast_response[0] = contrast_response[1]
            contrast_high[0] = contrast_high[1]
            contrast_low[0] = contrast_low[1]

        # Correct for Nyquist
#        fnyq = (0.5/self.dcmInfile.PixelSpacing[0])
        fnyq = (0.5/self.pix2phantomm(cs,1.))
        for vpi in range(0,num_freq):
            if contrast_high[vpi]+contrast_low[vpi]>1e-6:
                contrast_response[vpi] = (contrast_high[vpi]-contrast_low[vpi])/(contrast_high[vpi]+contrast_low[vpi])
            else:
                contrast_response[vpi] = 0.

            #check if not beyond Nyquist
            if contrast_freqs[vpi]>fnyq:
                contrast_response[vpi] = 0

        # Correct for consistency
        crLimitabs = self.crLimit*contrast_response[0]
        for vpi in range(1,num_freq-1):
            if(contrast_response[vpi]<crLimitabs and contrast_response[vpi-1]>crLimitabs and contrast_response[vpi+1]>crLimitabs):
                contrast_response[vpi] = (contrast_response[vpi-1] + contrast_response[vpi+1])/2.

        if contrast_response[0]<1.e-6:
            print("Error in MTF: Rotated image?")
            return error

        ctfmtf = self.CTFtoMTF(cs,contrast_freqs,contrast_response)

        if cs.verbose or self.bShowCTF:
            contrast = []
            for ca in contrast_response:
                contrast.append(ca/contrast_response[0])
            plt.figure()
            plt.plot(contrast_freqs,contrast)
            plt.xlabel("Frequency [lp/mm]")
            plt.ylabel("Relative contrast")
            plt.title("Contrast")
            cs.hasmadeplots = True

        #acos = np.abs(np.cos(rotanglerad))
        #mtf_aapm = acos*(0.5/self.pix2phantomm(1.))*.9 # should ignore acos, because I calculate after rotation, so there
        mtf_aapm = (0.5/self.pix2phantomm(cs,1.))*.9 # should ignore acos, because I calculate after rotation, so there

        # Test how well calculated frequencies match with given frequencies
        maxid = len(calc_freq)
        for id in reversed(range(2,len(calc_freq))):
            #if( (calc_freq[id]<1e-6 or np.abs(calc_freq[id]- calc_freq[id-1])<1.e-6) ):
            if calc_freq[id]<1e-6:
                maxid = id
        print("maxid:",maxid)
        if maxid<5 and not self.bIgnoreMTFError:
            print("Error in MTF: Rotated image?")
            return error
        slope, intercept, r_value, p_value, std_err = stats.linregress(contrast_freqs[0:maxid],calc_freq[0:maxid])
        if r_value**2<0.7:
            print("maxid:",maxid)
            for co,ca in zip(contrast_freqs,calc_freq):
                print(co,ca)
        # To get coefficient of determination (r_squared)
#        print("slope:",slope)
#        print("intercept:",intercept)
#        print("maxid:",maxid)
#        print("r-squared:", r_value**2)
        mtf_freq_confidence = 1.*min(slope,1.)/max(slope,1.)*r_value**2
        # at least 10 freqs must be found
        needfound = 10
        found  = 0
        for t,b in zip (contrast_tops,contrast_bots):
            # at least 10 freqs to be found
            if (t == 3 and b == 2):
                found += 1
        mtf_freq_confidence *= 1.*min(found,needfound)/needfound
        # no holes in freq list before the 10th
        first_error = len(contrast_freqs)-1
        for i in range(1,len(contrast_freqs)-1):
            t0 = contrast_tops[i]
            t1 = contrast_tops[i+1]
            b0 = contrast_bots[i]
            b1 = contrast_bots[i+1]
            if( (t0!=3 or b0!=2) and (t1==3 and b1==2) ):
                first_error = i
                break
        mtf_freq_confidence *= 1.*min(first_error,needfound)/needfound
        print("mtf_freq_confidence:",mtf_freq_confidence)
        if mtf_freq_confidence<.7:
            print("found/first_error/needfound:",found,first_error,needfound)
            print("slope/r2:",slope,r_value**2)
            if cs.verbose:
                plt.figure()
                plt.plot(contrast_freqs[0:maxid],calc_freq[0:maxid],'bo')
                plt.title("found vs given freq")
                cs.hasmadeplots = True

    #       print("confid:",mtf_found_confidence)
        cs.mtf.mtf_aapm = mtf_aapm
        cs.mtf.contrast_freqs    = copy.deepcopy(contrast_freqs)
        cs.mtf.contrast_response = copy.deepcopy(contrast_response)
        cs.mtf.contrast_high     = copy.deepcopy(contrast_high)
        cs.mtf.contrast_low      = copy.deepcopy(contrast_low)
        cs.mtf.contrast_tops     = copy.deepcopy(contrast_tops)
        cs.mtf.contrast_bots     = copy.deepcopy(contrast_bots)
        cs.mtf.ctfmtf            = copy.deepcopy(ctfmtf)
        cs.mtf.calc_freqs        = copy.deepcopy(calc_freq)
        cs.mtf.freq_confidence   = mtf_freq_confidence
        error = False
        return error

    def FillMTFBarDetails(self,cs,smallimage):
        """
        Make lists of all start/end positions of each linepair
        Just scale everything to a manual measurement
        Must already get inverted smallimage
        """
        startpos = []
        endpos = []

        wid = smallimage.shape[0]
        hei = smallimage.shape[1]

        washeightLo = 256.
        waswidthLo  = 438.

        # fill positions with manual values; x is for later
        startpos.append([0,1])   # absolute
        startpos.append([0,36])  # relative
        startpos.append([0,74])  # relative
        startpos.append([0,109]) # relative
        startpos.append([0,148]) # relative
        startpos.append([0,185]) # relative
        startpos.append([0,213]) # relative
        startpos.append([0,239]) # relative
        startpos.append([0,264]) # relative

        startpos.append([0,15])  # relative
        startpos.append([0,41])  # relative
        startpos.append([0,66])  # relative
        startpos.append([0,99])  # relative
        startpos.append([0,131]) # relative
        startpos.append([0,151]) # relative
        startpos.append([0,170]) # relative
        startpos.append([0,197]) # relative
        startpos.append([0,224]) # relative
        startpos.append([0,244]) # relative
        startpos.append([0,264]) # relative
        startpos.append([0,hei-1]) # absolute

        for sp in startpos:
            endpos.append([0,0])

        # start with low frequencies
        ipbar = smallimage[int(0.75*wid):int(0.85*wid),0:hei]
        if cs.verbose:
            plt.figure()
            plt.imshow(ipbar)
            plt.title('low frequencies')

        # project all and find xmin and xmax
        pwid = ipbar.shape[0]
        phei = ipbar.shape[1]
        pattern = np.zeros(phei,dtype=float)
        for y in range(0,phei):
            for x in range(0,pwid):
                pattern[y]+= ipbar[x,y]
            pattern[y] /= pwid

        Imin = min(pattern)
        Imax = max(pattern)

        threshold = .5*(Imin+Imax)
        ymin = 1
        for y in range(2,phei-1):
            if pattern[y] <= threshold and pattern[y+1]>threshold:
                ymin = y
                break
        ymax = phei-1
        for y in reversed(range(2,phei-1)):
            if pattern[y] <= threshold and pattern[y-1]>threshold:
                ymax = y
                break
        ymax = min(ymax,phei-2)

        #print("0000: w/ymin/ymax/yy = ", phei, "/",ymin,"/",ymax,"/",ymax-ymin)

        startpos[0][0] = int(.5+280./waswidthLo*wid)
        endpos[0][0]   = int(.5+(280.+70.)/waswidthLo*wid)
        # for higher: divide by %
        for vpi in range(1,8):
            startpos[vpi][1] = ymin+int(0.5+startpos[vpi][1]/washeightLo*(ymax-ymin))
            endpos[vpi-1][1] = startpos[vpi][1]
            startpos[vpi][0] = startpos[0][0]
            endpos[vpi][0]   = endpos[0][0]
        endpos[7][1] = ymin+int(0.5+startpos[8][1]/washeightLo*(ymax-ymin))

        # for low: divide by %
        startpos[8][0] = int(.5+100./waswidthLo*wid)
        endpos[8][0]   = int(.5+(100.+70.)/waswidthLo*wid)
        startpos[8][1] = startpos[0][1]
        for vpi in range(9,20):
            startpos[vpi][1] = ymin+int(0.5+startpos[vpi][1]/washeightLo*(ymax-ymin))
            endpos[vpi-1][1] = startpos[vpi][1]
            startpos[vpi][0] = startpos[8][0]
            endpos[vpi][0]   = endpos[8][0]
        endpos[19][1] = startpos[20][1]

        return startpos,endpos

    def AnalyseMTF_Part(self,cs,smallimage, startpos,endpos, vpi):
        """
        Determine contrast in one bar pattern element
            contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = self.AnalyseMTF_Part(smallimage, startpos[vpi-1],endpos[vpi-1], vpi)
        Already inverted smallimage as input
        """
        ipbar = smallimage[startpos[0]:endpos[0],startpos[1]:endpos[1]]

        contrast_response = 0.
        contrast_high     = 0.
        contrast_low      = 0.
        contrast_tops     = 0
        contrast_bots     = 0
        calc_freq = 0

        pwid = ipbar.shape[0]
        phei = ipbar.shape[1]

        pattern = np.zeros(phei,dtype=float)
        for y in range(0,phei):
            for x in range(0,pwid):
                pattern[y]+= ipbar[x,y]
            pattern[y] /= pwid

        if pattern.shape[0]<2:
            print("[AnalyseMTF_Part] SKIP: no pattern left")
            return contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq

        # 1. find abs min and max
        # 2. between (max+min)/2 point left and right
        # 3. find all tops and valleys
        mustplot,tops,bots = self.FindExtrema(cs,pattern)

        # determine 100%
        if vpi == 1 and len(tops)>0:
            xmin = endpos[0]
            length = smallimage.shape[0]-xmin
            ytop = tops[0]+startpos[1]
            hsize = max(1,int(.5+self.phantommm2pix(cs,0.75)/2.)  )
            if((ytop+hsize)>(smallimage.shape[1]-1)):
                print("[AnalyseMTF_Part] ERROR: cannot find baseline")
                return contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq

            baseline = []
            for x in range(xmin,xmin+length):
                val =0.
                for y in range(ytop-hsize,ytop+hsize+1):
                    val += smallimage[x,y]
                val /= (2*hsize +1.)
                baseline.append(val)
            basemin = min(baseline)
            basemax = max(baseline)

            contrast_response = basemax-basemin
            contrast_high     = basemax
            contrast_low      = basemin
            contrast_tops     = 1
            contrast_bots     = 1

            if cs.verbose:
                plt.figure()
                plt.plot(baseline)
                plt.title("Baseline")

        if (mustplot == True or (len(tops)!=3 or len(bots)!=2) ) and cs.verbose:
            print("vpi=",vpi," length(pattern)=",len(pattern))
            ybots = []
            for b in bots:
                ybots.append(pattern[b])
            ytops = []
            for t in tops:
                ytops.append(pattern[t])
            plt.figure()
            plt.plot(pattern)
            plt.plot(tops,ytops,'r+')
            plt.plot(bots,ybots,'bo')
            plt.title("Nuthink_"+str(vpi)+"_"+str(len(bots))+"_"+str(len(tops)))
            cs.hasmadeplots = True

        # 4. calculate reponse
        if len(tops)>0 and len(bots)>0:
            hival = 0.
            for y in range(0,min(len(tops),3)):
                hival += pattern[tops[y]]
            hival /= min(len(tops),3)

            loval = 0.
            for y in range(0,min(len(bots),2)):
                loval += pattern[bots[y]]
            loval /= min(len(bots),2)

            contrast_response = hival-loval
            contrast_high = hival
            contrast_low = loval
            contrast_tops = len(tops)
            contrast_bots = len(bots)

        calc_freq = 0.
        if len(tops) == 3 and len(bots)==2: # check location of extrema
            order_fine = True
            if not(tops[0]<tops[1] and tops[1]<tops[2]):
                order_fine = False
            if not(bots[0]<bots[1]):
                order_fine = False
            if not (bots[0]>tops[0] and bots[0]<tops[1]):
                order_fine = False
            if not (bots[1]>tops[1] and bots[1]<tops[2]):
                order_fine = False
            if not (order_fine):
                pass
            else:
                halflambda = 0.25*(tops[2]-tops[0])
                calc_freq = .5/self.pix2phantomm(cs,halflambda)
        if cs.verbose:
            print("Found",len(tops)," tops and",len(bots)," bottoms. Contrast=",contrast_response)
#        print(vpi,contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq)
        return contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq

    def FindExtrema(self,cs,pattern):
        """"
        1. Find abs min and max
        2. Between (max+min)/2 point left and right
        3. find all tops and valleys
        """
        mustplot = False
        if(len(pattern)<2):
            print("SKIP: no pattern left")
            tops = []
            bots = []
            return mustplot,tops,bots

        phei = len(pattern)

        sigma = self.sigma_ext
        sigma_max = (len(pattern)-1)/6.
        sigma = min(sigma,sigma_max)
        xderiv1 = np.array([])
        if(sigma<0.45):
            xderiv1 = mymath.FiniteDifference1D(pattern,order=1)
        else:
            xderiv1 = scind.filters.gaussian_filter1d(pattern, sigma, order=1)

        Imin = min(pattern)
        Imax = max(pattern)

        threshold = .5*(Imin+Imax)
        ymin = 1
        for y in range(2,phei-1):
            if(pattern[y] <= threshold and pattern[y+1]>threshold):
                ymin = y
                break
        ymax = phei-1
        for y in reversed(range(2,phei-1)):
            if(pattern[y] <= threshold and pattern[y-1]>threshold):
                ymax = y
                break
        ymax = min(ymax,phei-2)
        yminmax = [ymin,ymax]
        mustplot,tops,bots = self.FindAllExtrema(cs,pattern,xderiv1,yminmax)

        # check if too many, then try with larger gaussian
        goon = False
        while( ( len(tops)>3 or len(bots)>2 ) and goon == False):
            if cs.verbose:
                print("Too many extrema; rerun with larger sigma")
            tops_bk = copy.deepcopy(tops)
            bots_bk = copy.deepcopy(bots)
            mustplot_bk = mustplot
            sigma *=2
            if(sigma > sigma_max):
                sigma = sigma_max
                goon = True  # break from loop

            xderiv1 = scind.filters.gaussian_filter1d(pattern, sigma, order=1)
            mustplot,tops,bots = self.FindAllExtrema(cs,pattern,xderiv1,yminmax)

            if(len(tops)<3 or len(bots)<2): # restore previous
                tops = copy.deepcopy(tops_bk)
                bots = copy.deepcopy(bots_bk)
                mustplot = mustplot_bk
                goon = True  # break from loop

        # check if too few, then smaller gaussian
        goon = False
        while( ( len(tops)<3 or len(bots)<2 ) and goon == False):
            if cs.verbose:
                print("Too few extrema; rerun with smaller sigma")
            tops_bk = copy.deepcopy(tops)
            bots_bk = copy.deepcopy(bots)
            mustplot_bk = mustplot

            sigma /=2.
            xderiv1 = np.array([])
            if(sigma<0.45):
                xderiv1 = mymath.FiniteDifference1D(pattern,order=1)
                goon = True
            else:
                xderiv1 = scind.filters.gaussian_filter1d(pattern, sigma, order=1)
            mustplot,tops,bots = self.FindAllExtrema(cs,pattern,xderiv1,yminmax)

            if(len(tops)>3 or len(bots)>2): # restore previous
                tops = copy.deepcopy(tops_bk)
                bots = copy.deepcopy(bots_bk)
                mustplot = mustplot_bk
                goon = True  # break from loop

        if cs.verbose:
            print("[FindExtrema]C ntops/nbots = ",len(tops),"/",len(bots))
        return mustplot,tops,bots

    def FindAllExtrema(self,cs,pattern,xderiv1,yminmax):
        mustplot = False

        ymin = yminmax[0]
        ymax = yminmax[1]

        # look for sign changes of deriv1
        tops = []
        bots = []
        for y in range(ymin,ymax+1):
            if(xderiv1[y]>=0. and xderiv1[y+1]<0.):
                if(pattern[y+1]>pattern[y]):
                    tops.append(y+1)
                else:
                    tops.append(y)
                if(len(tops)>3):
                    print("[FindAllExtrema] ntops>3! Using only first 3.")
                    mustplot = True
            if(xderiv1[y]<=0. and xderiv1[y+1]>0.):
                if(pattern[y+1]<pattern[y]):
                    if(len(tops)>0):
                        bots.append(y+1)
                else:
                    if(len(tops)>0):
                        bots.append(y)
                if(len(bots)>2):
                    print("[FindAllExtrema] nbots>2! Using only first 2.")
                    mustplot = True

        if mustplot: # SOMETHING WRONG, INGORE XDERIV AND JUST LOOK AT MIN/MAX
            tops = []
            bots = []
            mustplot = False

            for y in range(ymin+1,ymax):
                if(pattern[y]>pattern[y+1] and pattern[y]>pattern[y-1]):
                    tops.append(y)
                    if(len(tops)>3):
                        print("[FindAllExtrema] ntops2>3! Using only first 3.")
                        mustplot = True
                if(pattern[y]<pattern[y+1] and pattern[y]<pattern[y-1]):
                    if(len(tops)>0):
                        bots.append(y)
                    if(len(bots)>2):
                        print("[FindAllExtrema] nbots2>2! Using only first 2.")
                        mustplot = True

        return mustplot,tops,bots

    def CTFtoMTF(self,cs,freq_in, ctf):
        # Med Phys. 2008 Jan;35(1):270-9.
        # First remove double entries
        freq = []
        for f in range(0,len(freq_in)):
            freq.append(freq_in[f])
            if (f>0):
                if(freq[f] == freq[f-1]):
                    freq[f] += .01

        # here ctf is (max-min)/(max+min); need ctf/ctf(0)
        # with the latter ctf, mtf is given as mtf(u) = pi/4 [ctf(u)+1/3ctf(3u)-1/5ctf(5u)+1/7ctf(7u)-1/9ctf(9u)....]

        # fit a 3rd order polynomial to CTF:
        coeffs = np.polyfit(freq, ctf, deg=3)
        poly = np.poly1d(coeffs)
        #        print(poly)
        mtf = copy.deepcopy(ctf)
        zerocor = mtf[1]
        fnyq = (0.5/self.pix2phantomm(cs,1.)) # cut-off frequencies > Nyquist
        N = len(ctf)
        for i in range(1,N):
            # find out how many harmonics are available
            j = 1
            factor = 2*j+1
            prev_harm =mtf[0]
            while True:
                if cs.verbose:
                    print(freq[i],i,j,factor*freq[i],fnyq)
                harmonic = np.polyval(poly, factor*freq[i])
                if(harmonic <=0. or harmonic>prev_harm or factor*freq[i]>fnyq ):
                    break
                if(int(j/2)==int((j+1)/2)):
                    factor = -factor
                mtf[i] += harmonic/factor
                j += 1
                factor = 2*j+1
                prev_harm = harmonic

            if(ctf[0]>1e-6):
                mtf[i] *= np.pi/4./ctf[0]

        # mtf[0] will not get any harmonic corrections, this is a nice estimate
        zerocor = mtf[1]-zerocor
        mtf[0]+= zerocor
        return mtf

    def CTFtoMTFNoFit(self,freq_in, ctf):
        # First remove double entries
        freq = []
        for f in range(0,len(freq_in)):
            freq.append(freq_in[f])
            if (f>0):
                if(freq[f] == freq[f-1]):
                    freq[f] += .01

        # here ctf is (max-min)/(max+min); need ctf/ctf(0)
        # with the latter ctf, mtf is given as mtf(u) = pi/4 [ctf(u)+1/3ctf(3u)-1/5ctf(5u)+1/7ctf(7u)-1/9ctf(9u)....]
        mtf = copy.deepcopy(ctf)
        N = len(ctf)
        for i in range(0,N):
            # find out how many harmonics are available
            j = 1
            factor = 2*j+1
            while(factor*freq[i]<freq[N-1] and j<5):
                harmonic = mymath.linearInterExtrapolate(freq,ctf,factor*freq[i])
                if(int(j/2)==int((j+1)/2)):
                    factor = -factor
                mtf[i] += harmonic/factor
                j += 1
                factor = 2*j+1
            if ctf[0]>1e-6:
                mtf[i] *= np.pi/4./ctf[0]

        return mtf

#----------------------------------------------------------------------
    def LowContrast(self,cs,roipts_orig=None):
        """
        for each disc, calc av in box in center and in box 1 cm above it (without Cu grid)
        Calculate CNR: (av(Cu)-av(bk))/sqrt((sd(Cu)^2+sd(bk)^2)/2)
        """
        error =True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        invertmax = 2**(self.invertmaxval(cs))-1
        low_cnr = []
        mean_sg = []
        mean_bk = []
        sdev_sg = []
        sdev_bk = []

        # pehamed
        mmCu = [ [-25.,-35.], [-15.,-35.], [-05.,-35.], [05.,-35.] ]
        CuDimmm = 04.
        yrefmm = 10. # offset to background

        if cs.forceRoom.phantom == lit.stWellhofer:
            mmCu = [ [36.4,22.8], [31.,17.5], [25.,11.9], [19.3,06.] ]
            CuDimmm = 02.5
            yrefmm = 5.5
        else:
            peha = XRayStruct(cs.dcmInfile,cs.pixeldataIn,cs.forceRoom)

        cs.loco.lo_rois = []
        for cu in mmCu:
            # first calculate statistics for low contrast object
            llx,lly = self.phantomposmm2pix(roipts_orig,cu[0]-CuDimmm/2.,cu[1]-CuDimmm/2.)
            urx,ury = self.phantomposmm2pix(roipts_orig,cu[0]+CuDimmm/2.,cu[1]+CuDimmm/2.)
            lrx,lry = self.phantomposmm2pix(roipts_orig,cu[0]+CuDimmm/2.,cu[1]-CuDimmm/2.)
            ulx,uly = self.phantomposmm2pix(roipts_orig,cu[0]-CuDimmm/2.,cu[1]+CuDimmm/2.)
            xlo = int(max(llx,ulx)+.5)
            yhi = int(max(lly,lry)+.5)
            xhi = int(min(urx,lrx)+.5)
            ylo = int(max(uly,ury)+.5)
            roipts = [ [xlo,ylo],[xlo,yhi],[xhi,yhi],[xhi,ylo] ]
            cs.loco.lo_rois.append(roipts)
            if cs.mustbeinverted:
                smallimage = invertmax - cs.pixeldataIn[xlo:xhi+1,ylo:yhi+1].astype(float)
            else:
                smallimage = cs.pixeldataIn[xlo:xhi+1,ylo:yhi+1].astype(float)
            roi_mean_s = np.mean(smallimage)
            roi_sdev_s = np.std(smallimage)
            mean_sg.append(roi_mean_s)
            sdev_sg.append(roi_sdev_s)

            # next calculate statistics for background
            ylo -= int(.5+self.phantommm2pix(cs,yrefmm))
            yhi -= int(.5+self.phantommm2pix(cs,yrefmm))
            roipts = [ [xlo,ylo],[xlo,yhi],[xhi,yhi],[xhi,ylo] ]
            cs.loco.lo_rois.append(roipts)
            if cs.forceRoom.phantom == lit.stWellhofer:
                if cs.mustbeinverted:
                    smallimage = invertmax - cs.pixeldataIn[xlo:xhi+1,ylo:yhi+1].astype(float)
                else:
                    smallimage = cs.pixeldataIn[xlo:xhi+1,ylo:yhi+1].astype(float)
                roi_mean_bk = np.mean(smallimage)
                roi_sdev_bk = np.std(smallimage)
            else:
                self.Uniformity(peha,roipts,"lowcontrast")
                roi_mean_bk = peha.unif.BKmean
                roi_sdev_bk = peha.unif.BKsdev

            mean_bk.append(roi_mean_bk)
            sdev_bk.append(roi_sdev_bk)
            low_cnr.append((roi_mean_s-roi_mean_bk)/np.sqrt(0.5*(roi_sdev_s**2+roi_sdev_bk**2)))
            if cs.verbose:
                print("mean fg/bk=",roi_mean_s,"/",roi_mean_bk)
                print("sdev fg/bk=",roi_sdev_s,"/",roi_sdev_bk)

        cs.loco.low_cnr = copy.deepcopy(low_cnr)
        cs.loco.mean_sg = copy.deepcopy(mean_sg)
        cs.loco.mean_bk = copy.deepcopy(mean_bk)
        cs.loco.sdev_sg = copy.deepcopy(sdev_sg)
        cs.loco.sdev_bk = copy.deepcopy(sdev_bk)

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
            dicomfields = [
                #{'key':"0008,0021",'name','value':0, 'quantity','level':,'rank':},
                {'key':"0008,0021",  'name':"SeriesDate"},
                {'key':"0008,0031",  'name':"SeriesTime"},
                {'key':"0008,0070",  'name':"Manufacturer"},
                {'key':"0008,0080",  'name':"InstitutionName"},
                {'key':"0008,1010",  'name':"StationName"},
                {'key':"0008,1030",  'name':"StudyDescription"},
                {'key':"0008,103E",  'name':"SeriesDescription"},
                {'key':"0008,1070",  'name':"Operator's Name"},
                {'key':"0010,0020",  'name':"PatientID"},
                {'key':"0018,0015",  'name':"BodyPartExamined"},
                {'key':"0018,0060",  'name':"kVp"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,1004",  'name':"PlateID"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)"},
                {'key':"0018,1150",  'name':"ExposureTime (ms)"},
                {'key':"0018,1152",  'name':"Exposure (mAs)"},
                {'key':"0018,115E",  'name':"ImageAreaDoseProduct"},
                {'key':"0018,1160",  'name':"FilterType"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"0018,1166",  'name':"Grid"},
                {'key':"0018,1190",  'name':"FocalSpot(s)"},
                {'key':"0018,1200",  'name':"Date of Last Calibration"},
                {'key':"0018,1260",  'name':"PlateType"},
                {'key':"0018,1400",  'name':"AcquisitionDeviceProcessingDescription"},
                {'key':"0018,1401",  'name':"AcquisitionDeviceProcessingCode"},
                {'key':"0018,1403",  'name':"CassetteSize"},
                {'key':"0018,1404",  'name':"ExposuresOnPlate"},
                {'key':"0018,1508",  'name':"PositionerType"},
                {'key':"0018,6000",  'name':"Sensitivity"},
                {'key':"0020,4000",  'name':"ImageComments"},
                {'key':"0028,0006",  'name':"PlanarConfiguration"},
                {'key':"0028,0101",  'name':"BitsStored"}
            ]
            if 'RelativeXRayExposure' in cs.dcmInfile: #DX
                dicomfields.append({'key':"0018,1405",'name':"Relative Exposure"})

        elif(info == "qclight"):
            dicomfields = [
                {'key':"0008,0021",  'name':"SeriesDate"},
                {'key':"0008,0031",  'name':"SeriesTime"},
                {'key':"0008,1070",  'name':"Operator's Name"},
                {'key':"0018,0060",  'name':"kVp"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,1004",  'name':"Plate ID"},
                {'key':"0018,1401",  'name':"Acquisition Device Processing Code"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)"},
                {'key':"0018,1150",  'name':"ExposureTime (ms)"},
                {'key':"0018,1152",  'name':"Exposure (mAs)"},
                {'key':"0018,115E",  'name':"ImageAreaDoseProduct"},
                {'key':"0018,1160",  'name':"FilterType"},
                {'key':"0018,1190",  'name':"FocalSpot(s)"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"0018,1166",  'name':"Grid"},
                {'key':"0018,1200",  'name':"Date of Last Calibration"},
                {'key':"0018,5021",  'name':"Postprocessing"},
                {'key':"0018,6000",  'name':"Sensitivity"}
            ]
            if 'RelativeXRayExposure' in cs.dcmInfile: #DX
                dicomfields.append({'key':"0018,1405",'name':"Relative Exposure"})

        elif(info == "qcwad"):
            offset = -25 # rank must be negative, so recalc as offset+real position
            if not 'DistanceSourceToDetector' in cs.dcmInfile: # WKZ-like fcr
                dicomfields = [
                    {'key':"0008,0021",  'name':"SeriesDate"},
                    {'key':"0008,0031",  'name':"SeriesTime", 'quantity':'time', 'level':1, 'rank':offset+2},
                    {'key':"0008,1070",  'name':"Operator's Name"},
                    {'key':"0018,1004",  'name':"Plate ID"},
                    {'key':"0018,1401",  'name':"Acquisition Device Processing Code", 'quantity':'processing', 'level':1, 'rank':offset+9},
                    {'key':"0018,1020",  'name':"SoftwareVersions"},
                    {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                    {'key':"0018,6000",  'name':"Sensitivity",'quantity':'S','level':1,'rank':offset+12}
                ]
            else:
                dicomfields = [
                    {'key':"0008,0021",  'name':"SeriesDate"},
                    {'key':"0008,0031",  'name':"SeriesTime", 'quantity':'time', 'level':1, 'rank':offset+2}, # spot 1 reserved for stand
                    {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)", 'quantity':'distance', 'level':1, 'rank':offset+3},
                    {'key':"0018,0060",  'name':"kVp", 'level':1, 'rank':offset+4},
                    {'key':"0018,1160",  'name':"FilterType", 'quantity':'filter', 'level':1, 'rank':offset+5},
                    {'key':"0018,1190",  'name':"FocalSpot(s)", 'quantity':'focalspot', 'level':1, 'rank':offset+6},
                    {'key':"0018,1166",  'name':"Grid", 'quantity':'grid', 'level':1, 'rank':offset+7},
                    {'key':"0018,5021",  'name':"Postprocessing", 'quantity':'processing', 'level':1, 'rank':offset+9}, # spot 8 reserved for rotation

                    {'key':"0018,1150",  'name':"ExposureTime (ms)", 'quantity':'ms','level':1,'rank':offset+10},
                    {'key':"0018,1152",  'name':"Exposure (mAs)", 'quantity':'mAs','level':1,'rank':offset+11},
                    {'key':"0018,115E",  'name':"ImageAreaDoseProduct", 'quantity':'DAP','level':1,'rank':offset+12},
                
                    {'key':"0008,1070",  'name':"Operator's Name"},
                    {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                    {'key':"0018,1020",  'name':"SoftwareVersions"},
                    {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                    {'key':"0018,1200",  'name':"Date of Last Calibration"},
                    {'key':"0018,6000",  'name':"Sensitivity",'quantity':'S'}
                ]

        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','rank':missing or a number} )
        results = []
        for df in dicomfields:
            key = df['key']
            value = ""
            replaced = False
            if(key == "0018,1152"):
                key1 = "0018,1153"
                try:
                    value =  cs.dcmInfile[dicom.tag.Tag(key1.split(',')[0],key1.split(',')[1])].value / 1000.
                    replaced = True
                except:
                    value = ""
                    replaced = False
            if(replaced == False):
                try:
                    value =  cs.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
                except:
                    value = ""

            df['value'] = value
            results.append( df )

        return results
#----------------------------------------------------------------------
    def QC(self,cs):
        """
        all in one package!
        """
        error = True
        msg = ""
        """
        Outline:
        1. Determine grid scaling, rotation, location.

        step1: reorder such that UL,UR,LR,LL // was N, E, S W
        step2: travel straight along NS and find edge of x-ray; similar for EW
        step3: find horizontal uniformity
        step4: find Cu wedge stuff
        step5: find resolution stuff
        step6: low contrast stuff
        step7: put it in a report
        """

        error,msg = self.checkPhantomRotation(cs)
        if error:
            msg += "AlignROI(BoundingBox) "

        # 2: travel straight along NS and find edge of x-ray; similar for EW
        if not error:
            error = self.XRayField(cs)
            if error:
                msg += "XRayField "

        # 3: find horizontal uniformity
        if not error:
            cs.unif.verticalROI = False
            error = self.HorizontalUniformity(cs,cs.po_roi)
            if error:
                msg += "HorizontalUniformity "

        # 4: find Cu wedge stuff
        if not error:
            error = self.CuWedge(cs,cs.po_roi)
#            error = self.CuWedgeFlat(cs,cs.po_roi)
            if error:
                msg += "CuWedge "

        # 5: find resolution stuff
        if not error:
            error = self.MTF(cs,cs.po_roi)
            if error:
                msg += "MTF "

        # 6: low contrast stuff
        if not error:
            error = self.LowContrast(cs,cs.po_roi)
            if error:
                msg += "LowContrast "

        return error,msg

    def mAsCalc(self,cs):
        """
        Convenience function to calculate mAs from DOP if not DR
        """
        if not 'ImageAndFluoroscopyAreaDoseProduct' in cs.dcmInfile: # WKZ-like fcr 
            return 0. # no reported DOP

        stand = self.TableOrWall(cs)
        if not 'DistanceSourceToDetector' in cs.dcmInfile:
            if stand == lit.stTable:
                sidmm = cs.forceRoom.sidtablemm
            else:
                sidmm = cs.forceRoom.sidwallmm
        else:
            sidmm = cs.dcmInfile.DistanceSourceToDetector
            if sidmm<1.e-6:
                sidmm = cs.forceRoom.sidtablemm
                if stand == lit.stWall:
                    sidmm = cs.forceRoom.sidwallmm

        sid_m = sidmm/1000.
        kvp = cs.dcmInfile.KVP
        dop = cs.dcmInfile.ImageAndFluoroscopyAreaDoseProduct
        if sid_m<1.e-6 or kvp < 1.e-6 or dop <1.e-6:
            return 0.

        factor =  8900.
        if stand == lit.stWall: # should also adjust area 27^2 for table, but hey.
            factor = 16790. # depends on filter used
        xrayarea = max(1e-6,(cs.xrayNSWEmm[0]+cs.xrayNSWEmm[1])*(cs.xrayNSWEmm[2]+cs.xrayNSWEmm[3])/(1000.*1000.))

        mAs_calc = (dop/np.power(kvp,2.4))*sid_m**2/xrayarea*factor
        return mAs_calc

    def XRayDev(self,cs):
        """
        Convenience function to calculate max deviation of L/R in %
        """
        stand = self.TableOrWall(cs)
        minedge = min(cs.xrayNSWEmm)
        maxedge = max(cs.xrayNSWEmm)
        meanedge = np.mean(cs.xrayNSWEmm)
        if not 'DistanceSourceToDetector' in cs.dcmInfile:
            if stand == lit.stTable:
                sidmm = cs.forceRoom.sidtablemm
            else:
                sidmm = cs.forceRoom.sidwallmm
        else:
            sidmm = cs.dcmInfile.DistanceSourceToDetector
            if sidmm<1.e-6:
                sidmm = cs.forceRoom.sidtablemm
                if stand == lit.stWall:
                    sidmm = cs.forceRoom.sidwallmm
        devedge = 100.*max(np.abs(minedge-meanedge),np.abs(maxedge-meanedge))/sidmm
        if maxedge-meanedge < meanedge-minedge:
            devedge *= -1
        return devedge

    def ReportEntries(self,cs):
        """
        Convenience function to list all calculated items with names, only if cs completely filled
        """
        labvals = []
        stand = self.TableOrWall(cs)

        ## Phantom orientation; level 1 = show by default; level 2 = show in details
        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','rank':missing or a negative number} )
        # if no rank given, the order of addition will be used
        # if no quantity given, 'name' will be used
        # if no level given, the default will be used
        offset = -25 # rank must be negative, so recalc as offset+real position
        labvals.append( {'name':'PhantomOrientation','value':cs.po_rot, 'quantity':'rotate','level':1,'rank':offset+7} )
        labvals.append( {'name':'AlignConfidence','value':100.*cs.bbox_confidence, 'quantity':'aligned','level':2} )
        labvals.append( {'name':'xray[N]cm','value':cs.xrayNSWEmm[0]/10., 'quantity':'xrN','level':2} )
        labvals.append( {'name':'xray[E]cm','value':cs.xrayNSWEmm[3]/10., 'quantity':'xrE','level':2} )
        labvals.append( {'name':'xray[S]cm','value':cs.xrayNSWEmm[1]/10., 'quantity':'xrS','level':2} )
        labvals.append( {'name':'xray[W]cm','value':cs.xrayNSWEmm[2]/10., 'quantity':'xrW','level':2} )
        labvals.append( {'name':'xrayDev','value':self.XRayDev(cs), 'quantity':'xrayDev','level':1, 'rank':offset+12} )

        ## uniformity
        if stand == lit.stWall:
            label = 'LRUniformity'
            value = 100.*cs.unif.LRuniformity
        else:
            label = 'ROIUniformity'
            value = 100.*cs.unif.ROIuniformity
        labvals.append( {'name':label,'value':value, 'quantity':'Uniformity','level':1,'rank':offset+13} )

        ## cuwedge
        labvals.append( {'name':'CuConfidence','value':cs.cuwedge.wedge_confidence*100,'level':2} ) # Confidence in wedge finding
        # SNR max
        labvals.append( {'name':'CuSNR_'+str(cs.cuwedge.roi_mmcu[-1]),'value':cs.cuwedge.roi_snr[-1], 'quantity':'SNR','level':1,'rank':offset+15} )
        # CNR between steps all > 1
        minCNR = cs.cuwedge.roi_cnr[0]
        for i in range(1,len(cs.cuwedge.roi_cnr)-1):
            minCNR = min (minCNR,cs.cuwedge.roi_cnr[i])
        labvals.append( {'name':'CuCNRmin','value':minCNR, 'quantity':'CNRmin','level':2} )
        # Dynamic Range
        labvals.append( {'name':'CuDR'+str(cs.cuwedge.roi_mmcu[0])+'_'+str(cs.cuwedge.roi_mmcu[-1]),'value':cs.cuwedge.dynamicRange, 'quantity':'DynRange','level':1,'rank':offset+14} )
        # guesskVp
        labvals.append( {'name':'kVp_est','value':cs.cuwedge.guesskVp, 'quantity':'kVp_est','level':2} )
        # approximate mAs from DAP until DX info, better approx mAs
        labvals.append( {'name':'mAs_est','value':self.mAsCalc(cs), 'quantity':'mAs_est','level':2} )

        ## low contrast
        for i,cnr in enumerate(cs.loco.low_cnr):
            labvals.append( {'name':'lowCNR_'+str(i),'value':cnr,'level':2} )

        ## mtf
        labvals.append( {'name':'MTFPosConfidence','value':cs.mtf.pos_confidence*100.,'level':2} )
        labvals.append( {'name':'MTFFreqConfidence','value':cs.mtf.freq_confidence*100.,'level':2} )
        labvals.append( {'name':'AreaContrast5','value':mymath.AreaUnderCurve(cs.mtf.contrast_freqs,cs.mtf.contrast_response),'level':2} )
        labvals.append( {'name':'AreaMTF5','value':mymath.AreaUnderCurve(cs.mtf.contrast_freqs,cs.mtf.ctfmtf),'level':2} )
        labvals.append( {'name':'MTF10','value':mymath.MTF10pct(cs.mtf.contrast_freqs,cs.mtf.ctfmtf),'level':2} )

        return labvals

    def saveAnnotatedImage(self,cs,fname):
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        im = scipy.misc.toimage(cs.pixeldataIn.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

        # now draw all rois in reserved color
        rois = []
        rois.append(cs.po_roi) # phantom orientation box
        rois.append(cs.xr_roi) # xray edges
        rois.append(cs.unif.roi) # horizontal uniformity box
        rois.append(cs.cuwedge.roi) # Cu wedge box
        for r in cs.cuwedge.step_rois:
            rois.append(r) # Cu wedge element
        rois.append(cs.mtf.roi) # MTF box
        for r in cs.loco.lo_rois:
            rois.append(r) # low contrast element

        draw = ImageDraw.Draw(im)
        for r in rois:
            roi =[]
            for x,y in r:
                roi.append( (int(x+.5),int(y+.5)))
            draw.polygon(roi,outline=0)
        del draw

        # convert to RGB for JPG, cause JPG doesn't do PALETTE and PNG is much larger
        im = im.convert("RGB")

        imsi = im.size
        if max(imsi)>2048:
            ratio = 2048./max(imsi)
            im = im.resize( (int(imsi[0]*ratio+.5), int(imsi[1]*ratio+.5)),Image.ANTIALIAS)
        im.save(fname)

    def checkPhantomRotation(self,cs,ebbox=None):
        """
        Concept: Try to find out if phantom is rotated over 90, 180 or 270 degrees (WKZ mostly)
        Workflow:
        1. Find center and orientation
        2. Find most likely location of CuWedge and LinePairs
        """
        error = True
        msg = ""
        #ang = 1
        #cs.pixeldataIn = np.rot90(cs.pixeldataIn,ang)

        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        # 1. Find center and orientation
        ang0 = 0
        error = self.findPhantomOrientation(cs,ebbox)
        if error:
            ang0 = 180
            cs.pixeldataIn = np.rot90(cs.pixeldataIn,2) # try at 180 deg rotation
            error = self.findPhantomOrientation(cs,ebbox)

            if error:
                msg += "NoOrientation "
                return error,msg

        # 2. Make box at x-5cm and calculate avg
        roipts_orig = cs.po_roi
        cs.test_rois = []
        avgs = []
        stds = []
        # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)

        if cs.forceRoom.phantom == lit.stWellhofer: # WKZ
            x0 = -15.+4
            y0 = -43.-2.
            x1 =  15.-4
            y1 = -60.+3.
        else:
            bump = 2.
            x0 = -15. +bump
            y0 = -40. -bump
            x1 =  15. +bump
            y1 = -50. -bump

        lims = [ [[ x0, y0],[ x1, y1]], # S
                 [[ y0, x0],[ y1, x1]], #
                 [[ x0,-y1],[ x1,-y0]],
                 [[-y0, x0],[-y1, x1]]
                 ]
        for ix,(mmll,mmur) in enumerate(lims):
            xpxll,ypxll = self.phantomposmm2pix(roipts_orig,mmll[0],mmll[1]) # lowerlef
            xpxur,ypxur = self.phantomposmm2pix(roipts_orig,mmur[0],mmur[1]) # upperright
            xpxul,ypxul = self.phantomposmm2pix(roipts_orig,mmll[0],mmur[1]) # upperleft
            xpxlr,ypxlr = self.phantomposmm2pix(roipts_orig,mmur[0],mmll[1]) # lowerright

            # box needs to be horizontal
            xlo = int(.5+ max(xpxll,xpxul))
            xhi = int(.5+ min(xpxlr,xpxur))
            ylo = int(.5+ max(ypxlr,ypxll))
            yhi = int(.5+ min(ypxur,ypxul))
            if xhi < xlo:
                swap = xlo
                xlo = xhi
                xhi = swap
            if yhi < ylo:
                swap = ylo
                ylo = yhi
                yhi = swap

            # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)
            cs.test_rois.append([ [xlo,yhi],[xhi,yhi],[xhi,ylo],[xlo,ylo] ])
            if cs.mustbeinverted:
                smallimage = invertmax-cs.pixeldataIn[xlo:xhi,ylo:yhi]
            else:
                smallimage = cs.pixeldataIn[xlo:xhi,ylo:yhi]

            avg = np.mean(smallimage)
            std = np.std(smallimage)
            print(ix,avg,std,avg/std)
            avgs.append(avg)
            stds.append(std)

        ang = 0
        cs.po_rot = ang0
        msg = ""
        if min(stds[1],stds[2]) > max(stds[0],stds[3]): # pos was goed
            pass
        elif min(stds[2],stds[3]) > max(stds[1],stds[0]): # pos +90rot
            ang = 1
        elif min(stds[3],stds[0]) > max(stds[2],stds[1]): # pos +180rot
            ang = 2
        elif min(stds[0],stds[1]) > max(stds[3],stds[2]): # pos +270rot
            ang = 3
        else:
            # if heuristics did not work, make it more simple: cu-wedge should have smallest sd
            min_stds = min(stds)
            min_index = stds.index(min_stds)            
            if min_index == 0: # south, that's fine
                ang = 0
            elif min_index == 1: # west
                ang = 1
            elif min_index == 2: # north
                ang = 2
            else:
                ang = 3
#            print("[checkPhantomRotation] ERROR! Cannot find orientation",avgs)
#            return True,"NoOrientation "

        if ang>0:
            cs.pixeldataIn = np.rot90(cs.pixeldataIn,-ang)
            cs.po_rot = 90*ang+ang0
        if ang>0 or ang0 >0:
            msg = 'Rot'+str(cs.po_rot)+' '
        error = self.findPhantomOrientation(cs)
        return error,msg

