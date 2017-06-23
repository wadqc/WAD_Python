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
    20170621: increased uniformity box to deal with no pix found.
    20161220: removed testing stuff; removed class variabled
    20160902: sync with wad2.0; Unified pywad1.0 and wad2.0
    20150616: better orientation module
    20150609: Initial from QCXRay_Lib
"""
__version__ = '20170621'
__author__ = 'aschilham'

import numpy as np
import scipy.ndimage as scind
import operator

try:
    # wad2.0 runs each module stand alone
    import QCDDL_constants as lit
except ImportError:
    from . import QCDDL_constants as lit

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
if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
    raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")

class DDLStruct:
    class Room :
        def __init__ (self,_name, outvalue=-1, 
                      sHmm=-1, sVmm=-1, 
                      pH2mm=-1, pV2mm=-1, 
                      protocolH = '', protocolV = '',
                      phantom=lit.stPehamed,sdthresh=-1):
            self.name = _name        # identifier of room
            self.outvalue = outvalue # value of pixels outside x-ray field
            self.sidHorizontal_mm  = sHmm
            self.sidVertical_mm    = sVmm
            self.p2mm512Horizontal = pH2mm # factor to convert pix in 512 to mm on phantom; from manual fit
            self.p2mm512Vertical   = pV2mm
            self.protocolHorizontal = protocolH
            self.protocolVertical = protocolV
            self.phantom = phantom
            self.sdthresh = sdthresh # threshold on stdev for determin table or wall, now only for CALHOS, maybe WKZ?
            self.skipFFT = False
            if(phantom == lit.stWellhofer):
                self.skipFFT = True # only for wellhofer

    class UnifStruct :
        def __init__ (self):
            self.ROIuniformity = self.LRuniformity = LineUniformity = 0. # fraction
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

    #####################
    roomUnknown = Room(lit.stUnknown)

    def readDICOMtag(self,key): # slice=2 is image 3
        value = wadwrapper_lib.readDICOMtag(key,self.dcmInfileRaw,self.imslice)
        return value

    def maybeInvert(self):
        if self.dcmInfile is None:
            return
        self.mustbeinverted = False

        if self.dcmInfile.PhotometricInterpretation == "MONOCHROME2":
            self.mustbeinverted = True

        # somehow does not give the correct behaviour for DDL
        # hence a more pragmatic approach: 
        # if the mode of the hostogram of pixelvalues of a large part of the image > middle bin, then do invert
        widthpx = np.shape(self.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(self.pixeldataIn)[1]
        sqsize = min(widthpx,heightpx)
        midx = int(.5*(widthpx-1)+.5)
        midy = int(.5*(heightpx-1)+.5)
        sqpart = 5 # was 3 for CR/DR
        smallimage = wadwrapper_lib.extract(self.pixeldataIn, [int(sqsize/sqpart), int(sqsize/sqpart)],[midx,midy])
        hist = np.histogram(smallimage,bins=256)
        peakpos = hist[1][np.argmax(hist[0])]
        if hist[1][-1]-peakpos<peakpos-hist[1][0]:
            self.mustbeinverted = False
        else:
            self.mustbeinverted = True

        print("Must be Inverted", self.mustbeinverted)

    def determineDeviceID(self):
        if self.dcmInfile is None:
            return

        self.guessroom = self.roomUnknown

        stationnames = [
            ["AZUDDL",self.roomF7],
        ]
        dicomfields = [ ["0008,1010",  "Station Name"],["0018,1000","Device Serial Number"]]
        key = dicomfields[0][0]
        dicvalue = self.readDICOMtag(key)
        if dicvalue != "":
            dicvalue = dicvalue.upper()
            for sn in stationnames:
                if dicvalue.find(sn[0]) != -1:
                    self.guessroom = sn[1]
                    print("DetermineDeviceID:", self.guessroom.name, '(', dicvalue, ')')
                    return
        # try one of the DDLrooms;
        didiserialnumber = [
            ["24211830",self.roomF7],
        ]
        key = dicomfields[1][0]
        dicvalue = self.readDICOMtag(key)
        if dicvalue != "":
            for sn in didiserialnumber:
                if  dicvalue.find(sn[0]) != -1:
                    self.guessroom = sn[1]
                    print("DetermineDeviceID:", self.guessroom.name, '(', dicvalue, ')')
                    return

    def __init__ (self,dcmInfile,pixeldataIn,dicomMode):
        self.knownHorizontalOrVertical = None
        self.verbose = False

        # input image
        self.dcmInfileRaw = dcmInfile
        self.pixeldataInRaw = pixeldataIn
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.dicomMode = dicomMode
        self.imslice = 0
        # one name for 2D 3D
        if self.dicomMode != wadwrapper_lib.stMode2D:
            self.imslice = 0 if pixeldataIn is None else int(len(self.pixeldataInRaw)/2)
            self.dcmInfile = self.dcmInfileRaw.info
            self.pixeldataIn = None if pixeldataIn is None else self.pixeldataInRaw[self.imslice]

        self.expertOverridepixToGridScaleCm = None
        self.mustbeinverted = False # if pixval(Cu) = high and pixval(air) = low, then mustbeinverted
        self.lastimage = None # GUI feedback

        # phantom orientation
        self.test_rois = []
        self.po_center_roi = [] # xmid,ymid,rad
        self.po_roi = [] # list [ [x,y] ]
        self.po_rot = 0  # only is phantom is rotate multiple of 90
        self.po_angledeg = 0. # exact angle of rotation after multiple of 90
        self.po_fftroi = []

        # horz. uniformity
        self.unif = self.UnifStruct()

        # Cu Wedge
        self.cuwedge = self.CuStruct()

        # Low Contrast
        self.loco = self.LoCoStruct()
        self.bbox_confidence = 0.

        # for matlib plotting
        self.hasmadeplots = False

        # room names
        self.roomF7 = self.Room("AZUDDL", outvalue=0,sHmm=1500,sVmm=1500,pH2mm=0.6384,pV2mm=0.6384,protocolH='HSG',protocolV='Defaec')
        #roomF7.skipFFT = True

        self.guessroom = self.roomUnknown
        self.determineDeviceID()

        if not pixeldataIn is None:
            self.maybeInvert()

class DDLQC:
    def __init__(self):
        self.qcversion = __version__
        self.boxradmm = 110  # choose 11 cm or 8 cm for clean surroundings

    def readDICOMtag(self,cs,key): # slice=2 is image 3
        value = wadwrapper_lib.readDICOMtag(key,cs.dcmInfileRaw,cs.imslice)
        return value

    def pixToGridScaleCm(self,cs):
        if not cs.expertOverridepixToGridScaleCm is None:
            return cs.expertOverridepixToGridScaleCm

        pixel_spacing_x = 512./cs.dcmInfile.Rows
        stand = self.HorizontalOrVertical(cs)
        # calculate distance = L0-Thick
        if stand == lit.stHorizontal:
            return pixel_spacing_x*cs.guessroom.p2mm512Horizontal

        return pixel_spacing_x*cs.guessroom.p2mm512Vertical

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

        return int(xpx), int(ypx)

#----------------------------------------------------------------------
    def invertmaxval(self,cs):
        # Not stored properly in self.dcmfileIn.BitsStored
        dicomfields = [ ["0028,0101",  "Bits Stored"]]
        key = dicomfields[0][0]
        dicvalue = self.readDICOMtag(cs, key)
        if dicvalue != "":
            return dicvalue
        else:
            return 0

    def HorizontalOrVertical(self,cs):
        if cs.knownHorizontalOrVertical is not None:
            return cs.knownHorizontalOrVertical

        result = lit.stUnknown
        if cs.dcmInfile.ProtocolName == cs.guessroom.protocolHorizontal:
            result =  lit.stHorizontal
        elif cs.dcmInfile.ProtocolName == cs.guessroom.protocolVertical:
            result =  lit.stVertical

        if result != lit.stUnknown:
            cs.knownHorizontalOrVertical = result
        else:
            print('|%s|%s|%s'%(cs.guessroom.protocolHorizontal, cs.guessroom.protocolVertical, cs.dcmInfile.ProtocolName))
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
        sqpart = 5 # was 3 for CR/DR
        if cs.mustbeinverted:
            smallimage = invertmax - wadwrapper_lib.extract(cs.pixeldataIn, [int(sqsize/sqpart),int(sqsize/sqpart)],[midx,midy])
        else:
            smallimage = wadwrapper_lib.extract(cs.pixeldataIn, [int(sqsize/sqpart),int(sqsize/sqpart)],[midx,midy])
        smallimage = scind.gaussian_filter(smallimage, sigma=5)
        cs.lastimage = smallimage

        x0,y0 = np.unravel_index(smallimage.argmin(), smallimage.shape)
        immidx = int(midx-smallimage.shape[0]/2+x0)
        immidy = int(midy-smallimage.shape[1]/2+y0)
        rad=20
        cs.po_center_roi = [immidx,immidy,rad] # feedback for GUI

        # 1.2 position 10x10phantomcm box and manipulate corner positions
        pix2phantommm = self.pixToGridScaleCm(cs)
        # choose 11cm each direction or 8 cm
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
            for i in [80,90,70,60,110]:
                radtest = 2*i/pix2phantommm
                dist = np.abs(radtest-radpts)
                if dist < mindist:
                    mindist = dist
                    self.boxradmm = i
            print("found mindist=", mindist, "rad=", self.boxradmm)
        # try to find rotation of phantom
        if not cs.guessroom.skipFFT:
            error,roipts,rotangledeg = self.FieldRotationFFT(cs,roipts)
            if error: # probably too small angle, so try again
                error,roipts,rotangledeg = self.FieldRotationFFT(cs,roipts,30.)

        else:
            rotangledeg = 0.
            error = False
        cs.po_angledeg = rotangledeg

        cs.bbox_confidence = 0.
        if not error:
            last_attempt= False
            radmmtests = [80,70,90,110,60]#[80,110,90,70,60]
            for radmm in radmmtests:
                print('radmm=%d, '%radmm, end="")
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

        print("using rad,conf:", self.boxradmm, cs.bbox_confidence)
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
            roipts.append([int(xp), int(yp)])
        return roipts

    def _fieldRotationFFT(self,cs,smallimage,initangle=None):
        runmode = 2

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

            kwart = np.zeros((int(fwid/2), int(fhei/2)),dtype=float)
            for x in range(0,int(fwid/2)):
                for y in range (0,fhei/2):
                    kwart[x,y] += psd2D[int(fwid/2)+x,int(fhei/2)+y]
                    kwart[x,y] += psd2D[int(fhei/2)-y,int(fwid/2)+x]

            kwartmax = np.max(kwart)
            # Find local maxima
            posxmax = []
            posymax = []
            kwart2 = np.zeros((int(fwid/2),int(fhei/2)),dtype=float)
            while len(posxmax)<20 and np.max(kwart)>0.05:
                xx,yy = np.unravel_index(kwart.argmax(), kwart.shape)
                if( not(xx<3 and yy<3) ):
                    posxmax.append(xx)
                    posymax.append(yy)
                    kwart2[xx,yy] = kwart[xx,yy]
                kwart[xx,yy] = 0

            # a bit of sorting:
            index = range(len(posxmax))
            index.sort(key = posxmax.__getitem__)
            posxmax[:] = [posxmax[i] for i in index]
            posymax[:] = [posymax[i] for i in index]

            for i,x in enumerate(posxmax):
                if x>posxmax[-1]/10.: # skip low 10 pct
                    posxmax = posxmax[i:]
                    posymax = posymax[i:]
                    break
            maxoff = len(posxmax)-5#
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
                print('orientation fit: %f + %f*x; offset=%d'%(intercept, slope, dummy2))
                print('angledeg=%f,initangledeg=%f;'%(anglerad/np.pi*180., 0 if initangle is None else initangle))
                cs.hasmadeplots = True
        else: # 'new'
            # must yield anglerad and r2
            #Let's log-scale the image so we're dealing with something in uint8 range.
            # Calculate a 2D power spectrum
            abs_data = np.abs( F2 )#**2.
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
            mask[0:int(wid/2),0:int(hei/2)] = True
            mask[int(wid/2):wid,int(hei/2):hei] = True
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
                print('run %d angle: %f (%f), r2=%f'%(i, anglerad, anglerad/np.pi*180., r1_value**2.))
                r2anglenum.append( (r1_value**2.,anglerad,len(rows)))
                mask = ~mask
                if cs.verbose:
                    plt.figure()
                    plt.plot(rows,cols,'bo')
                    dafit = [intercept+slope*x for x in rows]
                    plt.plot(rows,dafit,'b-')
                    plt.title("Orientation")
                    print('orientation fit: %f + %f*x'%(intercept, slope))
                    print('angledeg=%f,initangledeg=%f;'%(anglerad/np.pi*180., 0 if initangle is None else initangle))
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
        if 1<0:
            import os.path
            dumppath = os.path.join(os.path.expanduser('~'),'temp')
            np.save(os.path.join(dumppath,'smallimage.npy'), smallimage)

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
            print("FieldRotationFFT:",label,"confidence too low:", confidence, off)
            #print(offanglerad)
            return error,roipts_orig,rotangledeg
        error = False

        print("rotangledegFFT:", rotangledeg, confidence, off)

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
        searchrad = max(1,searchrad)
        print("%s searchrad="%what, searchrad)
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
                minx = max(0,x0-searchrad)
                maxx = min(widthpx-2,x0+searchrad)
                y0 = rp[1]
                miny = max(0,y0-searchrad)
                maxy = min(heightpx-2,y0+searchrad)
                if cs.mustbeinverted:
                    cropped = workimage[minx:maxx+1,miny:maxy+1]
                else:
                    cropped = invertmax - workimage[minx:maxx+1,miny:maxy+1]
                #plt.figure()
                #plt.imshow(cropped)
                #plt.show()
                # smoothing to get rid of noise
                sigma = max(1.5,cropped.shape[0]/8.)
                cropped = scind.gaussian_filter(cropped, sigma)
                # find maximum location in projections
                xflat = np.mean(cropped,axis=0)
                yflat = np.mean(cropped,axis=1)
                x1 = minx+np.unravel_index(yflat.argmax(), yflat.shape)[0]
                y1 = miny+np.unravel_index(xflat.argmax(), xflat.shape)[0]

                rp[0] = x1
                rp[1] = y1
            roipts = self.ConsistencyAlign(cs,roipts0,roipts,what)
            conf_pts.append((self.ROIConfidence(cs,roipts,what),copy.deepcopy(roipts)))

        # Just take the last one; should be best!
#        conf_pts = sorted(conf_pts) #sorting orientation/distance box does worsen results

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
            print("AlignRoi (",what,"):", label, ", confidence too low:", confidence)

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
		     need to check only is NS^2 == len0^2+len1^2 and EW^2=len2^2+len3^2
		     1cm dev must lead to fail. So lengths (12cm) minus 11cm. Scaling is from (sid+10)/sid.
		     find magnification
		     we want a cm dev to drop confidence to 0.5, so all lengths (emperically) need to be reduced by magnified 11cm
            """
            for (x0,y0) in roipts:
                if cs.pixeldataIn[x0][y0] == 0 or cs.pixeldataIn[x0][y0] ==invertmax: # on annotation
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
            confidence = confidence **2. # 6 is too strict for deformations in DDL; power 3 is max!

        print(what+"Confidence = ", (confidence*100.), "%")
        return confidence

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
        yhipx = ymidpx-int(self.phantommm2pix(cs,95.))
        ylopx = ymidpx-int(self.phantommm2pix(cs,65.))
        seppx = int(self.phantommm2pix(cs,10.)+.5)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels

        outvalue = cs.guessroom.outvalue
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

        meanval = np.mean(hiyarr)
        for p,v in zip(range(0,xmidpx),hiyarr):
            if (outvalue>meanval and v<outvalue) or (outvalue<meanval and v>outvalue):
                x1px = p
                break

        meanval = np.mean(loyarr)
        for p,v in zip(range(0,xmidpx),loyarr):
            if (outvalue>meanval and v<outvalue) or (outvalue<meanval and v>outvalue):
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

        meanval = np.mean(hiyarr)
        for p,v in zip(reversed(range(0,xmidpx)),reversed(hiyarr)):
            if (outvalue>meanval and v<outvalue) or (outvalue<meanval and v>outvalue):
                x1px = p
                break

        meanval = np.mean(loyarr)
        for p,v in zip(reversed(range(0,xmidpx)),reversed(loyarr)):
            if (outvalue>meanval and v<outvalue) or (outvalue<meanval and v>outvalue):
                x2px = p
                break
        xhipx = xmidpx+min(x1px,x2px)-seppx

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
        enhance = 3
        mask_offset = 10
        nbsize = int(self.phantommm2pix(cs,50))
        if what == "lowcontrast":
            option = 1
            if option == 0: # option0: use all, do not look for grid
                mask = ~ma.make_mask_none(np.shape(smallimage))
            elif option == 1:   # option1: ImageJ autothreshold
                histo = np.histogram(smallimage, bins=min(256,int(np.max(smallimage)-np.min(smallimage))))
                level = histo[1][0]+wadwrapper_lib.threshold_isodata2(histo[0])
                mask = smallimage>level
            else: # option2: what I used in QCXRay, but no erosion
                mask = wadwrapper_lib.threshold_adaptive(smallimage, nbsize, method='gaussian', offset=mask_offset, mode='reflect', param=None)
        else:
            mask = wadwrapper_lib.threshold_adaptive(smallimage, nbsize, method='gaussian', offset=mask_offset, mode='reflect', param=None)
            for i in range (0,enhance):
                mask = scind.binary_erosion(mask)

        if cs.verbose:
            plt.figure()
            plt.imshow(smallimage)
            plt.title('Uniformity image'+what)
            plt.figure()
            plt.imshow(mask)
            plt.title('Uniformity Mask'+what)
            cs.hasemadeplots = True

        # 3. calculate line profile
        posval = []
        intens = []
        mincount = min(9, int(hei/3))

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
            print(counts)
            print("ERROR: Uniformity %s: no valid pixels found."%what)
            plt.figure()
            plt.imshow(smallimage)
            plt.title('Uniformity image')
            plt.figure()
            plt.imshow(mask)
            plt.title('Uniformity Mask')
            cs.hasmadeplots = True

        if BKcount>0:
            BKmean /= BKcount
            BKsdev = np.sqrt(BKsumsq/BKcount-BKmean*BKmean)
        if what == "lowcontrast":
            cs.unif.BKmean = BKmean
            cs.unif.BKsdev = BKsdev
            return

        cufraction = 1.*(wid*hei-BKcount)/(wid*hei)
        if cufraction<.1 or cufraction>.9:
            print("ERROR: Uniformity: invalid Cu fraction ", cufraction)
            return error

        print(cufraction)
        if  cs.verbose or bshowplot==True:
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
        print("LineUniformity%=", 100.*overlengthuniformity)
        print("LRuniformity%=", 100.*LRuniformity)
        print("ROIuniformity%=", 100.*ROIuniformity)
        print("AAPMROIlimit%=", 10)

        cs.unif.LineUniformity = overlengthuniformity
        cs.unif.ROIuniformity = ROIuniformity
        cs.unif.LRuniformity = LRuniformity
        cs.unif.BKmean = BKmean
        cs.unif.BKsdev = BKsdev
        cs.unif.posval = copy.deepcopy(posval)
        cs.unif.intens = copy.deepcopy(intens)

        """
		 Report
		 note AAPM_39 states:
		   For film output (hard-copy), optical densities are measured in the center of each quadrant
		   of the film and in the center position to determine absolute density and spatial uniformity.
		   Central film density is acceptable if within 0.10 OD of the programmed OD value (usually
		   1.20). Spatial uniformity is acceptable when all measured OD values are within ~10% of the
		   average OD. For soft-copy evaluation of the images on a workstation, the average digital value
		   of each ROI should be within 10% of the global average. Standard deviation should also be sim-
		   ilar in each of the five ROIs.
		   SD/Av < 5%

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
        if cs.guessroom.phantom == lit.stWellhofer:
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
        # find peaks
        pfactor = 10
        xy_max = []
        while len(xy_max)<n_edges and np.max(profile)/pfactor>.1:
            pfactor += 1
            xy_max,xy_min = wadwrapper_lib.peakdet(profile, delta=np.max(profile)/pfactor) # returns indices for freqs
        while len(xy_max)>n_edges and pfactor>1:
            pfactor -= 1
            xy_max,xy_min = wadwrapper_lib.peakdet(profile, delta=np.max(profile)/pfactor) # returns indices for freqs

        posedges = sorted([x for x,y in xy_max])
        if len(posedges) < n_edges:
            print("[AnalyseWedge]: Error, cannot find 6 edges: %d"%len(posedges))
            return error

        cs.cuwedge.wedge_confidence = 1.
        cs.cuwedge.wedge_confidence *= (1.*min(n_edges,len(posedges))/max(n_edges,len(posedges)))**3. # must be 6! 
        avg_dist = 1.*(posedges[-1]-posedges[0])/(len(posedges)-1)
        for ix in range(1,len(posedges)):
            dist = 1.*(posedges[ix]-posedges[ix-1])
            cs.cuwedge.wedge_confidence *= min(avg_dist,dist)/max(avg_dist,dist)

        if cs.verbose:
            print("Edge 0 at ", posedges[0])
            for ix in range(1,len(posedges)):
                print("Edge ", ix, " at ", posedges[ix], " sep= ", posedges[ix]-posedges[ix-1])
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
        flatpix = int( self.phantommm2pix(cs,2.5)+.5)
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
            print("mmCu", "SNR", "CNR")
            for m,s,c in zip(cs.cuwedge.roi_mmcu,cs.cuwedge.roi_snr,cs.cuwedge.roi_cnr):
                print(m, s, c)

            #cu_cs.guesskVp = guesskVp
        """
	    guesskVp = EstimatekV(roi_mean, roi_mmcu,give_info,cs);

        """
        error = False

        return error
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

        if cs.guessroom.phantom == lit.stWellhofer:
            mmCu = [ [36.4,22.8], [31.,17.5], [25.,11.9], [19.3,06.] ]
            CuDimmm = 02.5
            yrefmm = 5.5
        else:
            peha = DDLStruct(cs.dcmInfileRaw,cs.pixeldataInRaw,cs.dicomMode)

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
            if cs.guessroom.phantom == lit.stWellhofer:
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
                print("mean fg/bk=", roi_mean_s, "/", roi_mean_bk)
                print("sdev fg/bk=", roi_sdev_s, "/", roi_sdev_bk)

        cs.loco.low_cnr = copy.deepcopy(low_cnr)
        cs.loco.mean_sg = copy.deepcopy(mean_sg)
        cs.loco.mean_bk = copy.deepcopy(mean_bk)
        cs.loco.sdev_sg = copy.deepcopy(sdev_sg)
        cs.loco.sdev_bk = copy.deepcopy(low_cnr)

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
                ["0008,0021",  "SeriesDate"],
                ["0008,0031",  "SeriesTime"],
                ["0008,0070",  "Manufacturer"],
                ["0008,0080",  "InstitutionName"],
                ["0008,1010",  "StationName"],
                ["0008,1030",  "StudyDescription"],
                ["0010,0020",  "PatientID"],
                ["0018,0015",  "BodyPartExamined"],
                ["0018,1000",  "DeviceSerialNumber"],
                ["0018,1020",  "SoftwareVersions"],
                ["0018,1030",  "ProtocolName"],
                ["0028,0006",  "PlanarConfiguration"],
                ["0028,0101",  "BitsStored"],
                ["0028,0010",  "Rows"],
                ["0028,0011",  "Columns"]
            ]
        elif(info == "qc"):
            dicomfields = [
                ["0008,0021",  "SeriesDate"],
                ["0008,0031",  "SeriesTime"],
                ["0018,1000",  "DeviceSerialNumber"],
                ["0018,1020",  "SoftwareVersions"],
                ["0018,1030",  "ProtocolName"],
                ["0028,0010",  "Rows"],
                ["0028,0011",  "Columns"]
            ]

        results = []
        for df in dicomfields:
            key = df[0]
            value = ""
            try:
                value = self.readDICOMtag(cs,key)
            except:
                value = ""

            results.append( (df[1],value) )

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

        # 6: low contrast stuff
        if not error:
            error = self.LowContrast(cs,cs.po_roi)
            if error:
                msg += "LowContrast "

        return error,msg

    def ReportEntries(self,cs):
        """
        Convenience function to list all calculated items with names, only if cs completely filled
        """
        labvals = []
        stand = self.HorizontalOrVertical(cs)

        ## Phantom orientation
        labvals.append( ('PhantomOrientation',cs.po_rot) )
        labvals.append( ('PhantomAngle',cs.po_angledeg) )
        labvals.append( ('AlignConfidence',100.*cs.bbox_confidence) )

        ## uniformity
        labvals.append( ('LineUniformity',100.*cs.unif.LineUniformity) )
        labvals.append( ('LRUniformity',100.*cs.unif.LRuniformity) )
        labvals.append( ('ROIUniformity',100.*cs.unif.ROIuniformity) )

        ## cuwedge
        labvals.append( ('CuConfidence',cs.cuwedge.wedge_confidence*100) ) # Confidence in wedge finding
        labvals.append( ('Cu0mean',cs.cuwedge.roi_mean[-1]) ) # mean
        labvals.append( ('Cu0sd',cs.cuwedge.roi_sdev[-1]) ) # sdev
        # SNR max
        labvals.append( ('CuSNR_'+str(cs.cuwedge.roi_mmcu[-1]),cs.cuwedge.roi_snr[-1]) )
        # CNR between steps all > 1
        minCNR = cs.cuwedge.roi_cnr[0]
        for i in range(1,len(cs.cuwedge.roi_cnr)-1):
            minCNR = min (minCNR,cs.cuwedge.roi_cnr[i])
        labvals.append( ('CuCNRmin',minCNR) )
        # Dynamic Range
        labvals.append( ('CuDR'+str(cs.cuwedge.roi_mmcu[0])+'_'+str(cs.cuwedge.roi_mmcu[-1]),cs.cuwedge.dynamicRange) )
        # guesskVp
        labvals.append( ('guesskVp',cs.cuwedge.guesskVp) )

        ## low contrast
        for i,cnr in enumerate(cs.loco.low_cnr):
            labvals.append( ('lowCNR_'+str(i),cnr) )

        return labvals

    def saveAnnotatedImage(self,cs,fname):
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        im = scipy.misc.toimage(cs.pixeldataIn.transpose()+1,low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

        # now draw all rois in reserved color
        rois = []
        rois.append(cs.po_roi) # phantom orientation box
        rois.append(cs.unif.roi) # horizontal uniformity box
        rois.append(cs.cuwedge.roi) # Cu wedge box
        for r in cs.cuwedge.step_rois:
            rois.append(r) # Cu wedge element
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
                #cs.pixeldataIn = np.rot90(cs.pixeldataIn,2) # reset
                msg += "NoOrientation "
                return error,msg

        # 2. Make box at x-5cm and calculate avg
        roipts_orig = cs.po_roi
        cs.test_rois = []
        avgs = []
        stds = []
        # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)

        if cs.guessroom.phantom == lit.stWellhofer: # WKZ
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
            print(ix, avg, std, avg/std)
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

