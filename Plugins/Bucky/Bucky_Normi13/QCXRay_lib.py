__author__ = 'aschilha'
"""
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED! And make sure rescaling is corrected!

TODO:
Changelog:
    20160205: Distinguish between linepairs insert typ38 and RXT02
    20160202: added uniformity
    20151109: start of new module, based on QCXRay_lib of Bucky_PEHAMED_Wellhofer of 20151029
"""
import dicom
import numpy as np
import scipy.ndimage as scind
import QCXRay_constants as lit
import QCXRay_math as mymath
from pyWADLib import wadwrapper_lib
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
if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
    raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")
import QCUniformity_lib

class Room :
    name = ""     # identifier of room
    outvalue = -1 # value of pixels outside x-ray field
    pidtablemm = -1    # distance between mid phantom and detector in mm
    pidwallmm = -1
    sidtablemm = -1
    sidwallmm = -1
    skipFFT = True # for now 

    xy06mm = [] # x,y position in mm of decimal dot in 0.6 lp/mm 
    xy10mm = [] # x,y position in mm of decimal dot in 1.0 lp/mm 
    
    def __init__ (self,_name, outvalue=-1, tablesid=-1, wallsid=-1, tablepid=-1, wallpid=-1,
                  linepairmarkers = {},artefactborderpx=[0,0,0,0],detectorname={}):
        self.name = _name
        self.outvalue = outvalue
        self.pidtablemm = tablepid
        self.pidwallmm  = wallpid
        self.sidtablemm = tablesid
        self.sidwallmm = wallsid
        self.artefactborderpx = artefactborderpx
        self.detector_name = detectorname # a dict of [detectorid] = name like 
        if len(linepairmarkers)>0:
            self.linepairmodel = linepairmarkers['type']
            if self.linepairmodel == 'RXT02':
                self.xy06mm = linepairmarkers['mm0.6']
                self.xy10mm = linepairmarkers['mm1.0']
            elif self.linepairmodel == '38':
                self.xy06mm = linepairmarkers['mm0.6']
                self.xy14mm = linepairmarkers['mm1.4']
                self.xy18mm = linepairmarkers['mm1.8']
                self.xy46mm = linepairmarkers['mm4.6']
            elif self.linepairmodel == 'None':
                pass
            else:
                raise ValueError('[Room] Unknown linepairmodel')
        
    def setPIDs(self,_pidtable, _pidwall):
        self.pidtablemm = _pidtable
        self.pidwallmm = _pidwall
    def setSIDS(self,_sidtable, _sidwall):
        self.sidtablemm = _sidtable
        self.sidwallmm = _sidwall

class XRayStruct:
    verbose = False
    knownTableOrWall = None

    roomUnknown = Room(lit.stUnknown)
    forceRoom = roomUnknown

    class CuStruct :
        roi_snr = []
        roi_cnr = []
        roi_mmcu = []
        roi_mean = []
        roi_sdev = []
        dynamicRange = -1
        guesskVp     = -1
        slope        = 0.
        roi = []
        step_rois = []
        wedge_confidence = -1.
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
        low_cnr = []
        mean_sg = []
        mean_bk = []
        sdev_sg = []
        sdev_bk = []
        lo_rois = []
        def __init__ (self):
            self.low_cnr = []
            self.mean_sg = []
            self.mean_bk = []
            self.sdev_sg = []
            self.sdev_bk = []
            self.lo_rois = []
            self.lo_rois_bku = []
            self.lo_rois_bkd = []

    class MTFStruct:
        dotxys = []
        contrast_response = []
        contrast_high = []
        contrast_low = []
        ctfmtf = []
        contrast_freqs = []
        contrast_tops = []
        contrast_bots = []
        calc_freqs = []  # frequencies in phantom units as measured from extremes
        mtf_aapm = -1.
        freq_confidence = -1.
        pos_confidence  = -1.
        roi = []
        def __init__ (self):
            self.dotxys = []
            self.contrast_response = []
            self.contrast_high = []
            self.contrast_low = []
            self.ctfmtf = []
            self.contrast_freqs = []
            self.contrast_tops = []
            self.contrast_bots = []
            self.mtf_aapm = -1.
            self.freq_confidence = -1.
            self.pos_confidence = -1.
            self.roi = []

    #####################
    # input image
    dcmInfile = None
    pixeldataIn = None
    mustbeinverted = False # if pixval(Cu) = high and pixval(air) = low, then mustbeinverted
    bbox_confidence = 0.

    # phantom orientation
    po_center_roi = [] # xmid,ymid,rad
    po_roi = [] # list [ [x,y] ]
    po_fftroi = []
    po_rot = 0
    test_rois = []

    # xray field
    xrayNSWEmm = []
    xr_roi = []

    # Cu Wedge
    cuwedge = None

    # MTF
    mtf = None

    # Low Contrast
    loco = None

    # Uniformity
    unif = None
    
    lastimage = None # GUI feedback

    # for matlib plotting
    hasmadeplots = False

    def maybeInvert(self):
        if self.dcmInfile is None:
            return
        self.mustbeinverted = False
        if self.dcmInfile.PhotometricInterpretation == "MONOCHROME2":
            self.mustbeinverted = True
        print "Must be Inverted",self.mustbeinverted


    def __init__ (self,dcmInfile,pixeldataIn,room):
        self.verbose = False
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.hasmadeplots = False
        self.mustbeinverted = False
        self.forceRoom = room
        self.lastimage = None
        self.po_center_roi = []
        self.po_roi = []
        self.po_rot = 0
        self.po_fftroi = []
        self.test_rois = []

        self.xrayNSWEmm = []
        self.xr_roi = []

        self.cuwedge = self.CuStruct()
        self.mtf = self.MTFStruct()
        self.loco = self.LoCoStruct()
        self.bbox_confidence = 0.

        self.maybeInvert()

class XRayQC:
    qcversion = 20160205

    boxradmm   = 90  #
    adjustmtfangledeg = 0. # if consistency check fails, add a little angle
    bShowMTFDetail = False
    bShowCTF = True
    bIgnoreMTFError = False # in rare cases MTF is "ok" but looks wrong

    crLimit = 0.1 # contrast limit for MTF
    sigma_ext = 1.5 # gaussian blur factor for extreme finder

    def __init__(self):
        self.adjustmtfangledeg = 0.
        self.bShowMTFDetail = False
        self.bShowCTF = False
        self.bIgnoreMTFError = False

    def pixToGridScale_mm(self,cs):
        try:
            pixel_spacing_x = cs.dcmInfile.ImagerPixelSpacing[0] # PixelSpacing already corrects for magnification in DX!
            # DX
            sid = cs.dcmInfile.DistanceSourceToDetector # source to image distance
            sip = cs.dcmInfile.DistanceSourceToPatient  # source to patient (=table top?!)
            return pixel_spacing_x*sip/sid
        except:
            pixel_spacing_x = cs.dcmInfile.PixelSpacing[0] # PixelSpacing already corrects for magnification in DX!
            return  pixel_spacing_x

    def pix2phantomm(self, cs, pix):
        pix2phantommm = self.pixToGridScale_mm(cs)
        return pix*pix2phantommm

    def phantommm2pix(self, cs,mm):
        pix2phantommm = self.pixToGridScale_mm(cs)
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

    def TableOrWall(self,cs):
        # find out (based on SID) if Table or Wall stand is used, and which detector (could be wireless)
        detectorid = cs.dcmInfile.DetectorID # 0018,700A
        if cs.knownTableOrWall is not None:
            return cs.knownTableOrWall,detectorid

        if detectorid in cs.forceRoom.detector_name:
            cs.knownTableOrWall = cs.forceRoom.detector_name[detectorid]
            return cs.knownTableOrWall,detectorid
            
        result = ""
        if cs.forceRoom.sidtablemm < 0. and cs.forceRoom.sidwallmm > 0.:
            result =  lit.stWall
        elif cs.forceRoom.sidtablemm > 0. and cs.forceRoom.sidwallmm < 0.:
            result =  lit.stTable

        else:
            try:
                sid = cs.dcmInfile.DistanceSourceToDetector
            except:
                return lit.stUnknown

            if sid>1600.:
                result = lit.stWall
            else:
                result = lit.stTable
        if result != lit.stUnknown:
            cs.knownTableOrWall = result
        else:
            print "ERROR! Cannot determine if Wall or Table is used!"
        return result,detectorid

#----------------------------------------------------------------------
    def findPhantomOrientation(self,cs,ebbox=None):
        # Assume nicely centered image. Draw line from center to outwards:
        #  1 right: find 2 narrow (about 1 mm) peaks within 10 mm of each other. First = 90 mm on EAST
        #  2 left: find wide increase, but no 2 narrow (about 1 mm) peaks within 10 mm of each other. Must be WEST
        #  3 up: optionally find wide increase and afterwards 2 narrow (about 1 mm) peaks within 10 mm of each other. Must be North
        #  4 down: no wide increase and afterwards 2 narrow (about 1 mm) peaks within 10 mm of each other. Could be South
        error = True

        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        # 1.1 find approx center of phantom (screw)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]
        midx = int(.5*(widthpx-1)+.5)
        midy = int(.5*(heightpx-1)+.5)
        
        hlinepx = int(max(2,self.phantommm2pix(cs, 1))+.5) # don't expect a gridline to be larger than 2mm
        # make rois starting from the center, running to N,E,S, stopping 2cm short of edge
        seppx = self.phantommm2pix(cs, 20)# 2 cm away from edges
        blockheight = 10 # px
        rois = []
        rois.append([midx,midx+blockheight,midy,seppx]) # North
        rois.append([midx,widthpx-1-seppx,midy,midy+blockheight]) # East
        rois.append([midx,midx+blockheight,midy,heightpx-1-seppx]) # South
        lines = []
        edgepos = []
        for r in rois:
            stepx = 1 if r[1]>r[0] else -1
            stepy = 1 if r[3]>r[2] else -1
            if cs.mustbeinverted:
                smallimage = invertmax - cs.pixeldataIn[r[0]:r[1]:stepx,r[2]:r[3]:stepy]
            else:
                smallimage = cs.pixeldataIn[r[0]:r[1]:stepx,r[2]:r[3]:stepy]
            axis = 0 if np.shape(smallimage)[0] < np.shape(smallimage)[1] else 1
            
            # turn each roi into a average line
            line = np.average(smallimage,axis=axis)

            # find first position of drop
            threshold = (np.amax(line)+np.amin(line))/2
            threshold = (line[0]+np.amin(line))/2
            ep = -1
            for x in range(hlinepx,len(line)-hlinepx-1):
                if line[x-hlinepx]>threshold and line[x]<threshold and line[x+hlinepx]> threshold:
                    ep = x
                    break
            if ep == -1:
                print threshold,line
            edgepos.append(ep)
            lines.append(line)
        if cs.verbose:
            cs.hasmadeplots = True
            plt.figure()
            for line,lab in zip(lines,['N','E','S']):
                plt.plot(line,label=lab)
            plt.plot(edgepos,[threshold,threshold,threshold],'ro')
            print edgepos
            plt.legend()
        
        for p,k in zip(edgepos,['N','E','S']):
            print k,p,self.pix2phantomm(cs, p)   

        # the twice 90 mm distance is the northedge+southedge
        px90 = int( (edgepos[0]+edgepos[2])/2. +.5 )
        
        #adjust phantom center
        midypx = int(midy -(edgepos[0]-edgepos[2])/2.+.5)
        midxpx = int(midx +(edgepos[1]-px90))
        print midxpx,np.shape(cs.pixeldataIn)[0]/2.
        print midypx,np.shape(cs.pixeldataIn)[1]/2.

        # define cornerpoints of 90x90cm square
        roipts = [ 
            [midxpx-px90,midypx-px90], # ll
            [midxpx-px90,midypx+px90], # ul
            [midxpx+px90,midypx+px90], # ur
            [midxpx+px90,midypx-px90], # lr
        ]

        # true to properly locate the cornerpoints of the 90mm box
        rotangledeg = 0.
        error = False
        cs.bbox_confidence = 0.
        self.boxradmm = 90
        roipts = self.RotateBoundingBox(roipts,rotangledeg)
        error,cs.bbox_confidence = self.AlignROI(cs,roipts,"BoundingBox",True)
        for (xa,ya) in roipts:
            if xa<0 or ya<0:
                error = True
            elif xa>=widthpx or ya>=heightpx:
                error = True

        print "using rad,conf:",self.boxradmm,cs.bbox_confidence
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

    def AlignROI(self, cs, roipts, what,blast_attempt=True):
        error = True
        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        pix2phantommm = self.pixToGridScale_mm(cs)
        searchrad = int(2.5/pix2phantommm+.5) # 2 mm around location

        workimage = cs.pixeldataIn
        if 'MTF' in what: #
            searchrad = int(0.5/pix2phantommm+.5)

        searchrad = max(1,searchrad)
        print "%s searchrad="%what,searchrad
        widthpx = np.shape(workimage)[0] ## width/height in pixels
        heightpx = np.shape(workimage)[1]

        conf_pts = []
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
                    cropped = invertmax - workimage[minx:maxx+1,miny:maxy+1]
                else:
                    cropped = workimage[minx:maxx+1,miny:maxy+1]

                # smoothing to get rid of noise
                sigma = max(1.5,cropped.shape[0]/8.)
                cropped = scind.gaussian_filter(cropped, sigma)

                # find maximum location in projections
                xflat = np.mean(cropped,axis=0)
                yflat = np.mean(cropped,axis=1)
                if what == 'MTF_RXT02': # search for max
                    x1 = minx+np.unravel_index(yflat.argmax(), yflat.shape)[0]
                    y1 = miny+np.unravel_index(xflat.argmax(), xflat.shape)[0]
                else: # search for min
                    x1 = minx+np.unravel_index(yflat.argmin(), yflat.shape)[0]
                    y1 = miny+np.unravel_index(xflat.argmin(), xflat.shape)[0]
                if cs.verbose:
                    if 'MTF' in what:
                        plt.figure()
                        plt.imshow(cropped)
                        plt.title('Align '+str(kk)+ ' '+what+str(i))
                        print sigma,"shift ",i," of point ",kk," =(",x1-x0,",",y1-y0,")"
                        cs.hasmadeplots = True

                rp[0] = x1
                rp[1] = y1
            if not what == 'MTF_RXT02':
                roipts = self.ConsistencyAlign(cs,roipts0,roipts,what)
            conf_pts.append((self.ROIConfidence(cs,roipts,what),copy.deepcopy(roipts)))

        # Just take the last one; should be best!
        # conf_pts = sorted(conf_pts) #sorting orientation/distance box does worsen results
        if 'MTF' in what: # sorting MTF seems to help
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
            print "AlignRoi (",what,"):",label,", confidence too low:",confidence

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

        if what == "BoundingBox":
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

        if what == "MTF_38":
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

        if what == 'MTF_RXT02':
            l0mm = 23.32
            dx = roipts[0][0]-roipts[1][0]
            dy = roipts[0][1]-roipts[1][1]
            l1mm = self.pix2phantomm(cs,np.sqrt(dx*dx+dy*dy))

            return min(l1mm,l0mm)/max(l1mm,l0mm)
        
        # First the lengths of just the sides of the box
        pix2phantommm = self.pixToGridScale_mm(cs)
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

        if what=="BoundingBox":
            invertmax = 2**(self.invertmaxval(cs))-1
            """
            Calculate lengths and
		    // need to check only is NS^2 == len0^2+len1^2 and EW^2=len2^2+len3^2
		    // 1cm dev must lead to fail. So lengths (12cm) minus 11cm. Scaling is from (sid+10)/sid.
		    // find magnification
		    // we want a cm dev to drop confidence to 0.5, so all lengths (emperically) need to be reduced by magnified 11cm
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
            confidence = confidence **6.

        if what == 'MTF_38':
            #  colimit = .1 # just a scaling for confidence
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

        print what+"Confidence = ", (confidence*100.),"%"
        return confidence

    #----------------------------------------------------------------------
    def XRayField(self,cs,roipts=None):
        """
        Use either min across line between box and edge, or use corner values
        """
        workim = scind.gaussian_filter(cs.pixeldataIn, sigma=5)

        error = False
        if roipts is None:
            roipts = cs.po_roi

        cs.xrayNSWEmm = []
        # north- and southside
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'N',workim))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'S',workim))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'W',workim))
        cs.xrayNSWEmm.append(self.FindXRayEdge(cs,roipts,'E',workim))
        if min(cs.xrayNSWEmm)<1.:
            error = True
        else:
            error = False
        print 'Edge [N/S/W/E] cm = %.1f %.1f %.1f %.1f' % (cs.xrayNSWEmm[0]/10., cs.xrayNSWEmm[1]/10., cs.xrayNSWEmm[2]/10., cs.xrayNSWEmm[3]/10. )
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

    def FindXRayEdge(self,cs,roipts,side,workim=None):
        # travel from center to edges
        if workim is None:
            workim = cs.pixeldataIn
            
        invertmax = 2**(self.invertmaxval(cs))-1

        pix2phantommm = self.pixToGridScale_mm(cs)
        widthpx  = np.shape(workim)[0] ## width/height in pixels
        heightpx = np.shape(workim)[1]

        outvalue = cs.forceRoom.outvalue
        # for DiDi, just take the minimal corner value
        if outvalue<0:
            outvalue = min(workim[0][0], workim[-1][0],workim[0][-1],workim[-1][-1])
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
                        valvec.append(int(workim[xa,ya]))
                    break

                val00 = int(workim[x0,y0])
                val10 = int(workim[x1,y0])
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
                if v<threshLow:
                    found = True
                    edgemm.append( p*pix2phantommm+self.boxradmm )
                    if cs.verbose:
                        plt.plot(p,v,'bo')
                        cs.hasmadeplots = True
                    break

            if not found and threshHigh>threshLow:
                for p,v in zip(posvec,valvec):
                    if v>=threshHigh:
                        found = True
                        edgemm.append( p*pix2phantommm+self.boxradmm )
                        if cs.verbose:
                            plt.plot(p,v,'bo')
                            cs.hasmadeplots = True
                        break
                
            if not found:
                edgemm.append(0)
        
        return max(edgemm)

#----------------------------------------------------------------------
    def CuWedge(self,cs, roipts_orig=None):
        # 1. Make box around wedge
        # 2. Find all steps and do some statistcs
        error = True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)
        bump = 2.
        xmmll = -63.  +bump
        ymmur =  65.  +bump
        xmmur =  62.  -bump
        ymmll =  82.  -bump

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
            print "[CuWedge]: Error, phantom angle too large, cannot sample wedge"
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

        digi_13 = [2.30, 1.85, 1.40, 1.00, 0.65, 0.30, 0.00]   # mm Cu; according to DIN 6866-13; note this thickness includes backplate
        #    1. Cut out rectangular roi
        xmin = roipts_orig[0][0]
        xmax = roipts_orig[1][0]
        ymin = roipts_orig[2][1]
        ymax = roipts_orig[1][1]
        if ymin>ymax:
            print "[AnalyseWedge]: Error, phantom angle too large, cannot sample wedge"
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
            cs.hasmadeplots = True

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
            cs.hasmadeplots = True

        # 2.2 Find the edges between the steps
        n_edges = 6
        posedges = []
        flatpix = int( self.phantommm2pix(cs,3.0)+.5)
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
            print "Edge 0 at ",posedges[0]
            for ix in range(1,n_edges):
                print "Edge ",ix," at ",posedges[ix]," sep= ",posedges[ix]-posedges[ix-1]

        # 2.3 Calculate statistics for each step
        cs.cuwedge.roi_mean = []
        cs.cuwedge.roi_sdev = []

        xlo = 0
        ylo = 0   # flatpix
        yhi = hei # - flatpix
        if ylo>yhi:
            print "[AnalyseWedge]: Error, phantom angle too large, cannot sample wedge"
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
            cs.cuwedge.roi_mmcu.append( digi_13[ix] )

        cs.cuwedge.dynamicRange = max(cs.cuwedge.roi_mean[n_edges]/cs.cuwedge.roi_mean[0],cs.cuwedge.roi_mean[0]/cs.cuwedge.roi_mean[n_edges])

        if cs.verbose:
            print "mmCu","SNR","CNR"
            for m,s,c in zip(cs.cuwedge.roi_mmcu,cs.cuwedge.roi_snr,cs.cuwedge.roi_cnr):
                print m,s,c

            #cu_cs.guesskVp = guesskVp
        """
	    guesskVp = EstimatekV(roi_mean, roi_mmcu,give_info,cs);

        """
        error = False

        return error
    #----------------------------------------------------------------------
    def MTF(self,cs,roipts_orig=None):
        error,confid = True,0.
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        if cs.forceRoom.linepairmodel == 'RXT02':
            x06px,y06px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy06mm[0],cs.forceRoom.xy06mm[1])
            x10px,y10px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy10mm[0],cs.forceRoom.xy10mm[1])
            roipts = [ 
                [x06px,y06px],
                [x10px,y10px]
            ]
        elif cs.forceRoom.linepairmodel == '38':
            x18px,y18px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy18mm[0],cs.forceRoom.xy18mm[1])
            x06px,y06px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy06mm[0],cs.forceRoom.xy06mm[1])
            x46px,y46px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy46mm[0],cs.forceRoom.xy46mm[1])
            x14px,y14px = self.phantomposmm2pix(roipts_orig,cs.forceRoom.xy14mm[0],cs.forceRoom.xy14mm[1])
            roipts = [ 
                [x18px,y18px],
                [x06px,y06px],
                [x14px,y14px],
                [x46px,y46px] ]
        else:
            raise ValueError('[MTF] Unknown linepairmodel %s'%cs.forceRoom.linepairmodel)
        
        error, confid = self.AlignROI(cs,roipts,'MTF_'+cs.forceRoom.linepairmodel)

        cs.mtf.roi = roipts
        cs.mtf.pos_confidence = confid

        if not error:
            error = self.AnalyseMTF(cs,roipts)
        return error

    def AnalyseMTF(self,cs,roipts_orig=None):
        if cs.forceRoom.linepairmodel == 'RXT02':
            return self._AnalyseMTF_RXT02(cs,roipts_orig)
        else:
            return self._AnalyseMTF_38(cs,roipts_orig)
            
    def _AnalyseMTF_38(self,cs,roipts_orig=None):
        """
        For linepairs insert model 38.
        Find rotation angle of linepairs insert. Cut out the back rotated insert.
        For each linepair, do an analysis
        """
        error = True
        invertmax = 2**(self.invertmaxval(cs))-1
    
        if roipts_orig is None:
            roipts_orig = cs.mtf.roi
    
        extend18 = self.phantommm2pix(cs,2.8)  # extend bbox beyond dot in '1.8' [mm]
        extend46 = self.phantommm2pix(cs,3.2)  # extend beyond dot in '4.6' [mm]
        print "2.8",extend18
        print "3.2",extend46
        # First cut out rotated roi
        id18 = 0
        id06 = 1
        id14 = 2
        id46 = 3
        len1846 = np.sqrt((roipts_orig[id18][0]-roipts_orig[id46][0])**2+(roipts_orig[id18][1]-roipts_orig[id46][1])**2)
        print "1846=",len1846
        extend18 = 0.0786*len1846
        extend46 = 0.0898*len1846
        copyimage = cs.pixeldataIn.astype(float)
        rotanglerad = 3.*np.pi/2.-.5*(np.arctan2((roipts_orig[id18][1]-roipts_orig[id46][1]),(roipts_orig[id18][0]-roipts_orig[id46][0]))+np.arctan2((roipts_orig[id06][1]-roipts_orig[id14][1]),(roipts_orig[id06][0]-roipts_orig[id14][0])))
        rotanglerad += self.adjustmtfangledeg/180*np.pi
        rotangledeg = (rotanglerad/np.pi*180.)
        print "MTF at",rotangledeg, "degrees"
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
            smallimage = invertmax-rotimage[minxco:maxxco+1,minyco:maxyco+1]
        else:
            smallimage = rotimage[minxco:maxxco+1,minyco:maxyco+1]
        cs.lastimage = smallimage
        
        for rp in roipts:
            rp[0] -= minxco
            rp[1] -= minyco
    
        return self._MTF_smallimage(cs, smallimage)
        
    def _AnalyseMTF_RXT02(self,cs,roipts_orig=None):
        """
        For linepairs insert model RXT02.
        Find rotation angle of linepairs insert. Cut out the back rotated insert.
        For each linepair, do an analysis
        """
        error = True
        invertmax = 2**(self.invertmaxval(cs))-1

        if roipts_orig is None:
            roipts_orig = cs.mtf.roi
        print roipts_orig

        extendup = self.phantommm2pix(cs,1.5)  # extend bbox beyond dot in '1.8' [mm]
        extenddo = self.phantommm2pix(cs,15.)  # extend beyond dot in '4.6' [mm]

        # extend starting points
        dy = roipts_orig[1][1]-roipts_orig[0][1]
        dx = roipts_orig[1][0]-roipts_orig[0][0]

        l0 = np.sqrt(dx*dx+dy*dy)
        roipts_orig[0][0] -= extendup/l0*dx
        roipts_orig[0][1] -= extendup/l0*dy
        roipts_orig[1][0] += extenddo/l0*dx
        roipts_orig[1][1] += extenddo/l0*dy
        
        #append points
        l1 = self.phantommm2pix(cs,35.)
        px = roipts_orig[1][0] - l1/l0*dy
        py = roipts_orig[1][1] + l1/l0*dx
        roipts_orig.append([px,py])

        px = roipts_orig[0][0] - l1/l0*dy
        py = roipts_orig[0][1] + l1/l0*dx
        roipts_orig.append([px,py])
        
        cs.mtf.roi = roipts_orig
        
        # angle
        rotanglerad = np.arctan2(dx,dy)

        # First cut out rotated roi
        copyimage = cs.pixeldataIn.astype(float)
        rotangledeg = (rotanglerad/np.pi*180.)
        print "MTF at",rotangledeg, "degrees"
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
            smallimage = invertmax-rotimage[minxco:maxxco+1,minyco:maxyco+1]
        else:
            smallimage = rotimage[minxco:maxxco+1,minyco:maxyco+1]
        cs.lastimage = smallimage

        for rp in roipts:
            rp[0] -= minxco
            rp[1] -= minyco

        return self._MTF_smallimage(cs, smallimage)

    def _MTF_smallimage(self,cs,smallimage):
        """
        Analysis driver of image cropped to linepairs insert
        """
        error = True
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

        start_endpos = []
        start_endpos = self.FillMTFBarDetails(cs,smallimage)

        # Find high contrast response of line patterns
        for vpi in range(1,num_freq):
            contrast_response[vpi] = 0.
            contrast_high[vpi]     = 0.
            contrast_low[vpi]      = 0.
            contrast_tops[vpi]     = 0
            contrast_bots[vpi]     = 0
            calc_freq[vpi]         = 0.
            contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = self.AnalyseMTF_Part(cs,smallimage, start_endpos[vpi-1], vpi)

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
            print "Error in MTF: Rotated image?"
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
        print "maxid:",maxid
        if maxid<5 and not self.bIgnoreMTFError:
            print "Error in MTF: Rotated image?"
            return error
        slope, intercept, r_value, p_value, std_err = stats.linregress(contrast_freqs[0:maxid],calc_freq[0:maxid])
        if r_value**2<0.7:
            print "maxid:",maxid
            for co,ca in zip(contrast_freqs,calc_freq):
                print co,ca
        # To get coefficient of determination (r_squared)
#        print "slope:",slope
#        print "intercept:",intercept
#        print "maxid:",maxid
#        print "r-squared:", r_value**2
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
        print "mtf_freq_confidence:",mtf_freq_confidence
        if mtf_freq_confidence<.7:
            print "found/first_error/needfound:",found,first_error,needfound
            print "slope/r2:",slope,r_value**2
            if cs.verbose:
                plt.figure()
                plt.plot(contrast_freqs[0:maxid],calc_freq[0:maxid],'bo')
                plt.title("found vs given freq")
                cs.hasmadeplots = True

    #       print "confid:",mtf_found_confidence
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
        startendpos = [] # [ [x0,y0], [x1,y1] ]
        endpos = []

        wid = smallimage.shape[0]
        hei = smallimage.shape[1]

        if cs.forceRoom.linepairmodel == 'RXT02':
            washeightLo = 287
            waswidthLo  = 253
    
            # fill positions with manual values; x is for later
            startendpos.append([[0,  1],[]]) # absolute 0.6
            startendpos.append([[0, 45],[]]) # absolute 0.7
            startendpos.append([[0, 85],[]]) # absolute 0.8
            startendpos.append([[0,121],[]]) # absolute 0.9
            startendpos.append([[0,156],[]]) # absolute 1.0
            startendpos.append([[0,188],[]]) # absolute 1.2
            startendpos.append([[0,216],[]]) # absolute 1.4
            startendpos.append([[0,244],[0,268]]) # absolute 1.6
    
            startendpos.append([[0,  5],[]]) # absolute 1.8
            startendpos.append([[0, 31],[]]) # absolute 2.0
            startendpos.append([[0, 54],[]]) # absolute 2.2
            startendpos.append([[0, 76],[0, 98]]) # absolute 2.5
            startendpos.append([[0,120],[]]) # absolute 2.8
            startendpos.append([[0,140],[]]) # absolute 3.1
            startendpos.append([[0,160],[]]) # absolute 3.4
            startendpos.append([[0,181],[0,200]]) # absolute 3.7
    
            startendpos.append([[0,214],[]]) # absolute 4.0
            startendpos.append([[0,232],[]]) # absolute 4.3
            startendpos.append([[0,252],[]]) # absolute 4.6
            startendpos.append([[0,270],[0,286]]) # absolute 5.0 

            lo_xpos = [160,160+50]
            hi_xpos = [45,45+50]
        
        elif cs.forceRoom.linepairmodel == '38'  :
            washeightLo = 302
            waswidthLo  = 300
    
            # fill positions with manual values; x is for later
            startendpos.append([[0,  6],[]]) # absolute 0.6
            startendpos.append([[0, 48],[]]) # absolute 0.7
            startendpos.append([[0, 86],[]]) # absolute 0.8
            startendpos.append([[0,122],[]]) # absolute 0.9
            startendpos.append([[0,160],[]]) # absolute 1.0
            startendpos.append([[0,197],[]]) # absolute 1.2
            startendpos.append([[0,225],[]]) # absolute 1.4
            startendpos.append([[0,252],[0,273]]) # absolute 1.6
    
            startendpos.append([[0,  4],[]]) # absolute 1.8
            startendpos.append([[0, 28],[]]) # absolute 2.0
            startendpos.append([[0, 54],[]]) # absolute 2.2
            startendpos.append([[0, 79],[0, 102]]) # absolute 2.5
            startendpos.append([[0,125],[]]) # absolute 2.8
            startendpos.append([[0,144],[]]) # absolute 3.1
            startendpos.append([[0,163],[]]) # absolute 3.4
            startendpos.append([[0,182],[0,202]]) # absolute 3.7
    
            startendpos.append([[0,218],[]]) # absolute 4.0
            startendpos.append([[0,237],[]]) # absolute 4.3
            startendpos.append([[0,257],[]]) # absolute 4.6
            startendpos.append([[0,276],[0,295]]) # absolute 5.0 
            
            lo_xpos = [185,185+50]
            hi_xpos = [60,60+50]
            
        for i in range(len(startendpos)-1):
            if len(startendpos[i][1]) == 0:
                startendpos[i][1] = list(startendpos[i+1][0])

        # start with low frequencies
        ipbar = smallimage[int(.75*wid):int(.85*wid),0:hei]
        if cs.verbose:
            plt.figure()
            plt.imshow(ipbar)
            plt.title('low frequencies')
            cs.hasmadeplots = True

        # low freqs
        xlo0 = int(.5+1.*lo_xpos[0]/waswidthLo*wid)
        xlo1 = int(.5+1.*lo_xpos[1]/waswidthLo*wid)
        # high freqs
        xhi0 = int(.5+1.*hi_xpos[0]/waswidthLo*wid)
        xhi1 = int(.5+1.*hi_xpos[1]/waswidthLo*wid)
        for vpi in range(len(startendpos)):
            start = startendpos[vpi][0]
            end  = startendpos[vpi][1]
            if vpi < 8:
                start[0] = xlo0
                end[0]   = xlo1
            else:
                start[0] = xhi0
                end[0]   = xhi1
            start[1] = max(0,int(1.*hei/washeightLo*start[1]+.5))
            end[1]   = min(hei-1,int(1.*hei/washeightLo*end[1]+.5))
            startendpos[vpi] = [start,end]


        return startendpos

    def AnalyseMTF_Part(self,cs,smallimage, start_endpos, vpi):
        """
        Determine contrast in one bar pattern element
            contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = self.AnalyseMTF_Part(smallimage, startpos[vpi-1],endpos[vpi-1], vpi)
        Already inverted smallimage as input
        """
        startpos = start_endpos[0]
        endpos   = start_endpos[1]
        
        ipbar = smallimage[startpos[0]:endpos[0],startpos[1]:endpos[1]]
        print '[AnalyseMTF_Part]',vpi,startpos,endpos
#        cs.hasmadeplots = True
#        plt.figure()
#        plt.title(vpi)
#        plt.imshow(ipbar) 
        
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
            print "[AnalyseMTF_Part] SKIP: no pattern left"
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
                print "[AnalyseMTF_Part] ERROR: cannot find baseline"
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
                cs.hasmadeplots = True

        if (mustplot == True or (len(tops)!=3 or len(bots)!=2) ) and cs.verbose:
            print "vpi=",vpi," length(pattern)=",len(pattern)
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
            print "Found",len(tops)," tops and",len(bots)," bottoms. Contrast=",contrast_response
#        print vpi,contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq
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
                print "Too many extrema; rerun with larger sigma"
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
                print "Too few extrema; rerun with smaller sigma"
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
            print "[FindExtrema]C ntops/nbots = ",len(tops),"/",len(bots)
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
                    print "[FindAllExtrema] ntops>3! Using only first 3."
                    mustplot = True
            if(xderiv1[y]<=0. and xderiv1[y+1]>0.):
                if(pattern[y+1]<pattern[y]):
                    if(len(tops)>0):
                        bots.append(y+1)
                else:
                    if(len(tops)>0):
                        bots.append(y)
                if(len(bots)>2):
                    print "[FindAllExtrema] nbots>2! Using only first 2."
                    mustplot = True

        if mustplot: # SOMETHING WRONG, INGORE XDERIV AND JUST LOOK AT MIN/MAX
            tops = []
            bots = []
            mustplot = False

            for y in range(ymin+1,ymax):
                if(pattern[y]>pattern[y+1] and pattern[y]>pattern[y-1]):
                    tops.append(y)
                    if(len(tops)>3):
                        print "[FindAllExtrema] ntops2>3! Using only first 3."
                        mustplot = True
                if(pattern[y]<pattern[y+1] and pattern[y]<pattern[y-1]):
                    if(len(tops)>0):
                        bots.append(y)
                    if(len(bots)>2):
                        print "[FindAllExtrema] nbots2>2! Using only first 2."
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
        #        print poly
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
                    print freq[i],i,j,factor*freq[i],fnyq
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
        error = True
        if roipts_orig is None:
            roipts_orig = cs.po_roi

        invertmax = 2**(self.invertmaxval(cs))-1
        low_cnr = []
        mean_sg = []
        mean_bk = []
        sdev_sg = []
        sdev_bk = []

        # discs are 1 cm in diameter
        diampx = self.phantommm2pix(cs, 6.)
        radpx  = int(diampx/2.+.5)
        yrefmm = 10. # offset to background

        cs.loco.lo_rois = [] # [x0,y0,rad]
        cs.loco.lo_rois_bku = [] # [x0,y0,rad] # rois above for bk
        cs.loco.lo_rois_bkd = [] # [x0,y0,rad] # rois below for bk
        cs.loco.mean_sg = [] # average in object
        cs.loco.sdev_sg = [] # average in object
        
        dxmm = 17. # dist in mm between loco objects
        dypx = self.phantommm2pix(cs,11) # dist in mm between loco objects
        for i in range(6):
            px,py = self.phantomposmm2pix(roipts_orig,-42.5+i*dxmm,-73.)
            cs.loco.lo_rois.append([px,py,1.*radpx])
            cs.loco.lo_rois_bku.append([px,py-dypx,1.*radpx])
            cs.loco.lo_rois_bkd.append([px,py+dypx,1.*radpx])

        for r in cs.loco.lo_rois:
            [px,py,radpx] = r
            # select data inside circle only 
            if cs.mustbeinverted:
                smallimage = invertmax - cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            else:
                smallimage = cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            wid = 2*radpx
            x,y = np.indices((wid, wid))
            mid = wid/2.
            mask = ((x-mid)**2 + (y-mid)**2 ) < radpx**2
            cs.loco.mean_sg.append(np.mean(smallimage[mask]))
            cs.loco.sdev_sg.append(np.std(smallimage[mask]))

        mean_bku = []
        mean_bkd = []
        sdev_bku = []
        sdev_bkd = []
        for r in cs.loco.lo_rois_bku:
            [px,py,radpx] = r
            # select data inside circle only 
            if cs.mustbeinverted:
                smallimage = invertmax - cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            else:
                smallimage = cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            wid = 2*radpx
            x,y = np.indices((wid, wid))
            mid = wid/2.
            mask = ((x-mid)**2 + (y-mid)**2 ) < radpx**2
            mean_bku.append(np.mean(smallimage[mask]))
            sdev_bku.append(np.std(smallimage[mask]))

        for r in cs.loco.lo_rois_bkd:
            [px,py,radpx] = r
            # select data inside circle only 
            if cs.mustbeinverted:
                smallimage = invertmax - cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            else:
                smallimage = cs.pixeldataIn[px-radpx:px+radpx,py-radpx:py+radpx]
            wid = 2*radpx
            x,y = np.indices((wid, wid))
            mid = wid/2.
            mask = ((x-mid)**2 + (y-mid)**2 ) < radpx**2
            mean_bkd.append(np.mean(smallimage[mask]))
            sdev_bkd.append(np.std(smallimage[mask]))

        cs.loco.low_cnr = []
        cs.loco.mean_bk = []
        cs.loco.sdev_bk = []
        # calculate for each object the max contrast to noise wrt upper or lower background. Keep the values for the max cnr
        for i in range(len(cs.loco.lo_rois_bkd)):
            cnr_up = (mean_bku[i]-cs.loco.mean_sg[i])/np.sqrt(0.5*(cs.loco.sdev_sg[i]**2+sdev_bku[i]**2))
            cnr_do = (mean_bkd[i]-cs.loco.mean_sg[i])/np.sqrt(0.5*(cs.loco.sdev_sg[i]**2+sdev_bkd[i]**2))
            if cnr_up>cnr_do:
                cs.loco.low_cnr.append(cnr_up)
                cs.loco.mean_bk.append(mean_bku[i])
                cs.loco.sdev_bk.append(sdev_bku[i])
            else:
                cs.loco.low_cnr.append(cnr_do)
                cs.loco.mean_bk.append(mean_bkd[i])
                cs.loco.sdev_bk.append(sdev_bkd[i])

            if cs.verbose:
                print "mean fg/bk=",cs.loco.mean_sg[i],"/",cs.loco.mean_bk[i]
                print "sdev fg/bk=",cs.loco.sdev_sg[i],"/",cs.loco.sdev_bk[i]
                print 'cnr = ',cs.loco.low_cnr[i]

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
                {'key':"0008,0080",  'name':"InstitutionName"},
                {'key':"0008,1010",  'name':"StationName"},
                {'key':"0008,1030",  'name':"StudyDescription"},
                {'key':"0008,103E",  'name':"SeriesDescription"},
                {'key':"0008,1070",  'name':"Operator's Name"},
                {'key':"0008,1090",  'name':"Modelname"},
                {'key':"0010,0020",  'name':"PatientID"},
                {'key':"0018,0015",  'name':"BodyPartExamined"},
                {'key':"0018,0060",  'name':"kVp"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1030",  'name':"ProtocolName"},
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)"},
                {'key':"0018,1111",  'name':"DistanceSourceToPatient (mm)"},
                {'key':"0018,1153",  'name':"Exposure (uAs)"},
                {'key':"0018,115E",  'name':"ImageAreaDoseProduct"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"0018,1166",  'name':"Grid"},
                {'key':"0018,1405",  'name':"RelativeXRayExposure"},
                {'key':"0018,1702",  'name':"CollimatorLeft"},
                {'key':"0018,1704",  'name':"CollimatorRight"},
                {'key':"0018,1706",  'name':"CollimatorUp"},
                {'key':"0018,1708",  'name':"CollimatorDown"},
                {'key':"0018,7001",  'name':"DetectorTemperature"},
                {'key':"0018,700A",  'name':"DetectorID"},
                {'key':"0018,700C",  'name':"DateCalibration"},
                {'key':"0018,7050",  'name':"FilterMaterial"},
                {'key':"0018,7062",  'name':"ExposureControlMode"},
                {'key':"0018,8150",  'name':"ExposureTime (us)"},
                {'key':"0028,0101",  'name':"BitsStored"},
                {'key':"0040,8302",  'name':"EntranceDose_mGy"},
                {'key':"200B,70BA",  'name':"FocalSpot"},
                {'key':"200B,1028",  'name':"PrivValue1"},
                {'key':"200B,7096",  'name':"PrivStand"},
            ]


        elif(info == "qcwad"):
            offset = -25 # rank must be negative, so recalc as offset+real position
            dicomfields = [
                {'key':"0008,0021",  'name':"SeriesDate"},
                {'key':"0008,0031",  'name':"SeriesTime", 'quantity':'time', 'level':1, 'rank':offset+2}, # spot 1 reserved for stand
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)", 'quantity':'distance', 'level':1, 'rank':offset+3},
                {'key':"0018,0060",  'name':"kVp", 'level':1, 'rank':offset+4},
                {'key':"0018,7050",  'name':"FilterMaterial", 'quantity':'filter', 'level':1, 'rank':offset+5},
                {'key':"0018,1166",  'name':"Grid", 'quantity':'grid', 'level':1, 'rank':offset+6},
                {'key':"200B,70BA",  'name':"FocalSpot", 'quantity':'spot', 'level':1, 'rank':offset+7},
                {'key':"0018,1702",  'name':"CollimatorLeft",'quantity':'Left','level':1, 'rank':offset+8},
                {'key':"0018,1704",  'name':"CollimatorRight",'quantity':'Right','level':1, 'rank':offset+9},
                {'key':"0018,1706",  'name':"CollimatorUp",'quantity':'Up','level':1, 'rank':offset+10},
                {'key':"0018,1708",  'name':"CollimatorDown",'quantity':'Down','level':1, 'rank':offset+11},

                {'key':"0018,8150",  'name':"ExposureTime (us)", 'quantity':'us','level':1,'rank':offset+12},
                {'key':"0018,1153",  'name':"Exposure (uAs)", 'quantity':'uAs','level':1,'rank':offset+13},
                {'key':"0018,115E",  'name':"ImageAreaDoseProduct", 'quantity':'DAP','level':1,'rank':offset+14},
            
                {'key':"0008,103E",  'name':"SeriesDescription"},
                {'key':"0008,1010",  'name':"StationName"},
                {'key':"0008,1070",  'name':"Operator's Name"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,700A",  'name':"DetectorID"},
                {'key':"0018,700C",  'name':"DateCalibration"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"200B,7063",  'name':"Postprocessing", 'quantity':'processing'},
                {'key':"0018,1405",  'name':"RelativeXRayExposure"},
                {'key':"0040,8302",  'name':"EntranceDose_mGy"},
                {'key':"200B,1028",  'name':"PrivValue1"},
                {'key':"200B,7096",  'name':"PrivStand"},
                {'key':"0018,0015",  'name':"BodyPartExamined"},
                {'key':"0018,1030",  'name':"ProtocolName"},
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
    def QCUnif(self,cs):
        """
        Separate uniformity in one step!
        """

        print '[QCUnif]',cs.dcmInfile.SeriesDescription
        error,msg = self.Uniformity(cs)

        if error:
            msg += "Uniformity (try)) "
            print  msg
            return error,msg
            
        if cs.unif.expertmode:
            exproi=cs.unif.expert_roipts
            print '[Uniformity] Cropped to',exproi

        return error,msg
            

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
        print '[QCNormi13]',cs.dcmInfile.SeriesDescription
        self.CropNormi13(cs)
        error,msg = self.checkPhantomRotation(cs)
        ###
        if error:
            msg += "AlignROI(BoundingBox) "

        # 2: travel straight along NS and find edge of x-ray; similar for EW
        if not error:
            error = self.XRayField(cs)
            if error:
                msg += "XRayField "

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
        stand,detector = self.TableOrWall(cs)

        if not cs.unif is None:
            num_art = 0
            if cs.unif.art_clusters != None:
                num_art = len(cs.unif.art_clusters)

            offset = -25 # rank must be negative, so recalc as offset+real position

            labvals.append( {'name':'Uniformity_(%)','value':cs.unif.unif_pct, 'quantity':'Uniformity','level':1,'rank':offset+7} )
            labvals.append( {'name':'Artefacts','value':num_art, 'quantity':'Artefacts','level':1,'rank':offset+6} )


            for kk in range(0,len(cs.unif.means)):
                labvals.append( {'name':'avg_'+str(kk),'value':cs.unif.means[kk], 'quantity':'Mean','level':2} )
            for kk in range(0,len(cs.unif.stdevs)):
                labvals.append( {'name':'sd_'+str(kk),'value':cs.unif.stdevs[kk], 'quantity':'STDev','level':2} )
            for lab,val in zip(['xmin','xmax','ymin','ymax'],cs.unif.art_borderpx):
                labvals.append( {'name':'crop_%s'%lab,'value':val, 'quantity':'borderpx','level':2} )
   
            return labvals

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

        """
        ## uniformity
        if stand == lit.stWall:
            label = 'LRUniformity'
            value = 100.*cs.unif.LRuniformity
        else:
            label = 'ROIUniformity'
            value = 100.*cs.unif.ROIuniformity
        labvals.append( {'name':label,'value':value, 'quantity':'Uniformity','level':1,'rank':offset+13} )
        """
        
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

    def drawThickCircle(self,draw,x,y,rad,color,thick):
        for t in range(-(thick-1)/2,(thick+1)/2):
            r1 = rad+t
            draw.ellipse((x-r1,y-r1,x+r1,y+r1), outline=color)

    def drawThickRectangle(self,draw,roi,color,thick):
        #[(x0,y0),(x1,y1)]
        if len(roi) == 4:
            x0 = roi[0]
            y0 = roi[1]
            x1 = roi[2]
            y1 = roi[3]
        else:
            x0 = roi[0][0]
            y0 = roi[0][1]
            x1 = roi[1][0]
            y1 = roi[1][1]

        for t in range(-(thick-1)/2,(thick+1)/2):
            draw.rectangle([(x0+t,y0+t),(x1-t,y1-t)],outline=color)

    def saveAnnotatedImage(self,cs,fname):
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        rectrois = []
        polyrois = []
        circlerois = []

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        if cs.unif is None:
            # first the base image
            im = scipy.misc.toimage(cs.pixeldataIn.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            # now add all box rois
            if len(cs.po_roi) >0:
                polyrois.append(cs.po_roi) # phantom orientation box
            if len(cs.xr_roi) >0:
                polyrois.append(cs.xr_roi) # xray edges
            if len(cs.cuwedge.roi) >0:
                polyrois.append(cs.cuwedge.roi) # Cu wedge box
            for r in cs.cuwedge.step_rois:
                polyrois.append(r) # Cu wedge element
            if len(cs.mtf.roi) >0:
                polyrois.append(cs.mtf.roi) # MTF box

            # add circlerois
            for r in cs.loco.lo_rois: # low contrast elements
                circlerois.append(r)
            for r in cs.loco.lo_rois_bku:
                circlerois.append(r)
            for r in cs.loco.lo_rois_bkd:
                circlerois.append(r)
            
        else:
            # first find the offset of the boxes that were defined on uncropped image
            if cs.unif.expertmode:
                exproi=cs.unif.expert_roipts
            else:
                bpx = cs.forceRoom.artefactborderpx
                exproi=[bpx[0],np.shape(cs.pixeldataIn)[0]-1-bpx[1],bpx[2],np.shape(cs.pixeldataIn)[1]-1-bpx[3]]

            wid,hei = np.shape(cs.unif.art_image)
            # first copy the artimage into the base image
            pdCopy = np.zeros(np.shape(cs.pixeldataIn))
            pdCopy[exproi[0]:exproi[0]+wid,exproi[2]:exproi[2]+hei] = cs.unif.art_image
            im = scipy.misc.toimage(pdCopy.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            # now add all box rois

            # show the box of the whole image so we can see what has been cropped
            wid,hei = np.shape(cs.pixeldataIn)
            #ectrois.append([(-exproi[0],-exproi[2]),(-exproi[0]+wid,-exproi[2]+hei)]) #[(x0,y0),(x1,y1)]
            rectrois.append([(0,0),(wid-1,hei-1)]) #[(x0,y0),(x1,y1)]

            # uniformity boxes
            for r in cs.unif.unif_rois: #[x0,dx,y0,dy]
                #rectrois.append([(r[0]-exproi[0],r[2]-exproi[2]),(r[0]-exproi[0]+r[1],r[2]-exproi[2]+r[3])])#[(x0,y0),(x1,y1)]
                rectrois.append([(r[0],r[2]),(r[0]+r[1],r[2]+r[3])])#[(x0,y0),(x1,y1)]

            # artefact circles
            #  translate(exproi[0],exproi[2])
            for r in cs.unif.art_rois: # x,y,r
                circlerois.append([r[0]+exproi[0],r[1]+exproi[2],r[2]])


        # now draw all rois in reserved color
        draw = ImageDraw.Draw(im)
        for r in polyrois:
            roi =[]
            for x,y in r:
                roi.append( (int(x+.5),int(y+.5)))
            draw.polygon(roi,outline=0)
        
        for r in rectrois:
            #draw.rectangle(r,outline=0)
            self.drawThickRectangle(draw, r, 0, 3)
            
        # now draw all cirlerois in reserved color
        for x,y,r in circlerois: # low contrast elements
            draw.ellipse((x-r,y-r,x+r,y+r), outline=0)
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
        2. Find most likely location of CuWedge (thick)
        3. Use location wrt center to find and undo 90 degree rotations
        """
        error = True
        msg = ""
        
        # inversion
        invertmax = 2**(self.invertmaxval(cs))-1

        # 2. Make box at x-5cm and calculate avg
        seppx = self.phantommm2pix(cs, 20)# 2 cm away from edges
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]

        if cs.mustbeinverted:
            smallimage = invertmax - cs.pixeldataIn[seppx:widthpx-seppx:3,seppx:heightpx-seppx:3] #every third pixel
        else:
            smallimage = invertmax - cs.pixeldataIn[seppx:widthpx-seppx:3,seppx:heightpx-seppx:3] #every third pixel
        
        smallimage = scind.gaussian_filter(smallimage, sigma=5)
        # Find location of min transmission (max attenuation)    
        x0,y0 = np.unravel_index(smallimage.argmin(), smallimage.shape)
        immidx = smallimage.shape[0]/2
        immidy = smallimage.shape[1]/2
        cs.lastimage = smallimage
        print x0,y0,immidx,immidy
        ang = 0
        msg = ""
        if x0 < immidx: # correct side on x
            if y0 < immidy: # correct side on y
                pass
            else:
                ang = 3
        else:
            if y0 < immidy: # +90rot
                ang = 1
            else:
                ang = 2

        cs.po_rot = ang
        if ang>0:
            cs.pixeldataIn = np.rot90(cs.pixeldataIn,-ang)
            cs.po_rot = 90*ang

        error = self.findPhantomOrientation(cs)
        return error,msg

    def CropNormi13(self,cs):
        cs_unif = QCUniformity_lib.UnifStruct(cs.dcmInfile, cs.pixeldataIn)
        cs_unif.mustbeinverted = cs.mustbeinverted
        qc_unif = QCUniformity_lib.Uniformity_QC()

        if qc_unif.NeedsCropping2(cs_unif,mode='normi13'):
            print '*** Needs Cropping ***'
            
            widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
            heightpx = np.shape(cs.pixeldataIn)[1]

            bpx = 0#self.phantommm2pix(cs,10)
            qc_unif.RestrictROI(cs_unif)
            [xmin_px,xmax_px, ymin_px,ymax_px] = cs_unif.expert_roipts
            xmin_px = max(0,xmin_px-bpx)
            ymin_px = max(0,ymin_px-bpx)
            xmax_px = min(widthpx-1,xmax_px+bpx)
            ymax_px = min(heightpx-1,ymax_px+bpx)
            cs.pixeldataIn = cs.pixeldataIn[xmin_px:xmax_px+1,ymin_px:ymax_px+1]
            
    def Uniformity(self,cs):
        # run external QCUniformity
        usestructure = True
        
        cs.unif = QCUniformity_lib.UnifStruct(cs.dcmInfile, cs.pixeldataIn)
        cs.unif.mustbeinverted = cs.mustbeinverted
        qc_unif = QCUniformity_lib.Uniformity_QC()

        bpx = cs.forceRoom.artefactborderpx
        cropped = False
        wid,hei = np.shape(cs.pixeldataIn)
        if qc_unif.NeedsCropping2(cs.unif,mode='uniformity'):
            roipts = qc_unif.RestrictROI(cs.unif)

            """
            minx = min(r[0] for r in roipts)
            miny = min(r[1] for r in roipts)
            maxx = max(r[0] for r in roipts)
            maxy = max(r[1] for r in roipts)
            
            if minx > 1: 
                bpx[0] += 10 
            if miny > 1: 
                bpx[2] += 10 
            if maxx < wid-1: 
                bpx[1] += 10 
            if maxy < hei-1: 
                bpx[3] += 10 
            """
            cs.unif.expertmode = True
            cropped = True
        if usestructure:
            qc_unif.artefactDetectorParameters(UseStructure=False, bkscale=25, fgscale=5.0, threshold=15)
        else:
            qc_unif.artefactDetectorParameters(UseStructure=True, bkscale=25, fgscale=5.0, threshold=3000)

        error = qc_unif.Uniformity(cs.unif,cs.forceRoom.artefactborderpx)
        if error:
            return error,'error in uniformity'

        error = qc_unif.Artefacts(cs.unif,bpx)
        if error:
            return error,'error in artefacts'
        return error,''
