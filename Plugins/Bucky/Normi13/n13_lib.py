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
#
from __future__ import print_function
"""
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED! And make sure rescaling is corrected!

Note: comparison will be against lit.stTable, if not matched (eg. overwritten by config) assume wall

TODO:
Changelog:
    20180501: Detect infinite loop in CuWedge
    20180205: fix in n13_geometry to allow finding droplines at two heights; added extra param mustbeprecropped to room  
    20180124: fix in uniformity where border px was ignored if cropping detected
    20171116: fix scipy version 1.0
    20170825: added optional dicom header fields (should at some point replace the kludge of checking for modality)
    20170731: shrink xrayfield to exclude constant outside region; add param for mirroring of images; 
              if crop_frac>0.98, likely invert ratio (swapped fore and background);
              fix cropping of images where original seach area is empty (XA all 0.)
    20170629: decreased threshold for orientation for 90% to 66%
    20170623: fix for finding wrong grid line in 1 direction; bot more robust MTF
    20170619: fix for missing step in Cu Wedge due to noise; fix for artefacts if treshold close to 0; 
              fix for low crontrast out-of-image; increase box search range; fix xray edge found too early;
              fix for wrong orientation if outside phantom included
    20170518: _findDropLine now uses median which is more robust for noise
    20170324: made Geometry._FineTunePhantomBox bit more robust (prefer small shift); also do not quit if phantomGrid 
              found with too little confidence
    20170310: add override params for inversion and pixmm; geometry changed xray-edge finding logic; draw uniformity_crop on image; 
    20161220: Removed class variables; removed testing stuff
    20160816: split in separate files for each block of analysis
    20160812: another attempt at consistent box finding
    20160811: fix bad align for too low contrast; trendremoval; hessian; bugfix rotated phantom; fixes for small detector
    20160802: sync with wad2.0
    20160701: Fix gridscale for none DX; fixed TableWall none DX; laxer confidence; fix for invert; fix ints
    20160205: Distinguish between linepairs insert typ38 and RXT02
    20160202: added uniformity
    20151109: start of new module, based on QCXRay_lib of Bucky_PEHAMED_Wellhofer of 20151029
"""
__version__ = '20180501'
__author__ = 'aschilham'

try:
    import pydicom as dicom
except ImportError:
    import dicom
import numpy as np
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

try:
    # wad2.0 runs each module stand alone
    import n13_constants as lit
    import n13_math as mymath
    import n13_geometry as Geometry
    import n13_cuwedge as CuWedge
    import n13_lowcontrast as LowContrast
    import n13_resolution as Resolution
    import unif_lib
except ImportError:
    from . import n13_constants as lit
    from . import n13_math as mymath
    from . import n13_geometry as Geometry
    from . import n13_cuwedge as CuWedge
    from . import n13_lowcontrast as LowContrast
    from . import n13_resolution as Resolution
    from . import unif_lib
    
class Room:
    def __init__(self,_name, outvalue=-1, pid_tw=[-1,-1], sid_tw=[-1,-1],
                  linepairmarkers = {}, artefactborderpx=[0,0,0,0], detectorname={}, auto_suffix=False):
        self.name = _name        # identifier of room
        self.pidmm = {}
        self.sidmm = {}

        # xray field edge detection
        self.outvalue = outvalue # value of pixels outside x-ray field
        self.skip_cropping = False # For images with a circular FOV, this should be set to True

        # hard overrides for use_ params
        self.pixmm = None          # allow hard override of pixmm, for example is ImagerPixelSpacing does not exist
        self.mustbeinverted = None # allow hard override of auto invert
        self.mustbemirrored = False # by default do not mirror image; must be hard overridden if to do
        self.mustbeprecropped = None # allow start with a hard crop of the image; for example if the auto crop fails, box = [xmin_px,xmax_px, ymin_px,ymax_px]

        # 
        if len(pid_tw) == 1: # forced
            self.pidmm[lit.stForced] = pid_tw[0] 
            self.sidmm[lit.stForced] = sid_tw[0]
        else:
            self.pidmm[lit.stTable]    = pid_tw[0] # distance between mid phantom and detector in mm
            self.pidmm[lit.stWall]     = pid_tw[1]
            self.sidmm[lit.stTable]    = sid_tw[0] # distance between source and detector in mm
            self.sidmm[lit.stWall]     = sid_tw[1]

        # uniformity
        self.artefactborderpx = artefactborderpx
        self.artefactborder_is_circle = False # No, the artefactborder is not circular
        # MTF
        if len(linepairmarkers)>0:
            self.linepairmodel = linepairmarkers['type']
            if self.linepairmodel == 'RXT02':
                self.xy06mm = linepairmarkers['xymm0.6']
                self.xy10mm = linepairmarkers['xymm1.0'] # x,y position in mm of decimal dot in 1.0 lp/mm 
            elif self.linepairmodel == 'typ38':
                self.xy06mm = linepairmarkers['xymm0.6'] # x,y position in mm of decimal dot in 0.6 lp/mm 
                self.xy14mm = linepairmarkers['xymm1.4'] # x,y position in mm of decimal dot in 1.4 lp/mm 
                self.xy18mm = linepairmarkers['xymm1.8'] # x,y position in mm of decimal dot in 1.8 lp/mm 
                self.xy46mm = linepairmarkers['xymm4.6'] # x,y position in mm of decimal dot in 4.6 lp/mm 
            elif self.linepairmodel == 'None':
                pass
            else:
                raise ValueError('[Room] Unknown linepairmodel')
            
        # for auto_suffix
        self.detector_name = detectorname # a dict of [detectorid] = name like 
        self.auto_suffix = auto_suffix # DetectorSuffix will return None if not set to True

class XRayStruct:
    ###
    # class variables
    roomUnknown = Room(lit.stUnknown)

    def PreCropImage(self):
        """
        Apply a precropping, needed for example if auto cropping fails.
        expects param room.mustbeprecropped = [xmin_px, xmax_px, ymin_px, ymax_px]
        """
        error = False
        if self.forceRoom.mustbeprecropped is None:
            return error
        if self.pixeldataIn is None:
            return error

        [xmin_px,xmax_px, ymin_px,ymax_px] = self.forceRoom.mustbeprecropped
        widthpx, heightpx = np.shape(self.pixeldataIn)
    
        if xmin_px<0 or ymin_px<0:
            return True
    
        if xmax_px>=widthpx or ymax_px>=heightpx:
            return True
        
        self.pixeldataIn = self.pixeldataIn[xmin_px:xmax_px+1,ymin_px:ymax_px+1]
        return error

    def FixInvertedImage(self):
        """
        DICOM can store images with inverted scale. We demand high values for no attenuation (air).
        Check the dicomtags and invert the pixels in needed.
        """
        error = True

        if self.dcmInfile is None:
            return error

        if self._fixed_inversion: # already fixed
            return False
        
        if not self._fixed_mirror: # already fixed, next time skip this step
            if self.forceRoom.mustbemirrored:
                print("Must be Mirrored (use)",self.forceRoom.mustbemirrored)
                self.original_mirrored = self.forceRoom.mustbemirrored
                if not self.pixeldataIn is None:
                    self.pixeldataIn = np.fliplr(self.pixeldataIn)
    
            self._fixed_mirror = True # already fixed, next time skip this step
            
            
        if not self.forceRoom.mustbeinverted is None:
            print("Must be Inverted (use)",self.forceRoom.mustbeinverted)
            self.original_inverted = self.forceRoom.mustbeinverted
        else:
            # determine max allowed pixelvalue; as 2^bits_stored -1; note this value is not properly stored in dcmfileIn.BitsStored!
            dicomfields = [ ["0028,0101",  "Bits Stored"]]
            key = dicomfields[0][0]
            dicvalue = self.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
            if dicvalue == "":
                return error
            self.max_pixel_value = (2**dicvalue)-1
    
            # now check storage format
            self.original_inverted = False
            if self.dcmInfile.PhotometricInterpretation == "MONOCHROME2":
                self.original_inverted = True

        if self.original_inverted:
            # inversion
            if not self.pixeldataIn is None:
                self.pixeldataIn = self.max_pixel_value-self.pixeldataIn

        self._fixed_inversion = True # already fixed, next time skip this step

        error = False
        print("DICOM in inverted format: ",self.original_inverted)
        return error

    def pixToGridScale_mm(self):
        # determine scaling of pixel to mm on phantom scale
        if not self.forceRoom.pixmm is None:
            return self.forceRoom.pixmm
        
        pixel_spacing_x = self.dcmInfile.ImagerPixelSpacing[0] # PixelSpacing already corrects for magnification in DX!
        stand = self.DetectorStand()
        # DX
        try:
            sid = self.dcmInfile.DistanceSourceToDetector # source to image distance
        except:
            sid = self.forceRoom.sidmm[stand]
        try:
            sip = self.dcmInfile.DistanceSourceToPatient  # source to patient (=table top?!)
        except:
            pid = self.forceRoom.pidmm[stand]
            sip = sid-pid
        return pixel_spacing_x*sip/sid

    def DetectorSuffix(self):
        # return a suffix for results. This is deprecated in wad2.0. instead make multiple selectors, and combine later.
        # see if we can define by detectorID (DX only)
        if not self.forceRoom.auto_suffix:
            return None
        
        if len(self.forceRoom.detector_name)>0: # if this list is defined, then we must use it
            if hasattr(self.dcmInfile, 'DetectorID'):
                detectorid = self.dcmInfile.DetectorID # 0018,700A
            else:
                raise ValueError('Unknown detector: detectorID not defined')

            if detectorid in self.forceRoom.detector_name:
                return self.forceRoom.detector_name[detectorid]
            else:
                raise ValueError('[DetectorSuffix] Unknown detector %s'%detectorid)
        else:
            return self.knownDetectorStand
        
    def DetectorStand(self):
        # find out (based on SID) if Table or Wall stand is used, and which detector (could be wireless)
        if self.knownDetectorStand is not None:
            return self.knownDetectorStand

        # find out if user forces a pid
        if lit.stForced in self.forceRoom.pidmm:
            self.knownDetectorStand = lit.stForced
            return self.knownDetectorStand
          
        # some scheme for auto detection of detector size could be done by ....
        if hasattr(self.dcmInfile, 'DetectorManufacturerModelName'):
            #PIXIUM3543EZ = Large Portable
            #PIXIUM2430EZ = Small Portable
            #PIXIUM4343RC = Fixed detector
            pass
        
        # it is not a forced detector stand, try to distinguish wall of table by distance source to detector
        self.knownDetectorStand = lit.stUnknown

        try:
            sid = self.dcmInfile.DistanceSourceToDetector
            if sid>1600.:
                self.knownDetectorStand = lit.stWall
            else:
                self.knownDetectorStand = lit.stTable
        except:
            raise ValueError('Unknown detector. Cannot determine by distance')

        return self.knownDetectorStand


    def pix2phantommm(self, px):
        # convenience function
        return px*self.phantom_px_in_mm

    def phantommm2pix(self, px):
        # convenience function
        return px/self.phantom_px_in_mm

    def __init__ (self, dcmInfile, pixeldataIn, room):
        ###
        self.verbose = False
        self.knownDetectorStand = None

        # input image
        self.dcmInfile   = dcmInfile
        self.pixeldataIn = pixeldataIn

        self.original_inverted = False # if pixval(Cu) = high and pixval(air) = low, then mustbeinverted
        self._fixed_inversion = False # already fixed inversion, skip checks
        self._fixed_mirror = False # already fixed mirror, skip checks
        self.max_pixel_value = 0 # max allowed value to store per pixel

        self.forceRoom = room

        # first correct image represenation
        error = self.FixInvertedImage()
        if error:
            # not a valid dicom image
            print('ERROR. Not a valid DICOM image')
            return None

        # apply precropping if needed
        error = self.PreCropImage()
        if error:
            print('ERROR. PreCrop parameters not valid')
            return None
        
        self.phantom_px_in_mm = self.pixToGridScale_mm()

        # all geometry related stuff in geom
        self.geom = Geometry.GeomStruct()

        # Cu Wedge
        self.cuwedge = CuWedge.CuStruct()

        # Low Contrast
        self.loco = LowContrast.LoCoStruct()        

        # MTF
        self.mtf = Resolution.MTFStruct()

        # Uniformity # actually not part of Normi13
        self.unif = None

        # for matlib plotting
        self.hasmadeplots = False

        # GUI feedback
        self.lastimage = None 


class XRayQC:
    def __init__(self):
        self.qcversion = __version__
        pass

    def drawThickCircle(self,draw,x,y,rad,color,thick):
        for t in range(-int((thick-1)/2),int((thick+1)/2)):
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

        for t in range(-int((thick-1)/2),int((thick+1)/2)):
            draw.rectangle([(x0+t,y0+t),(x1-t,y1-t)],outline=color)

    def saveAnnotatedImage(self, cs, fname, what):
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        rectrois = []
        polyrois = []
        circlerois = []

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        if what == 'normi13':
            # first the base image
            im = scipy.misc.toimage(cs.pixeldataIn.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            # now add all box rois
            if len(cs.geom.box_roi) >0:
                polyrois.append(cs.geom.box_roi) # phantom orientation box
            if len(cs.geom.xr_roi) >0:
                polyrois.append(cs.geom.xr_roi) # xray edges
            if len(cs.cuwedge.box_roi) >0:
                polyrois.append(cs.cuwedge.box_roi) # Cu wedge box
            for r in cs.cuwedge.step_rois:
                polyrois.append(r) # Cu wedge element
            if len(cs.mtf.roi) >0:
                polyrois.append(cs.mtf.roi) # MTF box

            # add circlerois
            for r in cs.loco.low_rois: # low contrast elements
                circlerois.append(r)
            for r in cs.loco.low_rois_bku:
                circlerois.append(r)
            for r in cs.loco.low_rois_bkd:
                circlerois.append(r)

        else:
            unif_crop = cs.unif.unif_crop # xmin, xmax, ymin, ymax
            art_crop  = cs.unif.art_crop # xmin, xmax, ymin, ymax
            if what == 'uniformity':
                # first the base image
                # remove bk trend
                wid,hei = np.shape(cs.pixeldataIn)
                wid = int(wid/3)
                hei = int(hei/3)
                pdCopy = cs.pixeldataIn.copy()
                mean = np.mean(cs.pixeldataIn[wid:2*wid:3, hei:2*hei:3])
                std  = np.std(cs.pixeldataIn[wid:2*wid:3, hei:2*hei:3])
                pdCopy = (cs.pixeldataIn-mean)/std
                cut = 2
                mask = pdCopy<-cut
                pdCopy[mask] = -cut
                mask = pdCopy>cut
                pdCopy[mask] = cut

                #im = scipy.misc.toimage(cs.pixeldataIn.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!
                im = scipy.misc.toimage(pdCopy.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            elif what == 'artefacts':
                # first copy the artimage into the base image
                wid,hei = np.shape(cs.unif.art_image)
                pdCopy = np.zeros(np.shape(cs.pixeldataIn))
                pdCopy[art_crop[0]:art_crop[1],art_crop[2]:art_crop[3]] = cs.unif.art_image
                im = scipy.misc.toimage(pdCopy.transpose(),low=1,pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            # now add all box rois
            # add uniformity crop area
            rectrois.append([(unif_crop[0],unif_crop[2]),(unif_crop[1],unif_crop[3])]) #[(x0,y0),(x1,y1)]

            # show the box of the whole image so we can see what has been cropped
            wid,hei = np.shape(cs.pixeldataIn)
            rectrois.append([(0,0),(wid-1,hei-1)]) #[(x0,y0),(x1,y1)]

            # uniformity boxes
            for r in cs.unif.unif_rois: #[x0,dx,y0,dy]
                rectrois.append([(r[0],r[2]),(r[0]+r[1],r[2]+r[3])])#[(x0,y0),(x1,y1)]

            # artefact circles
            for r in cs.unif.art_rois: # x,y,r
                circlerois.append([r[0]+art_crop[0],r[1]+art_crop[2],r[2]])


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

    #------------------- new
    def MTF(self, cs):
        # Calulate MTF from line pairs element
        return Resolution.MTF(cs)

    def LowContrast(self, cs):
        # Calulate low contrast
        return LowContrast.LowContrast(cs)

    def CuWedge(self, cs):
        # Calulate stuff from the Cu wedge, like dynamical range
        return CuWedge.CuWedge(cs)

    def XRayField(self, cs):
        # Find edges of XRay exposure
        # note that this should be done on uncropped image!
        return Geometry.XRayField(cs)

    def XRayDev(self, cs):
        """
        Convenience function to calculate max deviation of xray edge with mean in %
        """
        stand = cs.DetectorStand()
        minedge = min(cs.geom.xr_NSWEmm)
        maxedge = max(cs.geom.xr_NSWEmm)
        meanedge = np.mean(cs.geom.xr_NSWEmm)
        if not 'DistanceSourceToDetector' in cs.dcmInfile:
            sidmm = cs.forceRoom.sidmm[stand]
        else:
            sidmm = cs.dcmInfile.DistanceSourceToDetector
            if sidmm<1.e-6:
                sidmm = cs.forceRoom.sidmm[stand]

        devedge = 100.*max(np.abs(minedge-meanedge),np.abs(maxedge-meanedge))/sidmm
        if maxedge-meanedge < meanedge-minedge:
            devedge *= -1
        return devedge

    def FindPhantomGrid(self, cs):
        # Resolve Phantom Coordinates
        return Geometry.FindPhantomGrid(cs)

    def FixPhantomOrientation(self, cs):
        # Rotate over 90, 180, 270 if needed
        return Geometry.FixPhantomOrientation(cs)

    def CropNormi13(self, cs):
        # Crop image to phantom if needed
        return Geometry.CropPhantom(cs)

    def DICOMInfo(self, cs, info='dicom'):
        # Different from ImageJ version; tags "0008","0104" and "0054","0220"
        #  appear to be part of sequences. This gives problems (cannot be found
        #  or returning whole sequence blocks)
        # Possibly this can be solved by using if(type(value) == type(dicom.sequence.Sequence()))
        #  but I don't see the relevance of these tags anymore, so set them to NO

        if info == "dicom":
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
                {'key':"0008,1090",  'name':"Modelname"},
                {'key':"0010,0020",  'name':"PatientID"},
                {'key':"0018,0015",  'name':"BodyPartExamined"},
                {'key':"0018,0060",  'name':"kVp"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)"},
                {'key':"0018,115E",  'name':"ImageAreaDoseProduct"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"0018,1166",  'name':"Grid"},
                {'key':"0018,6000",  'name':"Sensitivity"},
            ]
            if cs.dcmInfile.Modality == 'CR':
                if not 'DistanceSourceToDetector' in cs.dcmInfile: # WKZ-like fcr
                    dicomfields.extend( [
                        {'key':"0018,1004",  'name':"PlateID"},
                        {'key':"0018,1401",  'name':"AcquisitionDeviceProcessingCode"},
                        {'key':"0018,1403",  'name':"CassetteSize"},
                        {'key':"0018,1404",  'name':"ExposuresOnPlate"},
                        {'key':"0018,1508",  'name':"PositionerType"},
                        {'key':"0028,0006",  'name':"PlanarConfiguration"},
                    ])

                dicomfields.extend( [
                    {'key':"0018,1150",  'name':"ExposureTime (ms)"},
                    {'key':"0018,1152",  'name':"Exposure (mAs)"},
                    {'key':"0018,1160",  'name':"FilterType"},
                    {'key':"0018,1190",  'name':"FocalSpot(s)"},
                    {'key':"0018,1260",  'name':"PlateType"},
                    {'key':"0018,1200",  'name':"Date of Last Calibration"},
                    {'key':"0018,1400",  'name':"AcquisitionDeviceProcessingDescription"},

                ])
            elif cs.dcmInfile.Modality == 'DX':
                dicomfields.extend( [
                    {'key':"0018,1030",  'name':"ProtocolName"},
                    {'key':"0018,1111",  'name':"DistanceSourceToPatient (mm)"},
                    {'key':"0018,1153",  'name':"Exposure (uAs)"},

                    {'key':"0018,1405",  'name':"RelativeXRayExposure"},
                    {'key':"0018,1702",  'name':"CollimatorLeft"},
                    {'key':"0018,1704",  'name':"CollimatorRight"},
                    {'key':"0018,1706",  'name':"CollimatorUp"},
                    {'key':"0018,1708",  'name':"CollimatorDown"},

                    {'key':"0018,7001",  'name':"DetectorTemperature"},
                    {'key':"0018,700A",  'name':"DetectorID"},
                    {'key':"0018,700C",  'name':"DateCalibration"},
                    {'key':"0018,702B",  'name':"ModelName"},
                    {'key':"0018,7050",  'name':"FilterMaterial"},
                    {'key':"0018,7062",  'name':"ExposureControlMode"},
                    {'key':"0018,8150",  'name':"ExposureTime (us)"},
                    {'key':"0028,0101",  'name':"BitsStored"},
                    {'key':"0040,8302",  'name':"EntranceDose_mGy"},
                    {'key':"200B,70BA",  'name':"FocalSpot"},
                    {'key':"200B,1028",  'name':"PrivValue1"},
                    {'key':"200B,7096",  'name':"PrivStand"},
                ])
                

        elif info == "qcwad":
            offset = -25 # rank must be negative, so recalc as offset+real position
            dicomfields = [
                {'key':"0008,0021",  'name':"SeriesDate"},
                {'key':"0008,0031",  'name':"SeriesTime", 'quantity':'time', 'level':1, 'rank':offset+2}, # spot 1 reserved for stand
                {'key':"0018,1110",  'name':"DistanceSourceToDetector (mm)", 'quantity':'distance', 'level':1, 'rank':offset+3},
                {'key':"0018,0060",  'name':"kVp", 'level':1, 'rank':offset+4},
                {'key':"0018,1166",  'name':"Grid", 'quantity':'grid', 'level':1, 'rank':offset+6},

                {'key':"0018,115E",  'name':"ImageAreaDoseProduct", 'quantity':'DAP','level':1,'rank':offset+14},
                {'key':"0018,6000",  'name':"Sensitivity",'quantity':'S'},

                {'key':"0008,103E",  'name':"SeriesDescription"},
                {'key':"0008,1010",  'name':"StationName"},
                {'key':"0008,1070",  'name':"Operator's Name"},
                {'key':"0018,1000",  'name':"DeviceSerialNumber"},
                {'key':"0018,1020",  'name':"SoftwareVersions"},
                {'key':"0018,1164",  'name':"ImagerPixelSpacing"},
                {'key':"0018,0015",  'name':"BodyPartExamined"},
                {'key':"0018,1030",  'name':"ProtocolName"},
            ]
            if cs.dcmInfile.Modality == 'CR':
                dicomfields.extend([
                ])
                if not 'DistanceSourceToDetector' in cs.dcmInfile: # WKZ-like fcr
                    dicomfields.extend( [
                        {'key':"0018,1004",  'name':"Plate ID"},
                        {'key':"0018,1401",  'name':"Acquisition Device Processing Code", 'quantity':'processing', 'level':1, 'rank':offset+9},
                    ])
                else:
                    dicomfields.extend( [
                        {'key':"0018,1160",  'name':"FilterType", 'quantity':'filter', 'level':1, 'rank':offset+5},
                        {'key':"0018,1190",  'name':"FocalSpot(s)", 'quantity':'focalspot', 'level':1, 'rank':offset+7},
                        {'key':"0018,5021",  'name':"Postprocessing", 'quantity':'processing', 'level':1, 'rank':offset+9}, # spot 8 reserved for rotation
                        {'key':"0018,1150",  'name':"ExposureTime (ms)", 'quantity':'ms','level':1,'rank':offset+12},
                        {'key':"0018,1152",  'name':"Exposure (mAs)", 'quantity':'mAs','level':1,'rank':offset+13},
                        {'key':"0018,1200",  'name':"Date of Last Calibration"},
                    ])

            elif cs.dcmInfile.Modality == 'DX':
                dicomfields.extend([
                    {'key':"0018,7050",  'name':"FilterMaterial", 'quantity':'filter', 'level':1, 'rank':offset+5},
                    {'key':"200B,70BA",  'name':"FocalSpot", 'quantity':'spot', 'level':1, 'rank':offset+7},
                    {'key':"0018,1702",  'name':"CollimatorLeft",'quantity':'Left','level':1, 'rank':offset+8},
                    {'key':"0018,1704",  'name':"CollimatorRight",'quantity':'Right','level':1, 'rank':offset+9},
                    {'key':"0018,1706",  'name':"CollimatorUp",'quantity':'Up','level':1, 'rank':offset+10},
                    {'key':"0018,1708",  'name':"CollimatorDown",'quantity':'Down','level':1, 'rank':offset+11},
                    {'key':"0018,8150",  'name':"ExposureTime (us)", 'quantity':'us','level':1,'rank':offset+12},
                    {'key':"0018,1153",  'name':"Exposure (uAs)", 'quantity':'uAs','level':1,'rank':offset+13},
                    {'key':"0018,700A",  'name':"DetectorID"},
                    {'key':"0018,700C",  'name':"DateCalibration"},
                    {'key':"200B,7063",  'name':"Postprocessing", 'quantity':'processing'},
                    {'key':"0018,1405",  'name':"RelativeXRayExposure"},
                    {'key':"0040,8302",  'name':"EntranceDose_mGy"},
                    {'key':"200B,1028",  'name':"PrivValue1"},
                    {'key':"200B,7096",  'name':"PrivStand"},
                ])

        # optional fields: if not already in dicomfields, and if exist in header, then include
        opt_fields = [
            {'key':"0018,1153",  'name':"Exposure (uAs)"},
            {'key':"0018,1405",  'name':"RelativeXRayExposure"},
            {'key':"0018,8150",  'name':"ExposureTime (us)"},
            {'key':"0018,1150",  'name':"ExposureTime"},
            {'key':"200B,70BA",  'name':"FocalSpot"},
            {'key':"0018,1190",  'name':"Focal Spot(s)"},
            {'key':"0028,0101",  'name':"BitsStored"},
            {'key':"0040,8302",  'name':"EntranceDose_mGy"},
            {'key':"0018,7050",  'name':"FilterMaterial"},
            {'key':"0018,1160",  'name':"FilterType"},
            {'key':"0018,702B",  'name':"ModelName"},
            {'key':"0018,7001",  'name':"DetectorTemperature"},
            {'key':"0018,700A",  'name':"DetectorID"},
            {'key':"0018,700C",  'name':"DateCalibration"},
            {'key':"0018,1200",  'name':"Date of Last Calibration"},

            {'key':"0018,1405",  'name':"RelativeXRayExposure"},
            {'key':"0018,1411",  'name':"ExposureIndex"},
            {'key':"0018,1412",  'name':"TargetExposureIndex"},
            {'key':"0018,1413",  'name':"DeviationIndex"},
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

        df_keys = [ k['key'].upper() for k in dicomfields ]
        for df in opt_fields:
            key = df['key'].upper()
            
            # skip if already in list
            if key in df_keys:
                continue
        
            # only include if present
            try:
                value =  cs.dcmInfile[dicom.tag.Tag(key.split(',')[0],key.split(',')[1])].value
                df['value'] = value
                results.append( df )
            except:
                continue
        
            
        return results

    def ReportEntries(self,cs):
        """
        Convenience function to list all calculated items with names, only if cs completely filled
        """
        labvals = []
        if not cs.unif is None:
            num_art = 0
            if cs.unif.art_clusters != None:
                num_art = len(cs.unif.art_clusters)

            offset = -5 # rank must be negative, so recalc as offset+real position

            labvals.append( {'name':'Uniformity_(%)','value':cs.unif.unif_pct, 'quantity':'Uniformity','level':1,'rank':offset} )
            labvals.append( {'name':'Artefacts','value':num_art, 'quantity':'Artefacts','level':1,'rank':offset+1} )

            for kk in range(0,len(cs.unif.means)):
                labvals.append( {'name':'avg_'+str(kk),'value':cs.unif.means[kk], 'quantity':'Mean','level':2} )
            for kk in range(0,len(cs.unif.stdevs)):
                labvals.append( {'name':'sd_'+str(kk),'value':cs.unif.stdevs[kk], 'quantity':'STDev','level':2} )
            for lab,val in zip(['xmin','xmax','ymin','ymax'],cs.unif.art_borderpx):
                labvals.append( {'name':'crop_%s'%lab,'value':val, 'quantity':'borderpx','level':2} )
            labvals.append( {'name':'Cropped_(%)','value':100*cs.unif.unif_crop_frac, 'quantity':'Cropped','level':2} )

            return labvals

        ## Phantom orientation; level 1 = show by default; level 2 = show in details
        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','rank':missing or a negative number} )
        # if no rank given, the order of addition will be used
        # if no quantity given, 'name' will be used
        # if no level given, the default will be used
        offset = -30 # rank must be negative, so recalc as offset+real position
        labvals.append( {'name':'PhantomOrientation','value':cs.geom.box_orientation, 'quantity':'angle','level':1,'rank':offset+7} )
        labvals.append( {'name':'AlignConfidence','value':100.*cs.geom.box_confidence, 'quantity':'aligned','level':2} )
        labvals.append( {'name':'xray[N]cm','value':cs.geom.xr_NSWEmm[0]/10., 'quantity':'xrN','level':2} )
        labvals.append( {'name':'xray[E]cm','value':cs.geom.xr_NSWEmm[3]/10., 'quantity':'xrE','level':2} )
        labvals.append( {'name':'xray[S]cm','value':cs.geom.xr_NSWEmm[1]/10., 'quantity':'xrS','level':2} )
        labvals.append( {'name':'xray[W]cm','value':cs.geom.xr_NSWEmm[2]/10., 'quantity':'xrW','level':2} )
        labvals.append( {'name':'xrayDev%','value':self.XRayDev(cs), 'quantity':'xrayDev','level':1, 'rank':offset+12} )
        labvals.append( {'name':'centershift_x','value':cs.geom.center_shiftxy[0], 'quantity':'shiftxpx','level':1} )
        labvals.append( {'name':'centershift_y','value':cs.geom.center_shiftxy[1], 'quantity':'shiftypx','level':1} )

        ## cuwedge
        labvals.append( {'name':'CuConfidence','value':cs.cuwedge.wedge_confidence*100,'level':2} ) # Confidence in wedge finding
        # SNR max
        labvals.append( {'name':'CuSNR_'+str(cs.cuwedge.step_mmcu[-1]),'value':cs.cuwedge.step_snr[-1], 'quantity':'SNR','level':1,'rank':offset+15} )
        # CNR between steps all > 1
        minCNR = cs.cuwedge.step_cnr[0]
        for i in range(1,len(cs.cuwedge.step_cnr)-1):
            minCNR = min (minCNR,cs.cuwedge.step_cnr[i])
        labvals.append( {'name':'CuCNRmin','value':minCNR, 'quantity':'CNRmin','level':2} )
        # Dynamic Range
        labvals.append( {'name':'CuDR'+str(cs.cuwedge.step_mmcu[0])+'_'+str(cs.cuwedge.step_mmcu[-1]),'value':cs.cuwedge.dynamicRange, 'quantity':'DynRange','level':1,'rank':offset+14} )

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

    def QC(self, cs):
        """
        Normi13 analysis all in one routine.

        Outline:
        (0. invert image or not; already in creation of qc struct!)
        1. geometry:
          1.1 Crop
          1.2 Phantom90Degrees
          1.3 PhantomCoordinates
          1.4 XRayEdges (check full field)
        2. Cu Wedge analysis
        3. Low Contrast analysis
        4. MTF analysis
        """
        error = True
        msg = ''

        #print('[QCNormi13]',cs.dcmInfile.SeriesDescription)

        # 1.1 geometry: crop
        error = self.CropNormi13(cs)
        if error:
            msg += 'Crop '
        
        # 1.2 geometry: orientation
        if not error:
            error = self.FixPhantomOrientation(cs)
            if error:
                msg += 'Orientation '

        # 1.3 geometry: phantom coordinates
        if not error:
            error = self.FindPhantomGrid(cs)
            if error:
                msg += 'Grid '
                error = False

        # 1.4: travel straight along NS and find edge of x-ray; similar for EW
        if not error:
            error = self.XRayField(cs)
            if error:
                msg += 'XRayField '

        # 2: find Cu wedge stuff
        if not error:
            error = self.CuWedge(cs)
            if error:
                msg += 'CuWedge '

        # 3: low contrast stuff
        if not error:
            error = self.LowContrast(cs)
            if error:
                msg += 'LowContrast '

        # 4: find resolution stuff
        if not error:
            error = self.MTF(cs)
            if error:
                msg += 'MTF '

        return error,msg

    # -------------------
    def Uniformity(self,cs):
        # run external QCUniformity
        usestructure = True

        cs.unif = unif_lib.UnifStruct(cs.dcmInfile, cs.pixeldataIn)
        cs.unif.pixmm = cs.forceRoom.pixmm
        
        qc_unif = unif_lib.Uniformity_QC()

        if not cs.forceRoom.skip_cropping and qc_unif.NeedsCroppingUnif(cs):
            # note: no cropping will occur, just the uniformity analysis will be restricted to crop_unif. rois are wrt original
            qc_unif.RestrictROIUniformity(cs.unif)
        else:
            borderpx=[0,0,0,0]
            widthpx  = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
            heightpx = np.shape(cs.pixeldataIn)[1]
            xmin0_px = borderpx[0]
            xmax0_px = widthpx-borderpx[1]-1
            ymin0_px = borderpx[2]# cs.border_offset_px
            ymax0_px = heightpx-borderpx[3] -1 #cs.border_offset_px -1
            cs.unif.unif_crop = [xmin0_px,xmax0_px, ymin0_px,ymax0_px]
            cs.unif.unif_crop_frac = 1.
            cs.unif_crop_inoutoverin = 1.


        if usestructure:
            qc_unif.artefactDetectorParameters(UseStructure=True, bkscale=25, fgscale=5.0, threshold=3000)
        else:
            qc_unif.artefactDetectorParameters(UseStructure=False, bkscale=25, fgscale=5.0, threshold=15)

        error = qc_unif.Uniformity(cs.unif, cs.forceRoom.artefactborderpx, cs.forceRoom.artefactborder_is_circle) # DOES NOT HAVE TO BE THE SAME BORDERPIX

        if error:
            return error,'error in uniformity'

        # here a little bit of cropping, but only for the cs.artefact_image. artefacts are wrt cropped cs.artefact_image
        error = qc_unif.Artefacts(cs.unif, borderpx=cs.forceRoom.artefactborderpx, 
                                  border_is_circle=cs.forceRoom.artefactborder_is_circle)
        if error:
            return error,'error in artefacts'

        return error,''

    def QCUnif(self, cs):
        """
        Separate uniformity in one step!
        """
        error = True
        msg = ''

        label = cs.dcmInfile.get('SeriesDescription', None)
        if label is None:
            label = cs.dcmInfile.get('BodyPartExamined', None)
            
        print('[QCUnif]', label)
        error, msg = self.Uniformity(cs)

        if error:
            msg += "Uniformity (try)) "
            print(msg)
            return error, msg

        return error, msg

