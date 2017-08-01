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
Normi13 analysis, geometry functions:
  o Cropping (restrict analysis field to part of image covered with Normi13 phantom)
  o 90 degree rotation of phantom (check and correct if phantom is rotated by 0, 90, 180, 270 degrees)
  o PhantomGrid (locate the phantom grid so that elements can be found by phantom coordinates)
  o XRayEdges (determine edges of XRayField in phantom coordinates (only meaningfull if not full field exposure))
"""
import copy
import numpy as np
import scipy.ndimage as scind
from scipy.signal import medfilt
import matplotlib.pyplot as plt

try:
    # wad2.0 runs each module stand alone
    import unif_lib
except ImportError:
    from . import unif_lib
    
class GeomStruct:
    def _roi_to_coords_transformation(self):
        """
        Use found phantombox to determine transformation from pos_xy in mm on phantom grid to pos_xy in px
        ymm is defined positive for top of phantom, and negative for bottom.
        box_roi = [ UL, LL, LR, UR]
        box_radmm = [half width, half height]
        """
        xmidpx = float(np.sum([x for x,y in self.box_roi]))/len(self.box_roi)
        ymidpx = float(np.sum([y for x,y in self.box_roi]))/len(self.box_roi)

        # calc difference in x, y in px per mm shift in x, y
        xmm_dxy = [0, 0]# 1 mm in x direction results in shift xy
        ymm_dxy = [0, 0]# 1 mm in y direction results in shift xy
        for xy in [0, 1]:
            dmm = 2*self.box_radmm[0] # distance travelled in mm along this dimension of box 
            dxpx = ( (self.box_roi[2][xy]-self.box_roi[1][xy])+(self.box_roi[3][xy]-self.box_roi[0][xy]) )/2.
            xmm_dxy[xy] = dxpx/dmm
            dmm = 2*self.box_radmm[1] # distance travelled in mm along this dimension of box 
            dypx = ( (self.box_roi[0][xy]-self.box_roi[1][xy])+(self.box_roi[3][xy]-self.box_roi[2][xy]) )/2.
            ymm_dxy[xy] = dypx/dmm

        return [xmidpx, ymidpx], xmm_dxy, ymm_dxy 

    def phantomposmm2pix(self, xmm, ymm):
        """
        Convert position as indicated on phantom to image pixels

        input: x and y coordinates in phantom grid mm
        output: x and y in pix
        """

        [xmid,ymid], xmm_d, ymm_d = self._roi_to_coords_transformation()

        xpx = xmid + xmm*xmm_d[0] + ymm*ymm_d[0]
        ypx = ymid + xmm*xmm_d[1] + ymm*ymm_d[1]

        return xpx ,ypx


    def __init__(self):
        self.crop_ranges = None    # [xmin_px,ymin_px, xmax_px,ymax_px]
        self.crop_frac = 1.        # area_in/area_all
        self.crop_inoutoverin = 0. #(mean_in-mean_out)/mean_in

        # phantombox  
        self.box_roi = []            # [UL, LL, LR, UR]
        self.box_radmm = []          # halfwidth, halfheight of phantombox in mm
        self.box_confidence = 0.
        self.box_orientation = 0     # rotation of box in degrees (0, 90, 180, 270)
        self.center_shiftxy = [0, 0] # shift of box center wrt to center of uncropped image in px
        self.xr_NSWEmm = []          # edges of xray exposure in mm
        self.xr_roi = []             # edges of xray exposure in px on cropped im
        self.orig_shape = []         # width,height in px of original image, so before cropping
        

def CropPhantom(cs):
    """
    Crop Image to Phantom if too much space around phantom
    """
    error = False

    cs_unif = unif_lib.UnifStruct(cs.dcmInfile, cs.pixeldataIn)
    cs_unif.pixmm = cs.forceRoom.pixmm
    qc_unif = unif_lib.Uniformity_QC()
    widthpx, heightpx = np.shape(cs.pixeldataIn)
    cs.geom.orig_shape = [widthpx, heightpx]

    if qc_unif.NeedsCropping2(cs_unif, mode='normi13'):
        print('*** Needs Cropping ***')        
        bpx = 0#self.phantommm2pix(cs,10)
        error = qc_unif.RestrictROINormi13(cs_unif)
        if error:
            return error, 'Could not crop to Normi13'

        [xmin_px,xmax_px, ymin_px,ymax_px] = cs_unif.unif_crop
        xmin_px = max(0,xmin_px-bpx)
        ymin_px = max(0,ymin_px-bpx)
        xmax_px = min(widthpx-1,xmax_px+bpx)
        ymax_px = min(heightpx-1,ymax_px+bpx)
        cs.pixeldataIn = cs.pixeldataIn[xmin_px:xmax_px+1,ymin_px:ymax_px+1]
        cs.geom.crop_ranges = [xmin_px,ymin_px, xmax_px,ymax_px]
        cs.geom.crop_frac = cs_unif.unif_crop_frac
        cs.geom.crop_inoutoverin = cs_unif.unif_crop_inoutoverin 
    return error

def FixPhantomOrientation(cs):
    """
    Concept: Try to find out if phantom is rotated over 90, 180 or 270 degrees (WKZ mostly)
    Workflow:
    1. Find center and orientation
    2. Find most likely location of CuWedge (thick)
    3. Use location wrt center to find and undo 90 degree rotations
    """
    error = False

    # 2. Make box at x-5cm and calculate avg
    seppx = int(cs.phantommm2pix(20))# 2 cm away from edges
    seppx = int(cs.phantommm2pix(10))# 2 cm away from edges

    widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
    heightpx = np.shape(cs.pixeldataIn)[1]

    still_in_edge = True
    while still_in_edge:
        seppx += int(cs.phantommm2pix(10))# 2 cm away from edges
        smallimage = cs.pixeldataIn[seppx:widthpx-seppx:3,seppx:heightpx-seppx:3] #every third pixel
        smallimage = scind.gaussian_filter(smallimage, sigma=5, )
    
        # find quandrant with min transmission; count which quad has most low transmission pix
        minim = smallimage.min()
        maxim = smallimage.max()
        imhei = smallimage.shape[1]
        imwid = smallimage.shape[0]
        imhei2 = int(imhei/2)
        imwid2 = int(imwid/2)
        thresmax = minim+.05*(maxim-minim) # 5% higher than minimum
        # make sure the outside of the phantom is not included
        if np.sum(smallimage[  0:imwid2,  0        ]<thresmax) > .66*imwid2:
            continue
        if np.sum(smallimage[  0:imwid2, -1        ]<thresmax) > .66*imwid2:
            continue
        if np.sum(smallimage[  0,        0:imhei2 ]<thresmax) > .66*imhei2:
            continue
        if np.sum(smallimage[ -1,        0:imhei2 ]<thresmax) > .66*imhei2:
            continue
        still_in_edge = False

        count = np.array([0,0,0,0])
        count[0] = np.sum(smallimage[     0:imwid2,      0:imhei2]<thresmax)
        count[1] = np.sum(smallimage[imwid2:imwid,       0:imhei2]<thresmax)
        count[2] = np.sum(smallimage[imwid2:imwid,  imhei2:imhei]<thresmax)
        count[3] = np.sum(smallimage[     0:imwid2, imhei2:imhei]<thresmax)

        count0 = count.argmax()

    if count0==0:
        # correct side on x
        # correct side on y
        ang = 0
    elif count0==1:
        # correct side on y
        ang = 1 # +90rot
    elif count0==2:
        ang = 2 # +180rot
    else:
        # correct side on x
        ang = 3 # +270rot
    
    # fix orientation
    box_orientation = ang
    if ang>0:
        cs.pixeldataIn = np.rot90(cs.pixeldataIn,-ang)
        box_orientation = 90*ang

    """
    error = self.findPhantomOrientation(cs)
    if error:
        # for GUI to show error image
        cs.lastimage = smallimage
    """

    # copy to geom
    cs.geom.box_orientation = box_orientation

    print('Phantom orientation: ',cs.geom.box_orientation)
    return error

def FindPhantomGrid(cs):
    """
    In order to locate objects on the phantom by the coordinates as shown on the phantom, a transformation 
    between phantom coordinates and pixel position should be found.
    Here, this is accomplished by finding the grid on the Normi13 phantom, working from the center of
    the phantom outwards.
    The shift of the center of the found box is stored in cs.geom.center_shiftxy

    The input image should be properly cropped (CropPhantom) and oriented (FixPhantomOrientation).
    """
    error = False
    # first try a to find 90x90 box
    roipts = _FindPhantomBox(cs, vertical=None, assumegood=False)

    # if it fails, try 70x90 box
    if roipts is None:
        print("Could not find 90x90 box, will try 90x70")
        roipts = _FindPhantomBox(cs, vertical=50, assumegood=False)

    # if it fails, try assume perfectly position and try to continue
    if roipts is None:
        print("Could not find 90x70 box, will assume 90x70")
        roipts = _FindPhantomBox(cs, vertical=50, assumegood=True)

    cs.geom.box_roi = roipts
    error, roipts, box_confidence = _FineTunePhantomBox(cs, roipts)

    cs.geom.box_roi = roipts
    cs.geom.box_confidence = box_confidence

    # 1.1 find approx center of phantom
    widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
    heightpx = np.shape(cs.pixeldataIn)[1]
    if cs.geom.crop_ranges is None:
        xmin_px = ymin_px = 0
    else:
        [xmin_px,ymin_px, xmax_px,ymax_px] = cs.geom.crop_ranges
    orig_midx = .5*(cs.geom.orig_shape[0]-1)
    orig_midy = .5*(cs.geom.orig_shape[1]-1)
    roi_midx = np.mean([x for x,y in roipts])
    roi_midy = np.mean([y for x,y in roipts])
    cs.geom.center_shiftxy = [roi_midx+xmin_px-orig_midx, roi_midy+ymin_px-orig_midy]
    return error

def _FindPhantomBox(cs, vertical=None, assumegood=False):
    """
    Draw lines from center outwards. Find first narrow valleys of 90 mm copper lines.
    North: must pass through Cu gap between steps; can be missing on small detector
    East: should always be findable
    South: can be missing on small detector
    West: will not be able to find because of line pairs element

    if vertical is given:
    Look at those y coordinates, and travel only form center to left and to right;
    this should be used in case the line finding fails for small detectors.
    """
    error = True

    # 1.1 find approx center of phantom
    widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
    heightpx = np.shape(cs.pixeldataIn)[1]
    midx = int(.5*(widthpx-1)+.5)
    midy = int(.5*(heightpx-1)+.5)

    # for testing a rotated input image can be forced
    # cs.pixeldataIn = scind.interpolation.rotate(cs.pixeldataIn, 2, axes=(1, 0), reshape=False, output=None, order=3, mode='constant', cval=0.0, prefilter=True)

    wid2_mm = 90 # half width of box in mm
    hei2_mm = 90 # half height of box in mm

    wid2_px = int(cs.phantommm2pix(wid2_mm) +.5) # half width of box in px
    hei2_px = int(cs.phantommm2pix(hei2_mm) +.5) # half height of box in px

    if assumegood: # if all fails, try to continue with a small detector anyway
        if not vertical is None:
            hei2_mm = vertical # half height of box in mm
            hei2_px = int(cs.phantommm2pix(hei2_mm) +.5) # half height of box in px

        midypx = int(midy)
        midxpx = int(midx)
    else:
        hlinepx = int(max(2,cs.phantommm2pix(1)+.5)) # don't expect a gridline to be larger than 2mm
        #seppx = int(cs.phantommm2pix(20))# stop 2 cm away from image edges
        seppx = int(cs.phantommm2pix(10))# stop 1 cm away from image edges (arbitrarily)
        blockheight = 10 # px number of lines to average (get id of noise!)

        if vertical is None:
            v_px = int(cs.phantommm2pix(70) +.5) # look for east/west at height of minus 70
            rois = []
            rois.append([midx,midx+blockheight, midy,seppx]) # North
            rois.append([midx,widthpx-1-seppx,  midy+v_px,midy+v_px+blockheight]) # East
            rois.append([midx,midx+blockheight, midy,heightpx-1-seppx]) # South
            rois.append([midx,seppx,  midy+v_px,midy+v_px+blockheight]) # West
            
            for factor in [2., 2.5, 3.]:
                lines = []
                edgepos = []
                for r in rois:
                    stepx = 1 if r[1]>r[0] else -1
                    stepy = 1 if r[3]>r[2] else -1
                    smallimage = cs.pixeldataIn[r[0]:r[1]:stepx,r[2]:r[3]:stepy]
    
                    ep, line, threshold = _findDropLine(smallimage, hlinepx, removeTrend=True, factor=factor)
                    if ep == -1:
                        print(threshold,line)
    
                    edgepos.append(ep)
                    lines.append(line)
                ns = np.sum( [cs.pix2phantommm(edgepos[i]) for i in [0,2]] )
                ew = np.sum( [cs.pix2phantommm(edgepos[i]) for i in [1,3]] )
                if np.abs(ew-ns)<5.: # if difference between ew and ns less than 0.5 cm, accept. else retry
                    break
                else:
                    print("Factor {}: Difference between NS and EW too large ({} mm). Will try again.".format(factor, ew-ns))

            if cs.verbose:
                cs.hasmadeplots = True
                plt.figure()
                for line,lab in zip(lines,['N','E','S','W']):
                    plt.plot(line,label=lab)

                plt.plot(edgepos,len(edgepos)*[threshold],'ro')
                print(edgepos)
                plt.legend()

            for p,k in zip(edgepos,['N','E','S','W']):
                print(k,p,cs.pix2phantommm(p))

            mp = cs.pix2phantommm(min(edgepos))
            if mp < 20: # edge was found within 2 cm of center... impossible
                return None

            # the twice 90 mm distance is the northedge+southedge
            hei2_px = int( (edgepos[0]+edgepos[2])/2. +.5 )
            wid2_px = int( (edgepos[1]+edgepos[3])/2. +.5 )
            #adjust phantom center
            midypx = int(midy -(edgepos[0]-edgepos[2])/2.+.5)
            midxpx = int(midx +(edgepos[1]-edgepos[3])/2.+.5)
        else:
            # travel only in horizontal direction, at plus and minus vertical grid mm from center
            # this does assume the phantom is properly centred vertically
            hei2_mm = vertical # half height of box in mm
            hei2_px = int(cs.phantommm2pix(hei2_mm) +.5) # half height of box in px

            rois = []
            rois.append([midx,seppx,            midy-hei2_px,midy-hei2_px+blockheight]) # top left
            rois.append([midx,widthpx-1-seppx,  midy-hei2_px,midy-hei2_px+blockheight]) # top right
            rois.append([midx,seppx,            midy+hei2_px,midy+hei2_px-blockheight]) # bottom left
            rois.append([midx,widthpx-1-seppx,  midy+hei2_px,midy+hei2_px-blockheight]) # bottom right
            while len(rois)<4:
                rois.append(rois[-1])
            lines = []
            edgepos = []
            for r in rois:
                stepx = 1 if r[1]>r[0] else -1
                stepy = 1 if r[3]>r[2] else -1
                smallimage = cs.pixeldataIn[r[0]:r[1]:stepx,r[2]:r[3]:stepy]
                ep, line, threshold = _findDropLine(smallimage, hlinepx, removeTrend=True)
                if ep == -1:
                    print(threshold,line)

                edgepos.append(ep)
                lines.append(line)
            if cs.verbose:
                cs.hasmadeplots = True
                plt.figure()
                for line,lab in zip(lines,['TL','TR','BL', 'BR']):
                    plt.plot(line,label=lab)
                plt.plot(edgepos,[threshold,threshold,threshold, threshold],'ro')
                print(edgepos)
                plt.legend()

            for p,k in zip(edgepos,['TL','TR','BL', 'BR']):
                print(k,p,cs.pix2phantommm(p))

            mp = cs.pix2phantommm(min(edgepos))
            if mp < 20: # edge was found within 2 cm of center... impossible
                print("Could not find 90 box, trying assumption of good position.")
                return None

            # the twice 90 mm distance is the northedge+southedge
            wid2_px = int( (edgepos[0]+edgepos[1]+edgepos[2]+edgepos[3])/4 +.5 )

            #adjust phantom center
            midypx = int(midy)
            midxpx = int(midx+ (edgepos[1]-edgepos[0]+edgepos[3]-edgepos[2])/4.)

    cs.geom.box_radmm = [wid2_mm, hei2_mm]
    # define cornerpoints of 90x90cm square
    roipts = [ 
        [midxpx-wid2_px,midypx-hei2_px], # ll
        [midxpx-wid2_px,midypx+hei2_px], # ul
        [midxpx+wid2_px,midypx+hei2_px], # ur
        [midxpx+wid2_px,midypx-hei2_px], # lr
    ]

    return roipts

def removeBKTrend(imageIn, sigma):
    # remove trend in background by subtraction wide median
    work = imageIn.astype(float) # must work in float images else gaussians are ints as well, which is bad for poor contrast
    sigma = max(3, 2*int(sigma/2)+1)
    work = 1000.*(work -medfilt(work, sigma)) 
    return work

def _findDropLine(smallimage, hlinepx, removeTrend=False, factor=2.):
    # find first position of drop for grid line
    # optionally remove a trend in the background

    axis = 0 if np.shape(smallimage)[0] < np.shape(smallimage)[1] else 1
    # turn each roi into a average line
    # line = np.average(smallimage,axis=axis)
    line = np.median(smallimage,axis=axis) # more robust for 1 px outlier
    
    # find first position of drop
    threshold = (np.amax(line)+np.amin(line))/2
    threshold = (line[0]+np.amin(line))/factor
    
    if removeTrend:
        line = removeBKTrend(line, 3*hlinepx)
        threshold = (0+np.amin(line))/factor # trend removal, so expect negative! sometimes 3!

    ep = -1
    for x in range(hlinepx,len(line)-hlinepx-1):
        if line[x-hlinepx]>threshold and line[x]<threshold and line[x+hlinepx]> threshold:
            ep = x
            break
    return ep, line, threshold


def ValidateROI(cs, roipts):
    # check if roi valid 

    widthpx, heightpx = np.shape(cs.pixeldataIn)
    # if any of the found corner points is outside the image, report error
    for (xa,ya) in roipts:
        if xa<0 or ya<0:
            print('ROI outside image: ', xa, ya)
            return False
        elif xa>=widthpx or ya>=heightpx:
            print('ROI outside image: ', xa, ya)
            return False
    return True

def _FineTunePhantomBox(cs, roipts):
    """
    Given an approximate phantom box, fine tune location of corner points to locations of grid crossings.
    Will try ScaleSpace Hessian extrema finding, and on fail plain minimum finding. 
    """
    # true to properly locate the cornerpoints of the 90mm box
    error = False
    roipts0 = copy.deepcopy(roipts)

    threshold = 0.65 # threshold for acceptance
    # using the determinant of the hessian gives nicest solution.
    # if it fails, try the dumb max finder
    best = (0, copy.deepcopy(roipts0))
    for radmm in [9.0, 6.0, 4.5]:
        error = False
        roipts = copy.deepcopy(roipts0)
        bbox_confidence, roipts = BBAlignROI(cs, roipts, radmm=radmm, useHessian=True)
        if bbox_confidence< threshold or ValidateROI(cs, roipts) == False:
            print('BBAlignROI failed! confidence:', bbox_confidence)
            error = True
            print('BBAlignROI try smaller searchradius')
        if bbox_confidence>best[0]:
            best = (bbox_confidence, copy.deepcopy(roipts)) 

    
    if error:
        print('BBAlignROI trying without Hessian')
        roipts = copy.deepcopy(roipts0)
        error = False
        bbox_confidence, roipts = BBAlignROI(cs, roipts, radmm=6.0, useHessian=False)

        if bbox_confidence< threshold or ValidateROI(cs, roipts) == False:
            print('BBAlignROI failed! confidence:', bbox_confidence)
            error = True
        if bbox_confidence>best[0]:
            best = (bbox_confidence, copy.deepcopy(roipts)) 

    bbox_confidence, roipts = best
    print('Using box with conf:', bbox_confidence)

    return error, roipts, bbox_confidence

def BBAlignROI(cs, roipts, radmm=9.0, useHessian=True):
    """
    Reposition corner points of roi to nearest likely position; check consistency of new roi; 
    calculate confidence in this new roi.

    """
    # candidate locations must be within 9 mm of starting location
    searchrad = int(cs.phantommm2pix(radmm)+.5)

    # for using blurring of sigma to get rid of noise
    sigma = max(1.5,searchrad/4.)

    # repeat repositioning ntries times
    ntries = 6

    print('BBAlignROI searchrad = ', searchrad)

    widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
    heightpx = np.shape(cs.pixeldataIn)[1]

    roipts0 = copy.deepcopy(roipts) # initial location needed for consistency check
    # find candidate positions with best confidence
    conf_pts = []
    conf_pts.append((0, copy.deepcopy(roipts))) # never select default pointset
    for itry in range(0, ntries):
        if itry >0:
            # adjust Laplace scale for Hessian of boundingbox; drill down in scale space
            sigma = max(cs.phantommm2pix(1),0.75*sigma) # never more accurate than 1 mm
            searchrad = int(4*sigma+1)

        for ipoint in range(0, len(roipts)):
            # restrict gaussian operations to small area around starting location,
            #  and checking boundaries of images
            rp = roipts[ipoint]
            x0 = int(rp[0])
            minx = max(0,x0-searchrad)
            maxx = min(widthpx-2,x0+searchrad)
            y0 = int(rp[1])
            miny = max(0,y0-searchrad)
            maxy = min(heightpx-2,y0+searchrad)

            # must work in float images else gaussians are ints as well, which is bad for poor contrast
            cropped = (cs.pixeldataIn[minx:maxx+1,miny:maxy+1]).astype(float)
            if useHessian == True:
                #Determinant of Hessian matrix 
                # maybe add requirement for minimum?
                Lxx = (scind.gaussian_filter(cropped,sigma,order=[2,0],mode='mirror'))
                Lyy = (scind.gaussian_filter(cropped,sigma,order=[0,2],mode='mirror'))
                Lxy = (scind.gaussian_filter(cropped,sigma,order=[1,1],mode='mirror'))
                cropped = Lxx*Lyy-Lxy**2
                Lxx = Lxy = Lyy = None

                # search for 2D max
                x1 = minx+np.unravel_index(cropped.argmax(), cropped.shape)[0]
                y1 = miny+np.unravel_index(cropped.argmax(), cropped.shape)[1]

            else: # sometimes Hessian fails to find good location; use following fall back
                # smoothing to get rid of noise, and subtract some average and multiply to enhance difference
                # remove background trend
                cropped = removeBKTrend(cropped, sigma)

                # remove noise
                cropped = scind.gaussian_filter(cropped, sigma)

                # search for min in projections
                xflat = np.mean(cropped,axis=0)
                yflat = np.mean(cropped,axis=1)
                x1 = minx+np.unravel_index(yflat.argmin(), yflat.shape)[0]
                y1 = miny+np.unravel_index(xflat.argmin(), xflat.shape)[0]

            # candidate location
            rp[0] = x1
            rp[1] = y1

        # force consistency with initial shape
        roiptsCA = BBConsistencyAlign(cs, roipts0, roipts)
        RCA = BBROIConfidence(cs, roiptsCA)
        RC  = BBROIConfidence(cs, roipts)
        if RCA> RC: # if consistency adjustment results in better confidence, accept it
            roipts = roiptsCA
            RC = RCA

        conf_pts.append((RC,copy.deepcopy(roipts)))

    # Just take best of the last two interations; prevent up down jumping
    #   generally sorting does worsen results somewhat if the starting point was already good
    if len(conf_pts)>1 and conf_pts[-2][0] > conf_pts[-1][0]:
        ci = -2
    else:
        ci = -1
    for i in range(0,len(roipts)):
        roipts[i][0] = conf_pts[ci][1][i][0]
        roipts[i][1] = conf_pts[ci][1][i][1]
    confidence = conf_pts[ci][0]

    return confidence, roipts

def BBConsistencyAlign(cs, roiptsold, roiptsnew):
    """
    4 types of motion:
        1. Shift (all points move in same direction)
        2. Scale (all distance to center change by same amount)
        3. Rotation (all points rotate along same angle around center)
        4. Shear (some points rotate more than others. Do not allow)

    Here just check if length of boxside agrees with given length. If not,
    try to find which point is most likely off. Try to fix that using the point
    which seems to be causing the least problems. Repeat twice.
    """
    midx0 = 0.
    midy0 = 0.
    for rp in roiptsold:
        midx0 += rp[0]
        midy0 += rp[1]
    midx0 /= len(roiptsold)
    midy0 /= len(roiptsold)

    shiftdist = []
    for rp0,rp1 in zip(roiptsold, roiptsnew):
        shiftdist.append(np.sqrt((rp1[0]-rp0[0])**2+(rp1[1]-rp0[1])**2))

    """
    [int(immidx-rad/2+.5),int(immidy-rad/2+.5)],
     [int(immidx-rad/2+.5),int(immidy+rad/2+.5)],
      [int(immidx+rad/2+.5),int(immidy+rad/2+.5)],
       [int(immidx+rad/2+.5),int(immidy-rad/2+.5)] ]

    0..3
    .  .
    .  .
    1..2
    """
    for irun in [0,1]:#,1]:
        #max line devs
        devs = [0,0,0,0]

        roiptsfix = copy.deepcopy(roiptsnew)
        idcmp = [ ([0,1], 1), ([3,2], 1), ([0,3], 0), ([1,2], 0) ] # [point0, point1], height/width

        # check if length of a side is not too different from what is expected
        for abc in idcmp:
            i0 = abc[0][0]
            i1 = abc[0][1]
            l1 = np.sqrt( (roiptsnew[i1][0]-roiptsnew[i0][0])**2+ (roiptsnew[i1][1]-roiptsnew[i0][1])**2)
            dev = np.abs(cs.pix2phantommm(l1)-cs.geom.box_radmm[abc[1]]*2)
            devs[i0] += dev*dev
            devs[i1] += dev*dev

        for i0,dev in enumerate(devs):
            if np.sqrt(dev)>5.: # this point causes to much deviation from expected value
                # use point that causes least troubles:
                bestpoint = 0
                for i1,d in enumerate(devs):
                    if d< devs[bestpoint]:
                        bestpoint = i1
                if np.sqrt(devs[bestpoint])>5.: 
                    return roiptsold # cannot find any good points

                [dx0, dy0] = [a - b for a,b in zip(roiptsold[bestpoint],[midx0, midy0])]
                len0 = np.sqrt(dx0*dx0+dy0*dy0)
                [dx1, dy1] = [a - b for a,b in zip(roiptsnew[bestpoint],[midx0, midy0])]
                len1 = np.sqrt(dx1*dx1+dy1*dy1)
                [dx, dy] = [a - b for a,b in zip(roiptsnew[bestpoint],roiptsold[bestpoint])]

                # translation of all points
                shiftx = (len1/len0-1.) * dx1
                shifty = (len1/len0-1.) * dy1

                # calc rotation angle
                da = np.arctan2(dy1-shifty, dx1-shiftx)-np.arctan2(dy0, dx0)

                a = np.arctan2(dy0,dx0)+da
                xx = shiftx+len0*np.cos(a)
                yy = shifty+len0*np.sin(a)

                # now apply to bad point
                [dxi, dyi] = [a - b for a,b in zip(roiptsold[i0],[midx0, midy0])]
                a = np.arctan2(dyi,dxi)+da
                leni = np.sqrt(dxi*dxi+dyi*dyi)

                xx = leni*np.cos(a)
                yy = leni*np.sin(a)

                roiptsfix[i0][0] = int(leni*np.cos(a)+midx0+shiftx)
                roiptsfix[i0][1] = int(leni*np.sin(a)+midy0+shifty)
                prev = roiptsnew[i0]
                repl = roiptsfix[i0]

        roiptsnew = roiptsfix
        shiftdist = []
        for rp0,rp1 in zip(roiptsold, roiptsnew):
            shiftdist.append(np.sqrt((rp1[0]-rp0[0])**2+(rp1[1]-rp0[1])**2))

    #sanity check
    if ValidateROI(cs, roiptsnew):
        return roiptsnew
    else:
        return roiptsold

def BBROIConfidence(cs, roipts):
    """
    Calculate confidence in found ROI: Check lenghts of sides, and diagonals, and compare 
    to given lengths, resulting in bad confidence if difference is more than 5 mm
    """
    confidence = 0.

    # First the lengths of just the sides of the box
    boxhwid = cs.phantommm2pix(cs.geom.box_radmm[0])
    boxhhei = cs.phantommm2pix(cs.geom.box_radmm[1])
    redlength = [boxhhei, boxhwid, boxhhei, boxhwid]
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

    """
    Calculate lengths and
            need to check only is NS^2 == len0^2+len1^2 and EW^2=len2^2+len3^2
            1cm dev must lead to fail. So lengths (12cm) minus 11cm. Scaling is from (sid+10)/sid.
            find magnification
            we want a cm dev to drop confidence to 0.5, so all lengths (emperically) need to be reduced by magnified 11cm
    """
    for (x0,y0) in roipts:
        if cs.pixeldataIn[int(x0)][int(y0)] == 0 or cs.pixeldataIn[int(x0)][int(y0)] == cs.max_pixel_value: # on annotation
            return 0

    confidence = 1.
    # ns via east+north
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
    for p,rl in enumerate(redlength):#range(0,4):
        confidence *= min(2*rl,lengths[p])/max(2*rl,lengths[p])
    confidence = confidence **6.

    print('BBConfidence =', (confidence*100.),'%')
    return confidence

def XRayField(cs):
    """
    Find edges of XRay exposure along the lines of the grid box.
    Use either min across line between box and edge, or use corner values
    """
    # remove noise; alternatively sample more pixels and average.
    workim = scind.gaussian_filter(cs.pixeldataIn, sigma=5)

    error = False

    xr_NSWEmm = []
    # north- and southside
    xr_NSWEmm.append(FindXRayEdge(cs, 'N', workim, rawim=cs.pixeldataIn))
    xr_NSWEmm.append(FindXRayEdge(cs, 'S', workim, rawim=cs.pixeldataIn))
    xr_NSWEmm.append(FindXRayEdge(cs, 'W', workim, rawim=cs.pixeldataIn))
    xr_NSWEmm.append(FindXRayEdge(cs, 'E', workim, rawim=cs.pixeldataIn))
    if min(xr_NSWEmm)<1.:
        error = True
    else:
        error = False
    print('Edge [N/S/W/E] cm = %.1f %.1f %.1f %.1f' % (xr_NSWEmm[0]/10., 
                                                       xr_NSWEmm[1]/10., 
                                                       xr_NSWEmm[2]/10., 
                                                       xr_NSWEmm[3]/10. ))
    xr_roi = [] #[ UL, LL, LR, UR]
    xco, yco = cs.geom.phantomposmm2pix(-xr_NSWEmm[2],  xr_NSWEmm[0])
    xr_roi.append([xco, yco])
    xco, yco = cs.geom.phantomposmm2pix( xr_NSWEmm[3],  xr_NSWEmm[0])
    xr_roi.append([xco, yco])
    xco, yco = cs.geom.phantomposmm2pix( xr_NSWEmm[3], -xr_NSWEmm[1])
    xr_roi.append([xco, yco])
    xco, yco = cs.geom.phantomposmm2pix(-xr_NSWEmm[2], -xr_NSWEmm[1])
    xr_roi.append([xco, yco])

    # copy to geom struct
    cs.geom.xr_roi = xr_roi
    cs.geom.xr_NSWEmm = xr_NSWEmm

    return error

def FindXRayEdge(cs, side, workim, rawim=None):
    # travel from center to edges
    # workim is a prepared image, getting rid of noise and other stuff
    widthpx, heightpx  = np.shape(workim) ## width/height in pixels

    outvalue = cs.forceRoom.outvalue
    # for DiDi, just take the minimal corner value
    if outvalue<0:
        outvalue = min(workim[0][0], workim[-1][0],workim[0][-1],workim[-1][-1])

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
            dypx = 1.*(cs.geom.box_roi[baseidmax][1]-cs.geom.box_roi[baseidmin][1])/(cs.geom.box_roi[baseidmax][0]-cs.geom.box_roi[baseidmin][0]) # diff in x if 1 px to y
            useboxradmm = cs.geom.box_radmm[0]
        else:
            dxpy = 1.*(cs.geom.box_roi[baseidmax][0]-cs.geom.box_roi[baseidmin][0])/(cs.geom.box_roi[baseidmax][1]-cs.geom.box_roi[baseidmin][1]) # diff in x if 1 px to y
            useboxradmm = cs.geom.box_radmm[1]
        posvec = []
        valvec = []
        rawvalvec = [] # on original, allow for shrinkage
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
                x0 = int(cs.geom.box_roi[baseidmax][0]) +int(xpos)
                y0 = int(cs.geom.box_roi[baseidmax][1]) +int(ypos)
            else:
                x0 = int(cs.geom.box_roi[baseidmin][0]) +int(xpos)
                y0 = int(cs.geom.box_roi[baseidmin][1]) +int(ypos)
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
                    if not rawim is None:
                        rawvalvec.append(int(rawim[xa,ya]))
                break

            if not rawim is None:
                val00 = int(rawim[x0,y0])
                val10 = int(rawim[x1,y0])
                val05 = 1.*val00+(xpos-(int)(xpos))*(val10-val00)
                rawvalvec.append(val05)

            val00 = int(workim[x0,y0])
            val10 = int(workim[x1,y0])
            val05 = 1.*val00+(xpos-(int)(xpos))*(val10-val00)
            posvec.append(pos)
            valvec.append(val05)
                
            if val00 == outvalue:
                break
            id += 1

        #minval = min(valvec)
        minval = min((min(valvec),outvalue))
        maxval = max((max(valvec),outvalue)) # outvalue might not be present in this direction
        meanval = np.mean(valvec)

        if meanval < outvalue: # do not include outvalue parts in mean
            npa = np.array(valvec)
            meanval = np.mean( npa[npa<outvalue] )
        elif meanval > outvalue:
            npa = np.array(valvec)
            meanval = np.mean( npa[npa>outvalue] )
            
        if outvalue < meanval: # looking for a low value
            threshLow = (9.*minval+meanval)/10.
            lab = "low"
        else:
            threshLow = (9.*maxval+meanval)/10.
            lab = "high"
        threshHigh = outvalue
            
        if cs.verbose:
            plt.figure()
            plt.plot(posvec,valvec)
            if not rawim is None:
                plt.plot(posvec,rawvalvec)
            plt.title(side+" "+str(threshLow)+" "+str(threshHigh)+" "+lab)
            cs.hasmadeplots = True

        found = False
        for ix, (p,v) in enumerate(zip(posvec,valvec)):
            if outvalue < meanval and v<threshLow:
                found = True
            elif outvalue >= meanval and v>threshLow:
                found = True
            if found:
                if not rawim is None: # try to exclude constant region outside which is now blurred.
                    usep = p
                    for _ix in reversed(range(0,ix)):
                        if rawvalvec[_ix] == rawvalvec[ix]:
                            usep = posvec[_ix]
                        else:
                            break
                    p = usep
                edgemm.append( cs.pix2phantommm(p)+useboxradmm )
                if cs.verbose:
                    plt.plot(p,v,'bo')
                    cs.hasmadeplots = True
                break

        if not found and threshHigh>threshLow:
            for ix, (p,v) in enumerate(zip(posvec,valvec)):
                if v>=threshHigh:
                    found = True
                    if not rawim is None: # try to exclude constant region outside which is now blurred.
                        usep = p
                        for _ix in reversed(range(0,ix)):
                            if rawvalvec[_ix] == rawvalvec[ix]:
                                usep = posvec[_ix]
                            else:
                                break
                        p = usep
                    edgemm.append( cs.pix2phantommm(p)+useboxradmm )
                    if cs.verbose:
                        plt.plot(p,v,'bo')
                        cs.hasmadeplots = True
                    break

        if not found: # add max edge pos
            edgemm.append( cs.pix2phantommm(max(posvec))+useboxradmm )


    return max(edgemm)

def MTFAlignROI(cs, roipts):
    """
    Reposition corner points of roi to nearest likely position; check consistency of new roi; 
    calculate confidence in this new roi.
    'MTF_'+cs.forceRoom.linepairmodel
    """
    # candidate locations must be within .5 mm of starting location
    searchrad = min(2,int(cs.phantommm2pix(0.5)+.5))

    # for using blurring of sigma to get rid of noise
    sigma = max(.8,searchrad/4.)

    # repeat repositioning ntries times
    ntries = 6

    print('MTFAlignROI searchrad = ', searchrad)

    widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
    heightpx = np.shape(cs.pixeldataIn)[1]

    roipts0 = copy.deepcopy(roipts) # initial location needed for consistency check
    # find candidate positions with best confidence
    conf_pts = []
    conf_pts.append((0, copy.deepcopy(roipts), 0)) # never select default pointset
    for itry in range(0, ntries):
        for ipoint in range(0, len(roipts)):
            # restrict gaussian operations to small area around starting location,
            #  and checking boundaries of images
            rp = roipts[ipoint]
            x0 = int(rp[0])
            minx = max(0,x0-searchrad)
            maxx = min(widthpx-2,x0+searchrad)
            y0 = int(rp[1])
            miny = max(0,y0-searchrad)
            maxy = min(heightpx-2,y0+searchrad)

            # must work in float images else gaussians are ints as well, which is bad for poor contrast
            cropped = (cs.pixeldataIn[minx:maxx+1,miny:maxy+1]).astype(float)

            # smoothing to get rid of noise, and subtract some average and multiply to enhance difference
            # remove background trend
            cropped = removeBKTrend(cropped, sigma)

            # remove noise
            cropped = scind.gaussian_filter(cropped, sigma)

            # search for min in projections
            xflat = np.mean(cropped,axis=0)
            yflat = np.mean(cropped,axis=1)
            if cs.forceRoom.linepairmodel == 'RXT02': # search for max
                x1 = minx+np.unravel_index(yflat.argmax(), yflat.shape)[0]
                y1 = miny+np.unravel_index(xflat.argmax(), xflat.shape)[0]
            elif cs.forceRoom.linepairmodel == 'typ38': # search for min
                x1 = minx+np.unravel_index(yflat.argmin(), yflat.shape)[0]
                y1 = miny+np.unravel_index(xflat.argmin(), xflat.shape)[0]
            else:
                raise ValueError('Unknown linepairmodel %s'%cs.forceRoom.linepairmodel)

            if cs.verbose:
                plt.figure()
                plt.imshow(cropped)
                plt.title('MTFAlign '+str(itry)+ ' '+str(ipoint))
                print(sigma,"shift ",itry," of point ",ipoint," =(",x1-x0,",",y1-y0,")")
                cs.hasmadeplots = True

            # candidate location
            rp[0] = x1
            rp[1] = y1

        # force consistency with initial shape
        roiptsCA, adjustmtftangledeg = MTFConsistencyAlign(cs, roipts0, roipts)
        RCA = MTFROIConfidence(cs, roiptsCA)
        RC  = MTFROIConfidence(cs, roipts)
        if RCA> RC: # if consistency adjustment results in better confidence, accept it
            roipts = roiptsCA
            RC = RCA
        else:
            adjustmtftangledeg = 0
        conf_pts.append((RC,copy.deepcopy(roipts), adjustmtftangledeg))

    # for MTF sorting seems to help:
    conf_pts = sorted(conf_pts)
    for i in range(0,len(roipts)):
        roipts[i][0] = conf_pts[-1][1][i][0]
        roipts[i][1] = conf_pts[-1][1][i][1]
    confidence = conf_pts[-1][0]
    adjustmtftangledeg = conf_pts[-1][2]

    return confidence, roipts, adjustmtftangledeg

def MTFConsistencyAlign(cs, roiptsold, roiptsnew):
    """
    4 types of motion:
        1. Shift (all points move in same direction)
        2. Scale (all distance to center change by same amount)
        3. Rotation (all points rotate along same angle around center)
        4. Shear (some points rotate more than others. Do not allow)

    Here just check if length of boxside agrees with given length. If not,
    try to find which point is most likely off. Try to fix that using the point
    which seems to be causing the least problems. Repeat twice.
    """
    midx0 = 0.
    midy0 = 0.
    for rp in roiptsold:
        midx0 += rp[0]
        midy0 += rp[1]
    midx0 /= len(roiptsold)
    midy0 /= len(roiptsold)

    shiftdist = []
    for rp0,rp1 in zip(roiptsold, roiptsnew):
        shiftdist.append(np.sqrt((rp1[0]-rp0[0])**2+(rp1[1]-rp0[1])**2))

    if cs.forceRoom.linepairmodel == 'typ38':
        """
        roipts = [ [x18px,y18px],[x06px,y06px],[x14px,y14px],[x46px,y46px] ]
        """
        adjustmtfangledeg = 0.
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
            if(np.abs(cs.pix2phantommm(l1)-l0mm)>.45):
                alter = i1
                if np.abs(np.median(shiftdist)-shiftdist[i0]) > np.abs(np.median(shiftdist)-shiftdist[i1]):
                #if(shiftdist[i0]>shiftdist[i1]):
                    alter = i0
                if (alter == id18): #18:
                    roiptsnew[alter][0] = int(0.5+roiptsnew[id46][0]-len18_46/len06_14*dd06_14[0])
                    roiptsnew[alter][1] = int(0.5+roiptsnew[id46][1]-len18_46/len06_14*dd06_14[1])
                    adjustmtfangledeg = 0.5
                elif(alter == id06): #06
                    roiptsnew[alter][0] = int(0.5+roiptsnew[id14][0]-len06_14/len18_46*dd18_46[0])
                    roiptsnew[alter][1] = int(0.5+roiptsnew[id14][1]-len06_14/len18_46*dd18_46[1])
                elif(alter == id14): #14
                    roiptsnew[alter][0] = int(0.5+roiptsnew[id06][0]+len06_14/len18_46*dd18_46[0])
                    roiptsnew[alter][1] = int(0.5+roiptsnew[id06][1]+len06_14/len18_46*dd18_46[1])
                elif(alter == id46): #46
                    roiptsnew[alter][0] = int(0.5+roiptsnew[id18][0]+len18_46/len06_14*dd06_14[0])
                    roiptsnew[alter][1] = int(0.5+roiptsnew[id18][1]+len18_46/len06_14*dd06_14[1])
                    adjustmtfangledeg = 0.5
        for abc,l0mm in zip(idcmp,lengthmms):
            i0 = abc[0]
            i1 = abc[1]
            l1 = np.sqrt( (roiptsnew[i1][0]-roiptsnew[i0][0])**2+ (roiptsnew[i1][1]-roiptsnew[i0][1])**2)

    if cs.forceRoom.linepairmodel == 'RXT02':
        adjustmtfangledeg = 0.

    #sanity check
    if ValidateROI(cs, roiptsnew):
        return roiptsnew, adjustmtfangledeg
    else:
        return roiptsold, adjustmtfangledeg

def MTFROIConfidence(cs, roipts):
    """
    Calculate confidence in found ROI: Check lenghts of sides, and diagonals, and compare 
    to given lengths, resulting in bad confidence if difference is more than 5 mm
    """
    confidence = 0.

    if cs.forceRoom.linepairmodel == 'RXT02':
        l0mm = 23.32
        dx = roipts[0][0]-roipts[1][0]
        dy = roipts[0][1]-roipts[1][1]
        l1mm = cs.pix2phantommm(np.sqrt(dx*dx+dy*dy))

        return min(l1mm,l0mm)/max(l1mm,l0mm)


    elif cs.forceRoom.linepairmodel == 'typ38':
        # First the lengths of just the sides of the box
        boxhwid = cs.phantommm2pix(cs.geom.box_radmm[0])
        boxhhei = cs.phantommm2pix(cs.geom.box_radmm[1])
        redlength = [boxhhei, boxhwid, boxhhei, boxhwid]
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
    else:
        raise ValueError('Unknown linepair model %s'%cs.forceRoom.linepairmodel)

    print('MTFConfidence =', (confidence*100.),'%')
    return confidence

