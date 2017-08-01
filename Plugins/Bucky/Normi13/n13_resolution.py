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
Normi13 analysis, resolution functions:
  o CTF
  o MTF
  
Horrible code; need to clean up. Lot's of unclear, redundant stuff

Changelog:
 20170731: MTF cutout rxt02 extended by 2.0 mm instead of 1.5 mm to fix no pattern found lowest freq in low res images; 
           MTF bar analysis shift shift half a bar two sides (instead of one)
"""

import numpy as np
import scipy.ndimage as scind
from scipy import stats
import copy
import matplotlib.pyplot as plt
try:
    # wad2.0 runs each module stand alone
    from n13_geometry import MTFAlignROI, ValidateROI
    import n13_math as mymath
except ImportError:
    from .n13_geometry import MTFAlignROI, ValidateROI
    from . import n13_math as mymath

class MTFStruct:
    def __init__ (self):
        self.crLimit   = 0.1 # contrast limit for MTF
        self.sigma_ext = 1.5 # gaussian blur factor for extreme finder

        self.bShowMTFDetail = False
        self.bShowCTF = True
        self.bIgnoreMTFError = False # in rare cases MTF is "ok" but looks wrong

        self.roi = []
        self.pos_confidence = -1.

        self.contrast_response = []
        self.contrast_high = []
        self.contrast_low = []
        self.ctfmtf = []
        self.contrast_freqs = []
        self.contrast_tops = []
        self.contrast_bots = []
        self.mtf_aapm = -1.
        self.freq_confidence = -1.
        self.calc_freqs = []  # frequencies in phantom units as measured from extremes


def MTF(cs):
    """
    Calculate MTF from CTF from line pairs module in Normi13
    """
    error = False

    if cs.forceRoom.linepairmodel == 'RXT02':
        x06px,y06px = cs.geom.phantomposmm2pix(cs.forceRoom.xy06mm[0], cs.forceRoom.xy06mm[1])
        x10px,y10px = cs.geom.phantomposmm2pix(cs.forceRoom.xy10mm[0], cs.forceRoom.xy10mm[1])
        roipts = [ 
            [x06px,y06px],
            [x10px,y10px]
        ]
    elif cs.forceRoom.linepairmodel == 'typ38':
        x18px,y18px = cs.geom.phantomposmm2pix(cs.forceRoom.xy18mm[0], cs.forceRoom.xy18mm[1])
        x06px,y06px = cs.geom.phantomposmm2pix(cs.forceRoom.xy06mm[0], cs.forceRoom.xy06mm[1])
        x46px,y46px = cs.geom.phantomposmm2pix(cs.forceRoom.xy46mm[0], cs.forceRoom.xy46mm[1])
        x14px,y14px = cs.geom.phantomposmm2pix(cs.forceRoom.xy14mm[0], cs.forceRoom.xy14mm[1])
        roipts = [ 
            [x18px,y18px],
            [x06px,y06px],
            [x14px,y14px],
            [x46px,y46px] ]
    else:
        raise ValueError('[MTF] Unknown linepairmodel %s'%cs.forceRoom.linepairmodel)

    pos_confidence, roipts, adjustmtftangledeg = MTFAlignROI(cs, roipts)
    threshold = 0.65 # threshold for acceptance

    if pos_confidence< threshold or ValidateROI(cs, roipts) == False:
        print('MTFAlignROI failed! confidence:', pos_confidence)
        error = True

    # copy to mtf struct
    cs.mtf.roi = roipts
    cs.mtf.pos_confidence = pos_confidence

    if not error:
        error = AnalyseMTF(cs, adjustmtftangledeg)
    return error

def AnalyseMTF(cs, adjustmtftangledeg=0):
    if cs.forceRoom.linepairmodel == 'RXT02':
        return _AnalyseMTF_RXT02(cs)
    elif cs.forceRoom.linepairmodel == 'typ38':
        return _AnalyseMTF_typ38(cs, adjustmtftangledeg)
    else:
        raise ValueError('[AnalyseMTF] Unknown linepairmodel %s'%cs.forceRoom.linepairmodel)

def _AnalyseMTF_typ38(cs, adjustmtfangledeg=0):
    """
    For linepairs insert model typ38.
    Find rotation angle of linepairs insert. Cut out the back rotated insert.
    For each linepair, do an analysis
    """
    error = True

    # first extend roi a little
    extend18 = cs.phantommm2pix(2.8)   # extend bbox beyond dot in '1.8' [mm]
    extend46 = cs.phantommm2pix(3.2)   # extend beyond dot in '4.6' [mm]
    print("2.8",extend18)
    print("3.2",extend46)
    # First cut out rotated roi
    id18 = 0
    id06 = 1
    id14 = 2
    id46 = 3
    len1846 = np.sqrt((cs.mtf.roi[id18][0]-cs.mtf.roi[id46][0])**2+(cs.mtf.roi[id18][1]-cs.mtf.roi[id46][1])**2)
    print("1846=",len1846)
    extend18 = 0.0786*len1846
    extend46 = 0.0898*len1846
    copyimage = cs.pixeldataIn.astype(float)
    rotanglerad = 3.*np.pi/2.-.5*(np.arctan2((cs.mtf.roi[id18][1]-cs.mtf.roi[id46][1]),(cs.mtf.roi[id18][0]-cs.mtf.roi[id46][0]))+np.arctan2((cs.mtf.roi[id06][1]-cs.mtf.roi[id14][1]),(cs.mtf.roi[id06][0]-cs.mtf.roi[id14][0])))
    rotanglerad += adjustmtfangledeg/180*np.pi # emperically found adjustment
    rotangledeg = (rotanglerad/np.pi*180.)
    print("MTF at",rotangledeg, "degrees")
    rotimage = scind.interpolation.rotate(copyimage, rotangledeg, axes=(1, 0), reshape=False, output=None, order=3, mode='constant', cval=0.0, prefilter=True)

    # rotate roi
    costerm = np.cos(rotanglerad)
    sinterm = np.sin(rotanglerad)
    xc = cs.pixeldataIn.shape[0]/2.
    yc = cs.pixeldataIn.shape[1]/2.
    roipts = []
    for rp in cs.mtf.roi:
        xp = xc +(rp[0]-xc)*costerm-(rp[1]-yc)*sinterm
        yp = yc +(rp[0]-xc)*sinterm+(rp[1]-yc)*costerm
        roipts.append([xp,yp])

    # extend roi
    roipts[id18][1] -= extend18
    roipts[id46][1] += extend46
    minxco = roipts[0][0]
    minyco = roipts[0][1]
    maxxco = roipts[0][0]
    maxyco = roipts[0][1]
    for rp in roipts:
        minxco = int(min(minxco,rp[0]))
        maxxco = int(max(maxxco,rp[0]))
        minyco = int(min(minyco,rp[1]))
        maxyco = int(max(maxyco,rp[1]))

    smallimage = rotimage[minxco:maxxco+1,minyco:maxyco+1]

    # copy to struct
    cs.lastimage = smallimage # for gui if error
    
    return _MTF_smallimage(cs, smallimage)

def _AnalyseMTF_RXT02(cs):
    """
    For linepairs insert model RXT02.
    Find rotation angle of linepairs insert. Cut out the back rotated insert.
    For each linepair, do an analysis
    """
    error = True

    roipts = copy.deepcopy(cs.mtf.roi)

    # need to extend roi
    extendup = cs.phantommm2pix(2.0)  # extend bbox beyond dot in '1.8' [mm]
    extenddo = cs.phantommm2pix(15.)  # extend beyond dot in '4.6' [mm]

    # extend starting points
    dy = roipts[1][1]-roipts[0][1]
    dx = roipts[1][0]-roipts[0][0]

    l0 = np.sqrt(dx*dx+dy*dy)
    roipts[0][0] -= extendup/l0*dx
    roipts[0][1] -= extendup/l0*dy
    roipts[1][0] += extenddo/l0*dx
    roipts[1][1] += extenddo/l0*dy

    #append points
    l1 = cs.phantommm2pix(35.)
    px = roipts[1][0] - l1/l0*dy
    py = roipts[1][1] + l1/l0*dx
    roipts.append([px,py])

    px = roipts[0][0] - l1/l0*dy
    py = roipts[0][1] + l1/l0*dx
    roipts.append([px,py])

    # angle
    rotanglerad = np.arctan2(dx,dy)

    # First cut out rotated roi
    copyimage = cs.pixeldataIn.astype(float)
    rotangledeg = (rotanglerad/np.pi*180.)
    print("MTF at",rotangledeg, "degrees")
    rotimage = scind.interpolation.rotate(copyimage, rotangledeg, axes=(1, 0), reshape=False, output=None, order=3, mode='constant', cval=0.0, prefilter=True)

    costerm = np.cos(rotanglerad)
    sinterm = np.sin(rotanglerad)
    xc = cs.pixeldataIn.shape[0]/2.
    yc = cs.pixeldataIn.shape[1]/2.

    new_roipts = []
    for rp in roipts:
        xp = int(xc +(rp[0]-xc)*costerm-(rp[1]-yc)*sinterm)
        yp = int(yc +(rp[0]-xc)*sinterm+(rp[1]-yc)*costerm)
        new_roipts.append([xp,yp])

    minxco = new_roipts[0][0]
    minyco = new_roipts[0][1]
    maxxco = new_roipts[0][0]
    maxyco = new_roipts[0][1]
    for rp in new_roipts:
        minxco = min(minxco,rp[0])
        maxxco = max(maxxco,rp[0])
        minyco = min(minyco,rp[1])
        maxyco = max(maxyco,rp[1])

    smallimage = rotimage[minxco:maxxco+1,minyco:maxyco+1]
    cs.lastimage = smallimage

    # copy to mtf struct
    cs.mtf.roi = roipts

    return _MTF_smallimage(cs, smallimage)

def _MTF_smallimage(cs, smallimage):
    """
    Analysis driver of image cropped to linepairs insert
    """
    error = True
    if cs.mtf.bShowMTFDetail or cs.verbose:
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
    start_endpos = FillMTFBarDetails(cs, smallimage)
    
    # Find high contrast response of line patterns
    for vpi in range(1,num_freq):
        contrast_response[vpi] = 0.
        contrast_high[vpi]     = 0.
        contrast_low[vpi]      = 0.
        contrast_tops[vpi]     = 0
        contrast_bots[vpi]     = 0
        calc_freq[vpi]         = 0.
        contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = AnalyseMTF_Part(cs,smallimage, start_endpos[vpi-1], vpi)

    if contrast_tops[0] == 0:
        contrast_response[0] = contrast_response[1]
        contrast_high[0] = contrast_high[1]
        contrast_low[0] = contrast_low[1]

    # Correct for Nyquist
    #fnyq = (0.5/self.dcmInfile.PixelSpacing[0])
    fnyq = 0.5/cs.pix2phantommm(1.)
    for vpi in range(0,num_freq):
        if contrast_high[vpi]+contrast_low[vpi]>1e-6:
            contrast_response[vpi] = (contrast_high[vpi]-contrast_low[vpi])/(contrast_high[vpi]+contrast_low[vpi])
        else:
            contrast_response[vpi] = 0.

        #check if not beyond Nyquist
        if contrast_freqs[vpi]>fnyq:
            contrast_response[vpi] = 0

    # Correct for consistency
    crLimitabs = cs.mtf.crLimit*contrast_response[0]
    for vpi in range(1,num_freq-1):
        if(contrast_response[vpi]<crLimitabs and contrast_response[vpi-1]>crLimitabs and contrast_response[vpi+1]>crLimitabs):
            contrast_response[vpi] = (contrast_response[vpi-1] + contrast_response[vpi+1])/2.

    if contrast_response[0]<1.e-6:
        print("Error in MTF: Rotated image?")
        return error

    ctfmtf = CTFtoMTF(cs, contrast_freqs,contrast_response)

    if cs.verbose or cs.mtf.bShowCTF:
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
    mtf_aapm = (0.5/cs.pix2phantommm(1.))*.9 # should ignore acos, because I calculate after rotation, so there

    # Test how well calculated frequencies match with given frequencies
    maxid = len(calc_freq)
    for id in reversed(range(2,len(calc_freq))):
        #if( (calc_freq[id]<1e-6 or np.abs(calc_freq[id]- calc_freq[id-1])<1.e-6) ):
        if calc_freq[id]<1e-6:
            maxid = id
    print("maxid:",maxid)
    if maxid<5 and not cs.mtf.bIgnoreMTFError:
        print("Error in MTF: Rotated image?")
        return error
    slope, intercept, r_value, p_value, std_err = stats.linregress(contrast_freqs[0:maxid],calc_freq[0:maxid])
    if r_value**2<0.7:
        print("maxid:",maxid)
        for co,ca in zip(contrast_freqs,calc_freq):
            print(co,ca)
    # To get coefficient of determination (r_squared)
    #print("slope:",slope)
    #print("intercept:",intercept)
    #print("maxid:",maxid)
    #print("r-squared:", r_value**2)
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

        #print("confid:",mtf_found_confidence)
        
    # copy to mtf struct
    cs.mtf.mtf_aapm          = mtf_aapm
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

def FillMTFBarDetails(cs, smallimage):
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

    elif cs.forceRoom.linepairmodel == 'typ38'  :
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
    else:
        raise ValueError('[FillMTFBarDetails] Unknown linepairmodel')

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

def AnalyseMTF_Part(cs,smallimage, start_endpos, vpi):
    """
    Determine contrast in one bar pattern element
        contrast_response[vpi],contrast_high[vpi],contrast_low[vpi],contrast_tops[vpi],contrast_bots[vpi],calc_freq[vpi] = self.AnalyseMTF_Part(smallimage, startpos[vpi-1],endpos[vpi-1], vpi)
    Already inverted smallimage as input
    """
    contrast_response = 0.
    contrast_high     = 0.
    contrast_low      = 0.
    contrast_tops     = 0
    contrast_bots     = 0
    calc_freq = 0

    startpos = start_endpos[0]
    endpos   = start_endpos[1]

    sw,wh = smallimage.shape
    offset = max(1, int((endpos[1]-startpos[1])/8.+.5)) # about half a linepair
    for oi in [(0,0), (0,1), (-1,0), (-1,1), (0,0) ]: # if it doesnt work, use 0
        #ipbar = smallimage[startpos[0]:endpos[0],max(0, startpos[1]+oi[0]*offset):min(endpos[1]+oi[1]*offset, wh)]
        ipbar = smallimage[max(0, startpos[0]+oi[0]*offset):min(endpos[0]+oi[1]*offset, wh),
                           max(0, startpos[1]+oi[0]*offset):min(endpos[1]+oi[1]*offset, wh)]
        print('[AnalyseMTF_Part]',vpi,startpos,endpos)
        
        pwid = ipbar.shape[0]
        phei = ipbar.shape[1]
    
        pattern = np.zeros(phei,dtype=float)
        for y in range(0,phei):
            for x in range(0,pwid):
                pattern[y]+= ipbar[x,y]
            pattern[y] /= pwid
        
        if pattern.shape[0]<2:
            print("[AnalyseMTF_Part] SKIP: no pattern left ({}:{}, {}:{})".format(startpos[0],endpos[0], startpos[1]+oi[0]*offset,endpos[1]+oi[1]*offset ))
            return contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq
    
        # 1. find abs min and max
        # 2. between (max+min)/2 point left and right
        # 3. find all tops and valleys
        mustplot,tops,bots = FindExtrema(cs, pattern)

        if len(tops) < 3 or len(bots)<2:
            print("[AnalyseMTF_Part] could not find 3 tops and 2 bots, trying again after shifting.")
            continue
        else:
            break

    # determine 100%
    if vpi == 1 and len(tops)>0:
        xmin = endpos[0]
        length = smallimage.shape[0]-xmin
        ytop = tops[0]+startpos[1]
        hsize = max( 1,int(.5+cs.phantommm2pix(0.75)/2.)  )
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
            cs.hasmadeplots = True

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
            calc_freq = .5/cs.pix2phantommm(halflambda)
    if cs.verbose:
        print("Found",len(tops)," tops and",len(bots)," bottoms. Contrast=",contrast_response)
    #print(vpi,contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq)
    return contrast_response,contrast_high,contrast_low,contrast_tops,contrast_bots,calc_freq

def FindExtrema(cs, pattern):
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

    sigma = cs.mtf.sigma_ext
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
    mustplot,tops,bots = FindAllExtrema(cs, pattern, xderiv1, yminmax)

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
        mustplot,tops,bots = FindAllExtrema(cs, pattern, xderiv1, yminmax)

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
        mustplot,tops,bots = FindAllExtrema(cs, pattern, xderiv1, yminmax)

        if(len(tops)>3 or len(bots)>2): # restore previous
            tops = copy.deepcopy(tops_bk)
            bots = copy.deepcopy(bots_bk)
            mustplot = mustplot_bk
            goon = True  # break from loop

    if cs.verbose:
        print("[FindExtrema]C ntops/nbots = ",len(tops),"/",len(bots))
    return mustplot,tops,bots

def FindAllExtrema(cs, pattern, xderiv1, yminmax):
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

def CTFtoMTF(cs, freq_in, ctf):
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
    fnyq = 0.5/cs.pix2phantommm(1.) # cut-off frequencies > Nyquist
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

def CTFtoMTFNoFit(freq_in, ctf):
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

