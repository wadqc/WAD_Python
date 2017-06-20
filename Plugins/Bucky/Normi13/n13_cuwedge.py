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
Normi13 analysis, Cu wedge functions:
  o Mean, SNR, STDEV
  o Dynamic Range (max step/min step)
"""
import numpy as np
import scipy.ndimage as scind
import matplotlib.pyplot as plt

class CuStruct:
    def __init__ (self):
        self.box_roi = [] # max horizontal box in complete edge

        # for each step: mmCu, roi, snr, cnr, mean, sdev
        self.step_mean = []
        self.step_sdev = []
        self.step_snr = []
        self.step_cnr = []
        self.step_mmcu = []

        self.wedge_confidence = -1.

        self.dynamicRange = -1
        self.step_rois = []

def CuWedge(cs):
    """
    Analyse the Cu wedge. First find all steps. Calculate Dynamic Range
    """
    error = True

    # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)

    bump = 2. # 2 mm offset from wedge edges
    xmmll = -63.  +bump
    ymmur =  65.  +bump
    xmmur =  62.  -bump
    ymmll =  82.  -bump

    xpxll,ypxll = cs.geom.phantomposmm2pix(xmmll, ymmll) # lowerlef
    xpxur,ypxur = cs.geom.phantomposmm2pix(xmmur, ymmur) # upperright
    xpxul,ypxul = cs.geom.phantomposmm2pix(xmmll, ymmur) # upperleft
    xpxlr,ypxlr = cs.geom.phantomposmm2pix(xmmur, ymmll) # lowerright

    # box needs to be horizontal, to minimize influence of Heel effect
    xlo = int(.5+ max(xpxll,xpxul))
    xhi = int(.5+ min(xpxlr,xpxur))
    ylo = int(.5+ max(ypxlr,ypxll))
    yhi = int(.5+ min(ypxur,ypxul))
    if ylo>yhi:
        print('[CuWedge]: Error, phantom angle too large, cannot sample wedge')
        return error

    # 1. Make box around wedge (-5.5; -4) to (+5.5; -5.5)
    box_roi = [ [xlo,yhi],[xhi,yhi],[xhi,ylo],[xlo,ylo] ]

    # copy to cu struct
    cs.cuwedge.box_roi = box_roi
    
    # real analysis
    error = _AnalyseWedge(cs)

    return error

def _AnalyseWedge(cs):
    """
    Concept:
        1. Cut out rectangular roi
        2. For each cu step, measure mean and sd
        3. Calculate SNR and CNR
        4. Calculate DynamicRange
    """
    error = True

    # mm Cu; according to DIN 6866-13; note this thickness includes backplate
    digi_13 = [2.30, 1.85, 1.40, 1.00, 0.65, 0.30, 0.00]   

    # 1. Cut out rectangular roi
    xmin = cs.cuwedge.box_roi[0][0]
    xmax = cs.cuwedge.box_roi[1][0]
    ymin = cs.cuwedge.box_roi[2][1]
    ymax = cs.cuwedge.box_roi[1][1]
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

    profile = scind.gaussian_filter1d(profile, sigma=3., order=1)
    if cs.verbose:
        plt.figure()
        plt.plot(profile)
        plt.title('Cu Wedge derivative')
        cs.hasmadeplots = True

    # 2.2 Find the edges between the steps
    n_edges = 6
    posedges = []
    flatpix = int( cs.phantommm2pix(3) +.5) # flat step at least 3mm wide
    mindist = 0.75*len(profile)/n_edges
    for ix in range(0,n_edges):
        accept = False
        while(not accept):
            """
            Look for the next biggest peak. 
            For simplicity, check if peak not too close to others (then prob. noisy)
            If acceptable, flatten part of profile
            """
            profminid = np.unravel_index(profile.argmax(), profile.shape)[0]
            miniy = max(0,profminid-flatpix)
            maxiy = min(wid-1,profminid+flatpix)
            flatval = .5*(profile[maxiy]+profile[miniy])
            for iy in range(miniy,maxiy+1):
                profile[iy] = flatval
                if len(posedges)>1:
                    mdist = np.min( [np.abs(profminid-pe) for pe in posedges ])
                    if mdist> mindist:
                        accept = True
                else:
                    accept = True

        posedges.append(profminid)

    posedges = sorted(posedges)

    # calculate a confidence by counting number of steps and comparing the sizes of the steps
    wedge_confidence = 1.
    wedge_confidence *= (1.*min(n_edges,len(posedges))/max(n_edges,len(posedges)))**3. # must be 6! 
    avg_dist = 1.*(posedges[-1]-posedges[0])/(len(posedges)-1)
    for ix in range(1,len(posedges)):
        dist = 1.*(posedges[ix]-posedges[ix-1])
        wedge_confidence *= min(avg_dist,dist)/max(avg_dist,dist)

    if cs.verbose:
        print("Edge 0 at ",posedges[0])
        for ix in range(1,n_edges):
            print("Edge ",ix," at ",posedges[ix]," sep= ",posedges[ix]-posedges[ix-1])

    # 2.3 Calculate statistics for each step
    xlo = 0
    ylo = 0   # flatpix
    yhi = hei # - flatpix

    step_rois = []
    step_mean = []
    step_sdev = []
    for p in posedges:
        xlo += flatpix
        xhi = p-flatpix
        value = np.mean(smallimage[xlo:xhi,ylo:yhi])
        step_mean.append( value )
        step_sdev.append( np.std(smallimage[xlo:xhi,ylo:yhi]) )
        roipts = [ [xlo+xmin,yhi+ymin],[xhi+xmin,yhi+ymin],[xhi+xmin,ylo+ymin],[xlo+xmin,ylo+ymin] ]
        step_rois.append(roipts)
        xlo = p
    xlo += flatpix
    xhi = wid-flatpix
    value = np.mean(smallimage[xlo:xhi,ylo:yhi])

    step_mean.append( value )
    step_sdev.append( np.std(smallimage[xlo:xhi,ylo:yhi]) )
    roipts = [ [xlo+xmin,yhi+ymin],[xhi+xmin,yhi+ymin],[xhi+xmin,ylo+ymin],[xlo+xmin,ylo+ymin] ]
    step_rois.append(roipts)

    step_mmcu = []
    step_snr = []
    step_cnr = []
    for ix in range(0,n_edges+1):
        step_snr.append( step_mean[ix]/step_sdev[ix] )
        #step_cnr.append( np.sqrt(2.*(step_mean[ix]-step_mean[n_edges])**2/(step_sdev[ix]**2+step_sdev[n_edges]**2) )
        if ix<n_edges:
            step_cnr.append( np.abs(step_mean[ix]-step_mean[ix+1])/np.sqrt(0.5*(step_sdev[ix]**2+step_sdev[ix+1]**2)) )
        else:
            step_cnr.append(0.)
        step_mmcu.append( digi_13[ix] )

    dynamicRange = max(step_mean[n_edges]/step_mean[0],step_mean[0]/step_mean[n_edges])

    if cs.verbose:
        print("mmCu","SNR","CNR")
        for m,s,c in zip(step_mmcu,step_snr,step_cnr):
            print(m,s,c)

    # copy to custruct
    cs.cuwedge.wedge_confidence = wedge_confidence

    cs.cuwedge.step_rois = step_rois
    cs.cuwedge.step_mean = step_mean
    cs.cuwedge.step_sdev = step_sdev
    cs.cuwedge.step_mmcu = step_mmcu
    cs.cuwedge.step_snr = step_snr
    cs.cuwedge.step_cnr = step_cnr
    cs.cuwedge.dynamicRange = dynamicRange

    error = False
    return error
