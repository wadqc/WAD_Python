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
CDMam_lib: classes and routines for calculation of CDMAM scores
Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

Some stuff that can be optimized:
 * choice for detection scale (radius around expectation locations)
 * choice for integration scale (sigma_bx in blurring/integration)
 * maybe exclude better influence of text next to boundary cells
 * maybe crop image to interesting size (before or first step of removal of grid)
 * speed up removal of grid (maybe with grid detection on lower scale?)

Changelog:
    20170802: added cdcom params; ensure cdcom.exe is executable
    20161220: removed class variables; removed testing stuff
    20160902: sync with wad2.0; Unified pywad1.0 and wad2.0;
    20150826: reshuffling in gridremove to demand less memory (windows cannot cope)
    20150213: first working version
"""
__version__ = '20170802'
__author__ = 'aschilham'

import subprocess
import os
import stat
import copy
import numpy as np
import scipy.ndimage as scind
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
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

try:
    # wad2.0 runs each module stand alone
    import CDMam_constants as lit
except:
    from . import CDMam_constants as lit

class CDMamPhantom:
    """
    Description of the CDMAM Phantom. For now version 3.2 and 3.4 are implemented
    """

    #
    #GT = 4 1
    #     3 2
    groundtruth = [ # for version 3.2 and 3.4
                    [0,3,2,1,4,3,2,1,2,3,0,0,0,0,0,0],#3.2
                    [1,1,2,3,2,3,3,1,2,1,4,0,0,0,0,0],#2.5
                    [1,3,1,1,4,4,1,4,3,3,4,2,0,0,0,0],#2
                    [2,4,3,4,3,2,1,3,2,3,1,4,2,0,0,0],#1.6
                    [1,1,3,1,2,1,4,2,4,2,2,3,3,1,0,0],#1.25
                    [4,2,1,4,2,1,4,1,3,4,1,2,2,4,1,0],#1
                    [1,4,2,4,4,2,1,3,2,2,3,1,3,1,4,2],#0.8
                    [2,1,1,2,3,3,4,1,4,4,3,1,2,1,4,1],#0.63
                    [3,3,1,4,1,4,2,2,3,2,4,3,2,3,1,2],#0.5
                    [0,2,4,2,1,3,3,1,1,4,3,3,4,4,1,1],#0.4
                    [0,0,1,3,4,2,4,3,2,2,3,4,3,3,4,4],#0.31
                    [0,0,0,4,1,4,2,2,4,4,2,4,2,1,3,1],#0.25
                    [0,0,0,0,2,4,4,3,4,2,3,1,2,4,2,1],#0.2
                    [0,0,0,0,0,2,1,2,3,2,3,4,3,4,3,4],#0.16
                    [0,0,0,0,0,0,3,2,3,3,1,1,4,3,3,2],#0.13
                    [0,0,0,0,0,0,0,4,2,2,4,2,2,3,4,0],#0.1
                    ]
    def __init__ (self, version):
        """
        Constructor of CDMamPhantom. 
        Difference between 3.2 and 3.4 is thickness/diammeter of gold discs
        """
        if version == "3.2":
            self.version = version
            self.thickness_um = [0.05, 0.06, 0.08, 0.10,0.13, 0.16, 0.20, 0.25,
                                 0.31, 0.40, 0.50, 0.63, 0.80, 1.00, 1.25, 1.60]
            self.diameter_mm = [3.20, 2.50, 2.00, 1.60,
                                1.25, 1.00, 0.80, 0.63,
                                0.50, 0.40, 0.31, 0.25,
                                0.20, 0.16, 0.13, 0.10]
            # fit: factor = 2.13036+0.195548*ln(diam-0.0573261)
            self.humanfactor = [
                2.35,2.31,2.26,2.22,2.16, # values from log fit to published values
                2.11,2.06,2.01,1.98,1.94,1.88,1.82,1.75,1.68,1.6,1.5] # published values for 1.00 to 0.08

            self.phantomlim = [ # diam_mm minthick maxthick
                                [3.20, 0.06, 0.40],
                                [2.50, 0.05, 0.50],
                                [2.00, 0.05, 0.63],
                                [1.60, 0.05, 0.80],
                                [1.25, 0.05, 1.00],
                                [1.00, 0.05, 1.25],
                                [0.80, 0.05, 1.60],
                                [0.63, 0.05, 1.60],
                                [0.50, 0.05, 1.60],
                                [0.40, 0.06, 1.60],
                                [0.31, 0.08, 1.60],
                                [0.25, 0.10, 1.60],
                                [0.20, 0.13, 1.60],
                                [0.16, 0.16, 1.60],
                                [0.13, 0.20, 1.60],
                                [0.10, 0.25, 1.25],
                                ]
        if version == "3.4":
            self.version = version
            self.thickness_um = [0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.13, 0.16,
                                 0.20, 0.25, 0.36, 0.50, 0.71, 1.00, 1.42, 2.00]
            self.diameter_mm = [2.00, 1.60, 1.25, 1.00, 0.80, 0.63, 0.50, 0.40,
                                0.31, 0.25, 0.20, 0.16, 0.13, 0.10, 0.08, 0.06]
            self.humanfactor = [
                2.26,2.22,2.16, # values from log fit to published values
                2.11,2.06,2.01,1.98,1.94,1.88,1.82,1.75,1.68,1.6,1.5,1.4, # published values for 1.00 to 0.08
                0.97] # values from log fit to published values
            self.phantomlim = [ # diam_mm minthick maxthick
                                [2.00, 0.04, 0.25],
                                [1.60, 0.03, 0.36],
                                [1.25, 0.03, 0.50],
                                [1.00, 0.03, 0.71],
                                [0.80, 0.03, 1.00],
                                [0.63, 0.03, 1.42],
                                [0.50, 0.03, 2.00],
                                [0.40, 0.03, 2.00],
                                [0.31, 0.03, 2.00],
                                [0.25, 0.04, 2.00],
                                [0.20, 0.05, 2.00],
                                [0.16, 0.06, 2.00],
                                [0.13, 0.08, 2.00],
                                [0.10, 0.10, 2.00],
                                [0.08, 0.13, 2.00],
                                [0.06, 0.16, 1.42],
                                ]

class CDMamStruct:
    """
    Class holding all information for IO between gui and Lib
    """
    def __init__ (self,dcmInfile,pixeldataIn,phantomversion):
        self.verbose = False

        # input image
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.imageFileName = None
        self.parsCDCOM = [] # optional extra runtime flags ('c','high')

        # for matlib plotting
        self.hasmadeplots = False

        self.phantom = CDMamPhantom(phantomversion)
        self.gridimage = None
        self.startingroi = [] # starting box in phantom to find grid coords
        self.gridrois = [] # cornerpoints of each grid cell
        self.plotcount = 0
        self.imnum = 0 # increasing number for id generation
        
        # main results
        self.diam_mm = [] # list of diameters in mm
        self.limit_um = [] # list of predicted human detectable thickness thresholds in um
        self.iqf = 0. # image quality factor
        self.threshold_fnames = [] # images of threshold_thickness
        self.fit_fnames = [] # images of fits

        # identifiers
        self.filtername = lit.stUnknown
        self.scannername = lit.stUnknown
        self.energypresentation = ""

class CDMam():
    """
    Class for calculation of CDMAM score
    """
    def __init__(self,guimode=False):
        self.qcversion = __version__
        self.guimode = guimode
        pass

    def _smoothHitMatrix(self,cs,hitmatrix):
        """
        Smoothing: kernel = 0.25 0.50 0.25
                            0.50 1.00 0.50
                            0.25 0.50 0.25
        and normalization = sum(kernel_used)
        at borders, just ignore, so north edge does not include top line of kernel
        """
    
        wid = len(hitmatrix[0])
        hei = len(hitmatrix)

        # 1. Smooth found values
        kernel = [
            [0.25, 0.50, 0.25],
            [0.50, 1.00, 0.50],
            [0.25, 0.50, 0.25]
        ]
        smoothhit = copy.deepcopy(hitmatrix)
        for y in range(hei):
            for x in range(wid):
                norm = 0
                val =0
                for yi in range(-1,2):
                    for xi in range(-1,2):
                        yy = y +yi
                        xx = x +xi
                        if(yy<0 or xx <0 or yy>(hei-1) or xx>(wid-1)):
                            continue
                        if(cs.phantom.groundtruth[yy][xx] == 0):
                            continue
                        val   += kernel[yi+1][xi+1]*hitmatrix[yy][xx]
                        norm  += kernel[yi+1][xi+1]
                if(norm>0):
                    smoothhit[y][x] = val/norm
        return smoothhit

    def _testline(self,line):
        """
        Test for changes from 0 to 1, return list of locations
        """
        left = []
        if max(line) == 0:
            return left
        for i in range(1,len(line)):
            if line[i]>line[i-1]:
                left.append(i)
        return left

    def _sigmoid(self, x, x0, k):
        """
        sigmoid function for fitting
        """
        y = 1. / (1. + np.exp(-k*(x-x0)))
        return y

    def _invsigmoid(self, y, x0, k):
        """
        inverse sigmoid function for evaluation
        """
        x = x0-np.log(1./y-1.)/k
        return x

    def thresholdThickness(self,cs,hitmatrix):
        """
        Calculation of threshold thickness
        """
        phantom    = cs.phantom

        diam_mm    = phantom.diameter_mm
        thick_um   = phantom.thickness_um
        gtmatrix   = phantom.groundtruth
        phantomlim = phantom.phantomlim

        wid = len(hitmatrix[0])
        hei = len(hitmatrix)
        dosmooth = True
        dofit    = True
        #  0 : not defined or no hit
        #  1 : hit
        #  0.5 : switched
        limit = []
        diams = []
        human = []
        if not dosmooth:
            for y in range(hei):
                for x in range(wid):
                    if (hitmatrix[y][x] > .45):
                        limit.append(thick_um[x])
                        diams.append(diam_mm[y])
                        human.append(phantom.humanfactor[y]*thick_um[x])
                        break
        else:
            smoothhit = self._smoothHitMatrix(cs,hitmatrix)

            # 2. show found values
            if cs.verbose:
                for y in range(hei):
                    for x in range(wid):
                        if (gtmatrix[y][x]==0):
                            print("", end=" ")
                        else:
                            print("%0.2f" % smoothhit[y][x], end=" ")
                    print("")
                print("")
  
            if(dofit == False):
                """
                direct threshold finding
                """
                for y in range(hei):
                    found = False
                    for x in range(wid):
                        if (smoothhit[y][x] > 0.62):
                            found = True
                            limit.append(thick_um[x])
                            break
                    if(found == False):
                        limit.append(thick_um[-1]+1.)
                    human.append(phantom.humanfactor[y]*limit[-1])
  
            else:
                """
                fit sigmoid threshold finding
                """
                thresh = 0.625
                minprop = 0.1
                maxprop = 0.999
  
                # 3. fill arrays of with acceptable datapoints for each diameter
                propsarr  = []
                thicksarr = []
                fitsarr   = []
                idxarr    = [] # to make sure we know to which diam it belongs
                fitidxs   = []
                for y in range(hei): # loop over diams
                    props  = []
                    thicks = []
                    for x in range(wid): # loop over thickness
                        if(gtmatrix[y][x] == 0): # only points that are in grid
                            continue
                        if(len(props) == 0):
                            if(smoothhit[y][x]>=thresh): # if first point is already larger, stop searching
                                props.append(smoothhit[y][x])
                                thicks.append(thick_um[x])
                                break
                            if(x<(wid-1) and smoothhit[y][x+1]<=smoothhit[y][x] and smoothhit[y][x]<1.):
                                continue # ensure we start with increasing props
                        if(smoothhit[y][x]>maxprop): # stop adding points after first 1 is reached
                            props.append(smoothhit[y][x])
                            thicks.append(thick_um[x])
                            break
                        if(smoothhit[y][x]>minprop): # only add points with an interesting threshold
                            props.append(smoothhit[y][x])
                            thicks.append(thick_um[x])
                    propsarr.append(copy.deepcopy(props))
                    thicksarr.append(copy.deepcopy(thicks))
                    idxarr.append(y)
                    diams.append(diam_mm[y])
  
  
                if cs.verbose:
                    print("thicks")
                    print(thicksarr)
                    print("props")
                    print(propsarr)
                    print("done")
  
                # 4. try to sigmoid fit each curve
                for y in range(len(propsarr)):
                    props  = propsarr[y]
                    thicks = thicksarr[y]
                    if(props[0]>=thresh): # only point already more than enough
                        limval = thicks[0]
                    else:
                        if(len(thicks)>1):
                            guess = (0.,10.)
                            try:
                                popt, pcov = curve_fit(self._sigmoid, thicks, props, guess,diag=(1./np.mean(thicks), 1./np.mean(props)))
                                if type(pcov)==float and pcov == np.inf:
                                    print("ERROR sigmacurve fitting (pcov = inf)")
                                else:
                                    fitsarr.append(popt)
                                    fitidxs.append(y)
                                    limval = self._invsigmoid(thresh,*popt)
                            except:
                                print("ERROR sigmacurve fitting (no convergence, try larger guess)")
  
                        elif(len(thicks)==1):
                            if(props[0]>=thresh): # only point already more than enough
                                limval = thicks[0]
                            else:
                                limval = phantomlim[idxarr[y]][2] # only point not enough, should put it to max
                        else:
                            if(minprop>=thresh):
                                limval = phantomlim[idxarr[y]][1] # minprop set to high
                            else:
                                limval = phantomlim[idxarr[y]][2] # should never be the case
                    limit.append(limval)
                    human.append(phantom.humanfactor[idxarr[y]]*limval)
  
        cs.iqf = 0.
        cs.diam_mm = []
        cs.limit_um = []
        for d,t in zip(diams,limit):
            cs.iqf += d*t
            cs.diam_mm.append(d)
            cs.limit_um.append(t)
            
        if cs.verbose:
            print("diam_mm","limit_um")
            for d,t in zip(diams,limit):
                print(d,t)
            print("iqf",cs.iqf)
  
        plt.figure()
        plt.plot(diams,limit,label='automatic')
        plt.plot(diams,human,'c',label='pred. human')
        plt.ylabel("threshold thickness [um]")
        plt.xlabel("diameter [mm]")
        cs.hasmadeplots = True
        # diammm accept achiev
        protlim = [
            [2.00, 0.0690, 0.0380],
            [1.00, 0.0910, 0.0560],
            [0.50, 0.1500, 0.1030],
            [0.25, 0.3520, 0.2440],
            [0.10, 1.6800, 1.1000]
        ]
  
        plt.plot(np.array(protlim)[:,0],np.array(protlim)[:,1],'r',label='acceptable')
        plt.plot(np.array(protlim)[:,0],np.array(protlim)[:,2],'g',label='achievable')
        plt.plot(np.array(phantomlim)[:,0],np.array(phantomlim)[:,1],'--k',label='phantom limit')
        plt.plot(np.array(phantomlim)[:,0],np.array(phantomlim)[:,2],'--k',label='phantom limit')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        if not self.guimode:
            fname = 'threshold_'+str(cs.imnum).zfill(3)+'.jpg'
            plt.savefig(fname)
            cs.threshold_fnames.append(fname)
        
        if(len(fitsarr)>0):
            clist = ['c','m','y','r','g','b','k']
            sym = 'o'
            plt.figure()
            plt.xlabel("thickness [um]")
            plt.ylabel("probability")
            cx = 0
            for y in range(len(fitsarr)):
                idx = fitidxs[y]
                xarr = []
                yarr = []
                thsim = (np.linspace(0.01, 1., 50)).tolist()
                for x in range(wid):
                    if(gtmatrix[idx][x]==0):
                        continue
                    xarr.append(smoothhit[idx][x])
                    yarr.append(thick_um[x])

                vals = self._sigmoid(thsim,*(fitsarr[y]))
                plt.plot(yarr,xarr,sym,label=str(diam_mm[idx]),color=clist[cx])
                plt.plot(thsim,vals,'-',color=clist[cx])
                cx += 1
                if cx == len(clist): 
                    cx = 0
                    sym = 'd'
                
            plt.plot([0.01,10.],[thresh,thresh],'--k')
            plt.xscale('log')
            plt.legend()
            plt.ylim(0,1)

        if not self.guimode:
            fname = 'fit_'+str(cs.imnum).zfill(3)+'.jpg'
            plt.savefig(fname)
            cs.fit_fnames.append(fname)
  
    def pixDim(self,cs):
        try:
            pix_mm = cs.dcmInfile.PixelSpacing[0]
        except:
            pix_mm = cs.dcmInfile.ImagerPixelSpacing[0]
        return pix_mm

    def removeGrid(self,cs,removeGrid):
        """
        Detect the grid of the phantom and remove it from the image
        """
        shift = int(1./self.pixDim(cs)+.5)
        
        # try to find a threshold on pixelvalue to define a value representing the grid
        maskval = 0.75*cs.pixeldataIn.mean()
        # make a mask of grid-like values
        mask = (cs.pixeldataIn < maskval)

        # hole closing of the mask
        mask = scind.binary_closing(mask,structure=np.ones((5,5)))
        mask = scind.binary_opening(mask,structure=np.ones((5,5)))
        mask = scind.binary_dilation(mask)
        # new since 20150211
        mask = scind.binary_dilation(mask)

        # fill the gridlines with the median values of the pixels around it
        medimage = np.roll(cs.pixeldataIn,shift,axis=0).astype(float)
        dest = cs.pixeldataIn+mask*(medimage-cs.pixeldataIn)
        # repeat to remove propagated mask # new since 20150211
        medimage = np.roll(dest,shift,axis=0).astype(float)
        dest = cs.pixeldataIn+mask*(medimage-cs.pixeldataIn)
        medimage = None
        cs.gridimage = mask.astype(float)
        mask = None
        
        # find gridobject
        gridobject         = scind.binary_fill_holes(cs.gridimage)
        label_im,nb_labels = scind.label(gridobject)
        sizes = scind.sum(gridobject, label_im, range(nb_labels + 1))
        gridobject = None
        
        #Clean up small connect components:
        mask_size = sizes < max(sizes) #(100/self.pixDim())**2
        remove_pixel = mask_size[label_im]
        label_im[remove_pixel] = 0

        # Now reassign labels with np.searchsorted:
        labels = np.unique(label_im)
        label_im = np.searchsorted(labels, label_im)

        medimage = np.roll(dest,shift,axis=0).astype(float)
        dest += cs.gridimage*(medimage-dest)
        medimage = None
        cs.gridimage *= label_im
        if -1>0: # remove everything outside grid
            wid = dest.shape[0]
            hei = dest.shape[1]
            mv = np.mean(dest[int(wid/4):int(3*wid/4),int(hei/4):int(3*hei/4)])
            dest = label_im*(dest-mv)+mv

        if removeGrid:
            cs.pixeldataIn = dest

    def findStartingCell(self,cs):
        """
        Find first gridcell for gridcell locating.
        Required: 
          1. in image onscreen, phantom is located at lefthand side, with radboud logo in lower right corner, 
            with the irradiated area which is not covered by the phantom on the righthand side
          2. uncovered detector area has higher pixel value than covered area
          3. Phantom 3.2 cell is 1cm by 1 cm, approx 45 deg rotated. 
        Workflow:
        1. Find in gridlines image the first row with four crossings from non-grid to grid.
        2. top corner of first cell is mid of the last two crossings
        3. Other points from  cell is 1cm by 1 cm, approx 45 deg rotated
        4. Look for max in small area around each of the four cornerpoints
        """
        pixel_spacing  = self.pixDim(cs) # spacing in mm
        hei = cs.pixeldataIn.shape[1]

        roipts = []
        for y in range(hei):
            line = self._testline(cs.gridimage[:,y])
            if len(line) == 4:
                roipts.append( [int((line[2]+line[3])/2.),y] )
                break
        if len(roipts) == 0:
            print("Error! Could not find starting cell. Either grid not properly detected, or phantom rotation too large!")
            return roipts

        # half the length of the hypothenusa of a 90 deg triangle of sides 1 cm
        boxstep = int(10./pixel_spacing*np.sqrt(2.)/2)
        x0 = roipts[0][0]
        y0 = roipts[0][1]
        roipts.append([x0+boxstep,y0+boxstep])
        roipts.append([x0,y0+2*boxstep])
        roipts.append([x0-boxstep,y0+boxstep])
        
        roipts = self.AlignROI(cs, roipts, distmm=4)
        return roipts
    
    def CDMamSingle(self,cs):
        """
        Calculate score for one CDMam image
        WorkFlow:
        1. Extract grid from image (will be used to define regions of interest for observer)
        2. Find starting grid cell
        3. Locate all gridcells
        4. Test observer on each cell
        """
        self.plotcount = 0
        angledeg = 45.
        fudge = 0.97
        cos_a = fudge*np.cos(angledeg/180.*np.pi)
        sin_a = fudge*np.sin(angledeg/180.*np.pi)
        offsetmm_x = 94.
        offsetmm_y = 15.

        #1. remove the grid for will be used to define regions of interest for observer
        removeGridFromPixelData = True # for displaying only on gridless phantom
        print("1/4: removing grid....", end=" ")
        self.removeGrid(cs,removeGridFromPixelData)
        print("done")

        print("2/4: find starting grid cell....", end=" ")
        #2. First find some coordinate system om phantom
        pixel_spacing  = self.pixDim(cs) # spacing in mm
        cs.startingroi = self.findStartingCell(cs)
        print("done")
            
        #3. Locate all gridcells
        print("3/4: locating all grid cells....", end=" ")
        self.locateGridCells(cs)
        print("done")

        #4. Test observer on each cell
        print("4/4: modelobserver for each grid cell....", end=" ")
        score = self.observerScore(cs)
        print("done")

        return score
    
    def CDCOMSingle(self,cs):
        """
        run cdcom.exe on given file
        input matrix.inp and matrix2.inp, and average
        (matrix2.inp should be matrix.inp after plain_nnc; avg shouldbe matrix.inp after my ncc, but this is not the case)
        run ncc
        calc threshold thickness
        """
        # first clean up output of previous runs
        fnames = ["matrix.inp","matrix2.inp"]
        for g in fnames:
            try:
                os.remove(os.path.join(os.getcwd(),g))
            except:
                pass

        # make sure cdcom.exe is executable
        cdcomexe = os.path.join(os.path.dirname(__file__), "cdcom.exe")
        try: # make module executable
            os.chmod(cdcomexe, os.stat(cdcomexe).st_mode | stat.S_IEXEC)
        except Exception as e:
            print('cannot make cdcom.exe executable')
        
        # run cdcom on new file
        cmd = [cdcomexe, cs.imageFileName]
        cmd.extend(cs.parsCDCOM)
        print(cmd)
        with open(os.devnull, "w") as fnull:
            subprocess.check_call(cmd, stdout = fnull, stderr = fnull)

        # import cdcom output
        score = (np.zeros(np.shape(cs.phantom.groundtruth),dtype=float)).tolist()
        for fn in fnames:
            data = np.genfromtxt(os.path.join(os.getcwd(),fn),dtype=int,autostrip=True,delimiter="  ")
            dascore = data.tolist()
            # transform score:
            # undef: input -1 -> output 0
            # miss : input  2 -> output 0
            # hit  : input  1 -> output 1
            for y in range(len(dascore)):
                for x in range(len(dascore[y])):
                    if(dascore[y][x] == 1):
                        score[y][x] += 0.5

            try: # clean up
                os.remove(os.path.join(os.getcwd(),fn))
            except:
                pass


        if cs.verbose:
            for di in range(0,len(cs.phantom.diameter_mm)):
                for th in range(0,len(cs.phantom.thickness_um)):
                    if (cs.phantom.groundtruth[di][th]==0):
                        print("    ", end=" ")
                    else:
                        print("%0.2f" % score[di][th], end=" ")
                print("")
            print("")


        return score


    def observerScore(self,cs):
        """
        Calculate score of model observer:
        Workflow: 
          1. is max response in cell in the correct corner?
          2. consistency fix of scoring sheet
        """
        if cs.verbose:
            print('diam_mm','scale_mm2pix', 'diam_mm*scale_mm2pix', 'sigma_px','rad_px')

        #1. is max response in cell in the correct corner?
        GT = cs.phantom.groundtruth
        score = copy.deepcopy(GT)
        for di in range(0,len(cs.phantom.diameter_mm)):
            for th in range(0,len(cs.phantom.thickness_um)):
                if(GT[di][th] >0):
                    score[di][th] = self.detectVarGauss(cs,di,th)
        
        if cs.verbose:
            for di in range(0,len(cs.phantom.diameter_mm)):
                for th in range(0,len(cs.phantom.thickness_um)):
                    if (GT[di][th]==0):
                        print(" ", end=" ")
                    else:
                        print(score[di][th], end=" ")
                print("")
            print("")

        # 2. consistency fix of scoring sheet
        score = self.nearestNeighborCorrection(cs,score)
        #score = (np.zeros(np.shape(cs.phantom.groundtruth),dtype=float)).tolist()
        return score
    
    
    def locateGridCells(self,cs):
        """
        Starting from one a given box, locate all gridcells, using the previous one as a starting point.
        This way local deformations will not lead to mismatches further away/
        Returns a matrix or roicornerpoints in cs.gridrois
        """
        roipts_orig = cs.startingroi

        xstep_dx = (roipts_orig[1][0]-roipts_orig[0][0])
        xstep_dy = (roipts_orig[1][1]-roipts_orig[0][1])
        ystep_dx = -xstep_dy
        ystep_dy =  xstep_dx

        x0 = roipts_orig[0][0]  - 1*xstep_dx-0*ystep_dx
        y0 = roipts_orig[0][1]  - 1*xstep_dx-0*ystep_dy

        #print(x0,y0,xstep_dx,xstep_dy,ystep_dx,ystep_dy)
        """
        GT = 4 1
             3 2
        """
        GT = cs.phantom.groundtruth

        cs.gridrois = []
        for di in range(0,len(cs.phantom.diameter_mm)):
            lineroi = []
            for th in range(0,len(cs.phantom.thickness_um)):
                if(GT[di][th] >0):
                    if(th >0 and di >0 and lineroi[-1] != None and cs.gridrois[di-1][th] != None):
                        nwpts = []
                        nwpts.append([lineroi[-1][2][0]+xstep_dx,lineroi[-1][2][1]+xstep_dy])
                        nwpts = self.AlignROI(cs,nwpts,1.)
                        roipts = []
                        roipts.append(lineroi[-1][1])
                        roipts.append(cs.gridrois[di-1][th][2])
                        roipts.append(nwpts[0])
                        roipts.append(lineroi[-1][2])
                    elif(th>0 and lineroi[-1] != None):
                        nwpts = []
                        nwpts.append([lineroi[-1][1][0]+xstep_dx,lineroi[-1][1][1]+xstep_dy])
                        nwpts.append([lineroi[-1][2][0]+xstep_dx,lineroi[-1][2][1]+xstep_dy])
                        nwpts = self.AlignROI(cs,nwpts,1.)
                        roipts = []
                        roipts.append(lineroi[-1][1])
                        roipts.append(nwpts[0])
                        roipts.append(nwpts[1])
                        roipts.append(lineroi[-1][2])
                    elif(di>0 and cs.gridrois[di-1][th] != None):
                        nwpts = []
                        nwpts.append([cs.gridrois[di-1][th][2][0]+ystep_dx,cs.gridrois[di-1][th][2][1]+ystep_dy])
                        nwpts.append([cs.gridrois[di-1][th][3][0]+ystep_dx,cs.gridrois[di-1][th][3][1]+ystep_dy])
                        nwpts = self.AlignROI(cs,nwpts,1.)
                        roipts = []
                        roipts.append(cs.gridrois[di-1][th][3])
                        roipts.append(cs.gridrois[di-1][th][2])
                        roipts.append(nwpts[0])
                        roipts.append(nwpts[1])
                    else:
                        roipts = []
                        roipts.append([int(x0+(th+0)*xstep_dx+(di+0)*ystep_dx+.5),int(y0+(th+0)*xstep_dy+(di+0)*ystep_dy)])
                        roipts.append([int(x0+(th+1)*xstep_dx+(di+0)*ystep_dx+.5),int(y0+(th+1)*xstep_dy+(di+0)*ystep_dy)])
                        roipts.append([int(x0+(th+1)*xstep_dx+(di+1)*ystep_dx+.5),int(y0+(th+1)*xstep_dy+(di+1)*ystep_dy)])
                        roipts.append([int(x0+(th+0)*xstep_dx+(di+1)*ystep_dx+.5),int(y0+(th+0)*xstep_dy+(di+1)*ystep_dy)])
                        roipts = self.AlignROI(cs,roipts,2)

                    lineroi.append(copy.deepcopy(roipts))
                else:
                    lineroi.append(None)
            cs.gridrois.append(copy.deepcopy(lineroi))


    def AlignROI(self, cs,roipts,distmm=5):
        """
        Moves the four cornerpoints of a ROI to the most likely crossing in the gridimage, within <distmm> distance from the starting point.
        Returns four adjusted cornerpoints.
        """
        pix2phantommm = self.pixDim(cs)
        searchrad = int(distmm/pix2phantommm+.5)

        searchrad = max(1,searchrad)
        widthpx = np.shape(cs.pixeldataIn)[0] ## width/height in pixels
        heightpx = np.shape(cs.pixeldataIn)[1]

        conf_pts = []
        conf_pts.append((0,copy.deepcopy(roipts)))
        for kk in range(0,3): #6
            for i in range(0, len(roipts)):
                rp = roipts[i]
                x0 = rp[0]
                minx = max(0,x0-searchrad)
                maxx = min(widthpx-2,x0+searchrad)
                y0 = rp[1]
                miny = max(0,y0-searchrad)
                maxy = min(heightpx-2,y0+searchrad)
                cropped = cs.gridimage[minx:maxx+1,miny:maxy+1]
                # smoothing to get rig of noise
                sigma = max(1.5,cropped.shape[0]/8.)
                cropped = scind.gaussian_filter(cropped.astype(float), sigma,mode='constant')
                (x1,y1) = np.unravel_index(cropped.argmax(),cropped.shape)
                x1 += minx
                y1 += miny
                rp[0] = x1
                rp[1] = y1
            searchrad = int(max(1,searchrad/2))
        return roipts

    def nearestNeighborCorrection(self,cs,scorematrix):
        """
        Nearest Neighbor Correction.
        From the CDMAM manual:
            1. A True needs 2 or more correctly indicated nearest neighbours to remain a True
            2. A False or Not indicated disk will be considered as True when it has 3 or 4 correctly indicated nearest neighbours.
            3. A True which has only 2 nearest neighbours (at the edges of the phantom) needs only 1 correctly indicated nearest neighbour to remain True.
            4. A False or Not indicated disk which has only 2 nearest neighbours will be regarded True if both nearest neighbours are correctly indicated.
            5. The absent corners of the phantom (0.03 mm/2.0 mm and 2.00mm/0.06 mm) will be regarded True when both nearest neighbours are correctly indicated.
        """
        gtmatrix = cs.phantom.groundtruth
        wid = len(gtmatrix[0])
        hei = len(gtmatrix)
        #  0 : not defined or no hit
        #  1 : hit
        #  0.5 : flipped
        hitmatrix = copy.deepcopy(scorematrix)
        nn = [[-1,0],[0,1],[1,0],[0,-1]]
        for y in range(hei):
            for x in range(wid):
                if(gtmatrix[y][x] == 0):
                    continue
                trues = 0
                count = 0
                for n in nn:
                    x1 = x+n[0]
                    y1 = y+n[1]
                    if(x1<0 or y1<0 or x1>(wid-1) or y1>(hei-1)):
                        continue
                    if(gtmatrix[y1][x1] == 0):
                        continue
                    count += 1
                    if(scorematrix[y1][x1] == 1):
                        trues += 1

                if(scorematrix[y][x] == 1): # True
                    if (count>2 and trues<2) or (count==2 and trues<1): # rule1: True needs at least 2 nn Trues
                        hitmatrix[y][x] = 0.5
                elif(scorematrix[y][x] == 0): # False
                    if (count>2 and trues>2) or (count==2 and trues>1): # rule2: False needs at least 3 nn Trues
                        hitmatrix[y][x] = 0.5

        # rule 3:
        if(scorematrix[0][1] == 1 and scorematrix[1][0] == 1):
            hitmatrix[0][0] = 0.5
        if(scorematrix[hei-2][wid-1] == 1 and scorematrix[hei-1][wid-2] == 1):
            hitmatrix[hei-1][wid-1] = 0.5

        if cs.verbose:
            for y in range(hei):
                for x in range(wid):
                    if (gtmatrix[y][x]==0):
                        print("    ", end=" ")
                    else:
                        print("%0.2f" % hitmatrix[y][x], end=" ")
                print("")
            print("")
        return hitmatrix


    def detectVarGauss(self,cs, di,th):
        """
        Blur with Gaussian of sigma matched to diameter of disc == averaging over disc size
        Detect max response within search radius around expected location
        """
        mode = 0 # Gaussian best sofar
        # mode = 1 # Trace Lxx+Lyy, look for max next best
        # mode = 2 # Det Lxx*Lyy-Lxy*Lyx, look for max

        bbox = cs.gridrois[di][th]
        gt   = cs.phantom.groundtruth[di][th]
        diam_mm = cs.phantom.diameter_mm[di]
        thick_um = cs.phantom.thickness_um[th]

        scale_mm2pix = 0
        for i in range(len(bbox)):
            scale_mm2pix += np.sqrt((bbox[i][0]-bbox[i-1][0])**2+(bbox[i][1]-bbox[i-1][1])**2)
        scale_mm2pix *= 0.25/10.
        rad_px = int(diam_mm/2.*scale_mm2pix+.5)

        minx = bbox[0][0]
        maxx = bbox[0][0]
        miny = bbox[0][1]
        maxy = bbox[0][1]
        for b in bbox:
            minx = min(minx,b[0])
            maxx = max(maxx,b[0])
            miny = min(miny,b[1])
            maxy = max(maxy,b[1])
        minx -= rad_px
        maxx += rad_px
        miny -= rad_px
        maxy += rad_px
        crop = cs.pixeldataIn[minx:maxx,miny:maxy].astype(float)

        if(mode == 0): # integrate /avg over area
            sigma_px = max(1.5,rad_px/2.) #max(1.5,rad_px/6.)
            crop = scind.gaussian_filter(crop,sigma_px,order=[0,0]) #,mode='constant')
        elif(mode == 1):
            sigma_px = max(3.5,rad_px/2.)
            Lxx = scind.gaussian_filter(crop,sigma_px,order=[2,0]) #,mode='constant')
            Lyy = scind.gaussian_filter(crop,sigma_px,order=[0,2])
            crop = Lxx+Lyy
        elif(mode == 2):
            sigma_px = max(3.5,rad_px/2.)
            Lxx = scind.gaussian_filter(crop,sigma_px,order=[2,0]) #,mode='constant')
            Lyy = scind.gaussian_filter(crop,sigma_px,order=[0,2])
            Lxy = scind.gaussian_filter(crop,sigma_px,order=[1,1])
            crop = Lxx*Lyy-Lxy*Lxy

        offset_mm = 2.85/np.sqrt(2.)
        xstep_dx = int(offset_mm/10.*(bbox[1][0]-bbox[0][0])+.5)
        xstep_dy = int(offset_mm/10.*(bbox[1][1]-bbox[0][1])+.5)
        ystep_dx = -xstep_dy
        ystep_dy =  xstep_dx
        # search radius
        rad_px = int(.25*scale_mm2pix+.5)#/2 #int(sigma_px+.5) #int(3./4.*scale_mm2pix+.5) # max diam is 3.2 mm
        rad2_px = rad_px**2
        pos = [
            [bbox[0][0]-minx+xstep_dx+ystep_dx,bbox[0][1]-miny+xstep_dy+ystep_dy], # GT =1
            [bbox[1][0]-minx-xstep_dx+ystep_dx,bbox[1][1]-miny-xstep_dy+ystep_dy], # GT =2
            [bbox[2][0]-minx-xstep_dx-ystep_dx,bbox[2][1]-miny-xstep_dy-ystep_dy], # GT =3
            [bbox[3][0]-minx+xstep_dx-ystep_dx,bbox[3][1]-miny+xstep_dy-ystep_dy], # GT =4
        ]
        if cs.verbose:
            print(diam_mm,scale_mm2pix, diam_mm*scale_mm2pix, sigma_px,rad_px)
        r0 = 0.
        if(mode == 0):
            r0= 1.e6
        else:
            r0= -1.e6# for trace or det
        response = [r0,r0,r0,r0]
        """
        GT = 4 1
             3 2

        """
        score = [4,1,2,3]
        for p,i in zip(pos,range(4)):
            for y in range(-rad_px,rad_px):
                r2my2 = rad2_px-y**2
                for x in range(-rad_px,rad_px):
                    if x**2 < r2my2:
                        if(mode == 0):
                            response[i] = min(response[i],crop[p[0]+x,p[1]+y])
                        else:
                            response[i] = max(response[i],crop[p[0]+x,p[1]+y]) # for Trace or Det

        if(mode == 0):
            rindex = np.array(response).argmin()
        else:
            rindex = np.array(response).argmax()
        if np.abs(response[rindex]-r0) <1.e-6: # check if we have an opinion
            guess = 0
        else:
            guess = score[rindex]

        if(1<0 and guess != gt and cs.plotcount<15 and self.guimode):
            print('guess',guess,gt,end=" ")
            for s,r in zip(score,response):
                print(s,':',r,end=" ")
            print
            plt.figure()
            if(mode == 0):
                plt.imshow(crop.transpose()).set_clim(min(response),max(response)) # (2600.,2800.)
            elif(mode == 1):
                plt.imshow(crop.transpose()).set_clim(-15.,15.)
            elif(mode == 2):
                plt.imshow(crop.transpose()).set_clim(-200.,200.)

            # black plussed for locations of centers of to look-for locations
            for pp in pos:
                plt.plot(pp[0],pp[1],'k+')
                circle2=plt.Circle((pp[0],pp[1]),rad_px,color='k',fill=False)
                plt.gcf().gca().add_artist(circle2)

            # white plus at location of truth
            pix = gt
            if gt == 4:
                pix = 0
            plt.plot(pos[pix][0],pos[pix][1],'w+')
            circle2=plt.Circle((pos[pix][0],pos[pix][1]),rad_px,color='w',fill=False)
            plt.gcf().gca().add_artist(circle2)


            # pink plus at location of max response
            nccloc = np.nonzero(crop == response[rindex])
            plt.plot(nccloc[0],nccloc[1],'m+')
            circle2=plt.Circle((nccloc[0],nccloc[1]),rad_px,color='m',fill=False)
            plt.gcf().gca().add_artist(circle2)
 
            # pink plus at location of max response
            plt.title(str(diam_mm)+" "+str(thick_um))
            plt.colorbar()
            cs.plotcount += 1
            cs.hasmadeplots = True

        # transform into hit or no hit
        if guess == gt:
            guess = 1
        else:
            guess = 0
        return guess

    def determineScannerID(self,cs):
        """
        Tries to identify machine model, used filter, and energy presentation
        """
        dicomfields = [
            ["0008,1090",  "ModelName"],
            ["0018,7050",  "FilterMaterialLT"],
            ["0019,10C1",  "MICRODOSE IMAGE CONTENT"] #---: SUM FOR PRESENTATION
        ]


        # 1. Try to id Scanner
        cs.scannername = lit.stUnknown
        dicomvalue = wadwrapper_lib.readDICOMtag(dicomfields[0][0],cs.dcmInfile) # Manufacturer's Model Name
        dicomvalue = str(dicomvalue).lower()
        if dicomvalue.find("l50")>-1:
            cs.scannername = lit.stL50
        elif dicomvalue.find("lorad selenia")>-1:
            cs.guessScanner = lit.stSelenia
        elif dicomvalue.find("dimensions")>-1:
            cs.guessScanner = lit.stDimensions

        # 2. Try to id Filter
        cs.filtername = wadwrapper_lib.readDICOMtag(dicomfields[1][0],cs.dcmInfile) # Filtername

        # 3. Try to id EnergyPresentation
        cs.energypresentation = wadwrapper_lib.readDICOMtag(dicomfields[2][0],cs.dcmInfile)


    def DICOMInfo(self,cs,info='dicom'):
        """
        Extract some info from the DICOM header. <info> is either 'dicom' (extensive list) or 'qc' (short list).
        Returns list of <tag description>,<tag value>
        """
        # Different from ImageJ version; tags "0008","0104" and "0054","0220"
        #  appear to be part of sequences. This gives problems (cannot be found
        #  or returning whole sequence blocks)
        # Possibly this can be solved by using if(type(value) == type(dicom.sequence.Sequence()))
        #  but I don't see the relevance of these tags anymore, so set them to NO

        self.determineScannerID(cs)
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

        results = []
        for df in dicomfields:
            key = df[0]
            value = ""
            try:
                value = wadwrapper_lib.readDICOMtag(key,cs.dcmInfile)
            except:
                value = ""

            results.append( (df[1],value) )

        return results

