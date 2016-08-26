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
Analysis of reverberations in Air pattern.
Workflow:
input: dicomobject,pixeldata  (the pixeldata is already cleaned; color is removed, leaving plain grayscale)
1. Isolate the reverberations pattern from the other information:
  a. Connected components analysis of pixelvalues >0
  b. Merge all not-to small components found along the y-axis in the middle of the imagedata, 
     excluding components fixed to the top of the image
  c. Use as a mask to isolate reverberations pattern
2. Straighten curve probe data:
  a. Detect points left and right from the middle, along the top of the reveberation pattern; if too few points
     are detected, the probe is not curved, so just crop the reveberation pattern -> rect_image
  b. Fit a circle through these points
  c. Interpolate data to straightend grid -> rect_image
3. Sensitivity analyis of rect_image to determine vertical extend of reverberation pattern:
  a. Make a profile along the vertical: average of cetral part of rect_image
  b. Extend: Determine noise_level, and calculate intersection of profile with noise_level
  c. Extend_fft: Do a fourier analysis to find the ringing pattern; remove it an find the minimum of the 
     background trend
4. Uniformity analysis of rect_image to find dips that suggest problems with the probe:
  a. Make a profile along the horizontal of the rect_image upto sensitivity_extend;
     FASTMODE: (experimental) use onle rect_image of first ring; should have all the info there already
  b. Detect dips that are deep enough, find the minimum in of the dips, and overal and the overal non-uniformity
  

Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

TODO:
    o Maybe put hard threshold on peak value for uniformity (is normalized, so why not?!)
Changelog:
    20160802: python3 compatible; fix ints
    20150626: remove division by zero for snr
    20150610: moved peakfit to wadwrapper_lib
    20150520: cleanup label
    20150429: Fixes reading of files wihtout NumberOfFrames tag
    20150422: Clean up; added documentation
    20150421: Added fft analysis; clean-up; added FAST_MODE: scipy.map_coordinates, x-prof only first ring; localnormalization
    20150420: Overal non-uniformity added; allow min/max at start/end
    20150417: Removed valleyfit; fixed radius range
    20150416: Sensitivity analysis and uniformity analysis
    20150410: Initial version
"""
__version__ = '20160802'
__author__ = 'aschilham'

import copy
try:
    from pyWADLib import wadwrapper_lib
except ImportError:
    import wadwrapper_lib
import numpy as np
import matplotlib.pyplot as plt
import operator
import scipy
import time
import os
import QCUS_math as mymath

try:
    import numpy.fft.rfftfreq as rfftfreq
except: # old numpy
    def rfftfreq(n, d=1.0):
        """
        Return the Discrete Fourier Transform sample frequencies
        (for usage with rfft, irfft).
        The returned float array `f` contains the frequency bin centers in cycles
        per unit of the sample spacing (with zero at the start).  For instance, if
        the sample spacing is in seconds, then the frequency unit is cycles/second.
        Given a window length `n` and a sample spacing `d`::
          f = [0, 1, ...,     n/2-1,     n/2] / (d*n)   if n is even
          f = [0, 1, ..., (n-1)/2-1, (n-1)/2] / (d*n)   if n is odd
        Unlike `fftfreq` (but like `scipy.fftpack.rfftfreq`)
        the Nyquist frequency component is considered to be positive.
        Parameters
        ----------
        n : int
            Window length.
        d : scalar, optional
            Sample spacing (inverse of the sampling rate). Defaults to 1.
        Returns
        -------
        f : ndarray
            Array of length ``n//2 + 1`` containing the sample frequencies.
        """

        val = 1.0/(n*d)
        N = n//2 + 1
        results = np.arange(0, N, dtype=int)
        return results * val    

def cropImage(image,offX,offY):
    # helper function for symmetric cropping of image
    wid,hei = np.shape(image)
    if len(offX)==1:
        xran = [offX[0],-offX[0]]
    else:
        xran = offX
    if xran[1] == 0:
        xran[1] = wid-1
        
    if len(offY)==1:
        yran = [offY[0],-offY[0]]
    else:
        yran = offY
    if yran[1] == 0:
        yran[1] = hei-1
    return image[xran[0]:xran[1],yran[0]:yran[1]]

class USStruct:
    verbose = False
    testing = False

    # input image
    dcmInfile = None
    pixeldataIn = None
    dicomMode = wadwrapper_lib.stMode2D

    
    # for matlib plotting
    hasmadeplots = False

    # reverbrations
    reverb_image  = None # isolated reverbrations (by largest connected component)
    rect_image    = None # straightend image for curvilinear probes
    curve_xyr     = None # center and radius [x,y,r] of curvilinear probe
    curve_residu  = None # residu of fitting of curvilinear probe
    curve_opendeg = None # opening angle in degrees
    unif_bots     = None # list of (pos,val) of kept valleys in normalized reverb profile
    unif_lowest   = 0 # lowest value of valleys
    unif_low      = 0 # lowest value 
    unif_line     = 0 # Overal non-uniformity over whole, unnormalized pofile
    
    rev_miny = None # (x,y) of position of minimum y of reverb image
    rev_maxy = None # (x,y) of position of maximum y of reverb image
    rev_maxx = None # (x,y) of position of maximum x of reverb image
    rev_minx = None # (x,y) of position of minimum x of reverb image
    rev_midx = None # estimated x of middle of reverb image
    rev_mask = None # mask of reverb image
    
    curve_radii  = None # [min, max] radius
    curve_angles = None # [min, max] angle

    # sensitivity
    sens_ylim   = None    # limit on y determined from peaks in y-dir 
    sens_ylim2  = None   # limit on y determined from fft analysis of peaks in y-dir 
    sens_noiseM = None # noise determined from (vertical) sensitivity profile
    sens_numtops = None # number of peaks found upto ylim
    sens_numbots = None # number of valleys found upto ylim
    sens_noiserange = None # number of lines used for calculation of noiselevel
    sens_basewavelength = None # basic wavelength i.e. distance between rings in mm
    
    # images
    image_fnames = [] # filenames of generated images
    
    # for gui
    cca_image = None # connected component analysis
    report = None # list ('description',value) of timing information
    
    def __init__ (self,dcmInfile,pixeldataIn,dicomMode):
        self.verbose = False
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.dicomMode = dicomMode
        self.hasmadeplots = False
        self.report = []
        
        self.cca_image = None
        self.reverb_image = None
        self.rect_image = None
        self.curve_xyr = [0,0,0]
        self.curve_residu = 0
        self.curve_opendeg = 0
        self.unif_bots = []
        self.unif_lowest = 0
        self.unif_low = 0
        self.unif_line = 1

        self.sens_ylim = 0
        self.sens_ylim2 = 0
        self.sens_noiseM = 0
        self.sens_numtops = 0
        self.sens_numbots = 0
        self.sens_noiserange = 2 
        self.sens_basewavelength = 0.
        
        # helper stuff
        rev_miny = None 
        rev_maxy = None 
        rev_maxx = None 
        rev_minx = None 
        rev_midx = None 
        curve_angles = []
        curve_radii  = []
        rev_mask = None
        
class US_QC:
    qcversion = __version__
    fastmode = False
    modeLocalNorm = False

    # parameters:
    smooth_uniformity = 1  # running average of this width before peak detection
    delta_uniformity = .15 # A dip in normalized reverb pattern must be at least <delta>
    smooth_sensitivity = 3 # running average of this width for sensitivity data
    fdelta_sensitivity = .10 # A peak in sensitivity profile must be at least <fdelta>*(max-noise)
    offsetVER = 0 # default lines to exclude from top and bottom when making profiles (10)
    offsetHOR = 0 # default rows to exclude from top and bottom when making profiles (10)
    fft_skip = 5 # in fft analysis, ignore peak if it is located within skip freqs from 0
                 # to ignore for peak finding (offset, trend)
    check_asym = False # do or do not check of strongly asym peaks
    sens_hicut = True # restrict peaks to max > noise (True) or min < noise

    def __init__(self,guimode=False,fastmode=False,modelocalnorm=False):
        self.guimode = guimode
        self.fastmode = fastmode
        self.modeLocalNorm = modelocalnorm # EXPERIMENTAL
        
    def readDICOMtag(self,cs,key,imslice=0): # slice=2 is image 3
        value = wadwrapper_lib.readDICOMtag(key,cs.dcmInfile,imslice)
        return value

    def DICOMInfo(self,cs,info='dicom'):
        # Different from ImageJ version; tags "0008","0104" and "0054","0220"
        #  appear to be part of sequences. This gives problems (cannot be found
        #  or returning whole sequence blocks)
        # Possibly this can be solved by using if(type(value) == type(dicom.sequence.Sequence()))
        #  but I don't see the relevance of these tags anymore, so set them to NO

        if info == "dicom":
            dicomfields = [
                ["0008,0012", "Instance Date"],
                ["0008,0013", "Instance Time"],
                ["0008,0060", "Modality"],
                ["0008,0070", "Manufacturer"],
                ["0008,1090", "Manufacturer Model Name"],
                ["0008,1010", "Station Name"],
                ["0008,1030", "Study Description"],
                ["0008,0068", "Presentation Intent Type"], 
                ["0018,1000", "Device Serial Number"],
                ["0018,1020", "Software Version(s)"],
                ["0018,1030", "Protocol Name"],
                ["0018,5010", "Transducer Data"],
                ["0018,5020", "Processing Function"],
                ["0028,0002", "Samples per Pixel"],
                ["0028,0101", "Bits Stored"],
                ["0028,2110", "Lossy Image Compression"],
                ["2050,0020", "Presentation LUT Shape"],
                ["0018,6011, 0018,6024", "Physical Units X Direction"],
                ["0018,6011, 0018,602c", "Physical Delta X"],
            ] # Philips

        elif info == "id":
            dicomfields = [
                ["0008,1010", "Station Name"],
                ["0018,5010", "Transducer"],
                ["0008,0012", "InstanceDate"],
                ["0008,0013", "InstanceTime"]
            ]
        elif info == "probe":
            dicomfields = [
                ["0018,5010", "Transducer"],
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
    def isolateReverbrations(self,cs,thresh=0):
        """
        Find reverbrations part of image.
        Workflow:
        1. Find reverberations as largest connected component != 0
        2. Return
        """
        error = True
        # cluster connected components with pixelvalues>0
        time0 = time.time()
        #thresh = 0
        #work = (cs.pixeldataIn>0) * (cs.pixeldataIn<255) # must give a signal, but 255 often reserved for overlay
        work = cs.pixeldataIn>thresh
        cca = wadwrapper_lib.connectedComponents()
        cs.cca_image,nb_labels = cca.run(work)
        mode = 'all_middle' #'largest_only'#'all_middle'
        if mode == 'largest_only': # model of Pepijn
            # select only largest cluster
            cluster_sizes = cca.clusterSizes()
            clus_val = np.argmax(cluster_sizes)
            cs.rev_mask = (cs.cca_image == clus_val)
            cs.report.append(('cca',time.time()-time0))
            clus = cca.indicesOfCluster(clus_val)
        else: 
            # select all clusters present in vertical middle area of image, excluding top and 0
            wid,hei = np.shape(cs.cca_image)
            # first remove very small clusters (can be located around top edge!)
            cca.removeSmallClusters(wid/10*hei/10)
            search = cs.cca_image[int(0.4*wid):int(0.6*wid),:]
            labs = []
            for ss in search.ravel():
                if ss>0 and not ss in labs:
                    labs.append(ss)
            # exclude labels in top rows (never part of imagedata, but full of annotations)
            search = cs.cca_image[:,0:5]
            notlabs = []
            for ss in search.ravel():
                if ss>0 and not ss in notlabs:
                    notlabs.append(ss)
            labs = [la for la in labs if la not in notlabs]

            cs.rev_mask = np.reshape(np.in1d(cs.cca_image,labs),np.shape(cs.cca_image))
            clus = np.where(cs.rev_mask)
            clus = [ (x,y) for x,y in zip(clus[0],clus[1]) ]
            
        # make an image of only largest cluster applied to original pixeldata
        cs.reverb_image = cs.pixeldataIn*cs.rev_mask
        cs.report.append(('mask',time.time()-time0))

        # bounds of reverb image
        time0 = time.time()
        cs.rev_miny = min(clus,key=operator.itemgetter(1))
        cs.rev_maxy = max(clus,key=operator.itemgetter(1))
        cs.rev_maxx = max(clus,key=operator.itemgetter(0))
        cs.rev_minx = min(clus,key=operator.itemgetter(0))
        cs.rev_midx = int(np.shape(cs.rev_mask)[0]/2)
        
        error = False
        return error

    def straightenCurve(self,cs):
        """
        Straighten curved reverb image if applicable
        Workflow:
        1. Fit circle through top of reverb data
        """
        error = True
        
        ## 2. Transform reverb data to rectangle if needed
        time0 = time.time()
        ## Fit a circle to top of reverb pattern
        # From top of reverb image down, look for pixels != 0 from mid to left and from mid to right;
        # if both found, add it to the list
        circL_xy = []
        circR_xy = []
        for y in range(cs.rev_miny[1],cs.rev_maxy[1]):
            for x in reversed(range(cs.rev_minx[0],cs.rev_midx)): #find left point
                xl = -1
                xr = -1
                if cs.rev_mask[x,y]:
                    xl = x
                    break
            for x in range(cs.rev_midx,cs.rev_maxx[0]): #find right point
                if cs.rev_mask[x,y]:
                    xr = x
                    break
            if xl>-1 and xr>-1:
                circL_xy.append((xl,y))
                circR_xy.append((xr,y))
                if xr-xl<10: # stop looking
                    break
        circ_xy = []
        circ_xy.extend(circL_xy)
        circ_xy.extend(circR_xy)
        circ_xy.sort(key=operator.itemgetter(1))
        if len(circ_xy)>10: # at least 10 point, else probably not a curved probe
            # use only central part for fitting, as deformations towards edges occur
            cf = mymath.CircleFit(circ_xy[int(2.*len(circ_xy)/3.):]) 
            fit = cf.fit()
            (xc,yc) = fit[0]
            Rc = fit[1]
            
            # calculate limiting angles and radii
            cs.curve_angles = [np.arctan2(circL_xy[0][0]-xc,circL_xy[0][1]-yc),np.arctan2(circR_xy[0][0]-xc,circR_xy[0][1]-yc)]
            maxrad = min( [
                (cs.rev_minx[0]-xc)/np.sin(cs.curve_angles[0]),
                (cs.rev_maxx[0]-xc)/np.sin(cs.curve_angles[1]),
                (cs.rev_maxy[1]-yc)
            ])
            cs.curve_radii = [Rc,maxrad ]
            cs.curve_xyr = [xc,yc,Rc]
            cs.curve_residu = cf.residu
            cs.curve_opendeg = (cs.curve_angles[1]-cs.curve_angles[0])/np.pi*180.
    
            cs.report.append(('fit',time.time()-time0))
            
            # transform reverb pattern to rectangle: interpolate at coords
            time0 = time.time()
            ang = np.linspace(cs.curve_angles[0],cs.curve_angles[1],cs.rev_maxx[0]-cs.rev_minx[0])
            rad = np.linspace(cs.curve_radii[0],cs.curve_radii[1],cs.rev_maxy[1]-cs.rev_miny[1])
            an,ra = scipy.meshgrid(ang,rad)
            xi = cs.curve_xyr[0]+ra*np.sin(an)
            yi = cs.curve_xyr[1]+ra*np.cos(an)
            coords = np.array( [xi,yi] )
            cs.rect_image = scipy.ndimage.map_coordinates(cs.reverb_image, coords).transpose()
            cs.report.append(('straighten',time.time()-time0))
        else: # not a curved probe
            cs.rect_image = cs.reverb_image[cs.rev_minx[0]:cs.rev_maxx[0],cs.rev_miny[1]:cs.rev_maxy[1]]
            
        error = False
        return error

    def _uniformityAnalysis(self,cs,normuniformity):
        # Helper function for analysis of normalized uniformity profile:
        error = True
        nx = len(normuniformity)
        time0 = time.time()
        self.fftAnalysis(cs,normuniformity,mode='uniformity')
        cs.report.append(('unif_fft',time.time()-time0))

        # peak detection with ugly hack to also allow peaks at index 0 and last
        #xy_max,xy_min = wadwrapper_lib.peakdet(normuniformity, delta=self.delta_uniformity,x=range(nx))
        p0 = normuniformity[::-1]
        p2 = np.append(np.append(p0[:-1],normuniformity),p0[1:]) 
        x2 = np.array(range(-(nx-1),2*nx-1))
        xy_max,xy_min = wadwrapper_lib.peakdet(p2, delta=self.delta_uniformity,x=x2)
        xy_max = [xy for xy in xy_max if xy[0]>-1 and xy[0]<nx]
        xy_min = [xy for xy in xy_min if xy[0]>-1 and xy[0]<nx]

        ## peak analysis
        # sorted list of extrema positions:
        minmax_pos = [da[0] for da in xy_min]
        minmax_pos.extend( [da[0] for da in xy_max] )
        minmax_pos = sorted(minmax_pos)

        # generate an image of non-uniformity
        x_max = np.array([xy[0] for xy in xy_max])
        y_max = np.array([xy[1] for xy in xy_max])
        x_min = np.array([xy[0] for xy in xy_min])
        y_min = np.array([xy[1] for xy in xy_min])

        yran = []
        if self.fastmode: # use only first ring
            print('Uniformity: FASTMODE',cs.sens_basewavelength)
            if cs.sens_basewavelength >0:
                ylim = int(cs.sens_basewavelength/self.pixels2mm(cs,1.)+.5)
                yran = [self.offsetVER,ylim]

        if len(yran) == 0:
            yran = [self.offsetVER,cs.sens_ylim] if cs.sens_ylim>0 else [self.offsetVER]

        plt.figure()
        crop = cropImage(cs.rect_image, [self.offsetHOR], yran)
        plt.imshow(crop.transpose(),cmap=plt.gray())
        x = np.array(range(nx))
        plt.plot(x,-100*normuniformity,'y-',label='100*normalized')
        plt.plot(x_min,-100*y_min,'ro',label='accepted dips')
        plt.plot(x_max,-100*y_max,'bo',label='accepted peaks')
        plt.xlim([0,nx])
        if not self.guimode:
            fname = 'uniformity_'+self.imageID(cs,probeonly=True)+'.jpg'
            plt.savefig(fname)
            cs.image_fnames.append(fname)
        else:
            cs.hasmadeplots = True

        if cs.verbose:
            plt.figure()
            plt.plot(x,normuniformity,'k-')
            plt.plot(x_min,y_min,'ro')
            plt.plot(x_max,y_max,'bo')
            cs.hasmadeplots = True
            
        # possible clean-up of peaks detected
        cs.unif_bots = []
        for damin in xy_min:
            # find range around each dip and remove strongly asymetric dips
            if self.check_asym:
                peak_pos = damin[0]
                peak_val = damin[1]
                asym = True
                delta = self.delta_uniformity
                while asym:
                    pos_left = peak_pos
                    for ix in reversed(range(0,peak_pos)):
                        if normuniformity[ix]>peak_val+delta:
                            pos_left= ix
                            break
                    pos_right = peak_pos
                    for ix in range(peak_pos,nx):
                        if normuniformity[ix]>peak_val+delta:
                            pos_right= ix
                            break
                    if pos_left == peak_pos:
                        print('ERROR! Could not find left point of peak',peak_ix,peak_pos,peak_val)
                        return error
                    if pos_right == peak_pos:
                        print('ERROR! Could not find left point of peak',peak_ix,peak_pos,peak_val)
                        return error
                    # check if we are looking at a strongly symmetric peak
                    if (pos_right-pos_left)>4*min(peak_pos-pos_left,pos_right-peak_pos):
                        delta -= .1*self.delta_uniformity
                        if delta < .2*self.delta_uniformity:
                            break
                    else:
                        asym = False
                if asym:
                    print('ERROR! Very asymmetric peak',peak_ix,peak_pos,peak_val,peak_pos-pos_left,pos_right-peak_pos)
                    continue

            # acceptable
            cs.unif_bots.append(damin)

        # some figures of merit to store: non-unif in deepest dip and overal
        cs.unif_lowest = 1. if len(cs.unif_bots)==0 else min(cs.unif_bots,key=operator.itemgetter(1))[1]
        cs.unif_low = min(normuniformity[1:-1]) # exclude bounds
        error = False
        return error

    def reverbUniformity(self,cs):
        """
        Analysis uniformity of reverb pattern
        Workflow:
        1. Define ROI of reverb (dependend on vertical profile)
        2. Average ROI horizontally and calculate COV
        3. Normalize profile as deviation wrt mean
        4. Find dips in profile
        """
        error = True
        
        ## 1. Define ROI of reverb (dependend on vertical profile)
        yran = []
        if self.fastmode:
            print('Uniformity: FASTMODE',cs.sens_basewavelength)
            if cs.sens_basewavelength >0:
                ylim = int(cs.sens_basewavelength/self.pixels2mm(cs,1.)+.5)
                yran = [self.offsetVER,ylim]
            print(yran)
        if len(yran) == 0:
            yran = [self.offsetVER,cs.sens_ylim] if cs.sens_ylim>0 else [self.offsetVER]

        crop = cropImage(cs.rect_image, [self.offsetHOR], yran)

        uniformity = np.mean(crop,axis=1)
        uniformitysmoothed = mymath.movingaverage(uniformity,self.smooth_uniformity)
        
        ## line uniformity
        intemin = np.min(uniformity)
        intemax = np.max(uniformity)
        inteavg = np.mean(uniformity)
        cs.unif_line = np.max([intemax-inteavg,inteavg-intemin])/inteavg

        # Output of COV 
        meanvalue = np.mean(uniformity) 
        stdvalue = np.std(uniformity)

        # normalized uniformity 
        if self.modeLocalNorm:
            frac = .1
            print("Uniformity, using Local Normalization",frac)
            normuniformity = mymath.localnormalization(uniformity, sigma=len(uniformity)*frac)
        else:
            normuniformity = ((uniformitysmoothed-meanvalue)/meanvalue)
            
        if cs.testing:
            try:
                import QCUS_testing as mytest
                mytest._exportProfile(normuniformity,fname='xprofile.tsv')
                mytest._exportNDArray(cs.rect_image)
            except:
                pass
            
        return self._uniformityAnalysis(cs, normuniformity)

    def fftAnalysis(self,cs,profile,mode='sensitivity'):
        """
        Fourier analysis: find main frequencies and power therein
        """
        if cs.verbose:
            try:
                import QCUS_testing as mytest
                if mode == 'sensitivity':
                    mytest._exportProfile( profile,fname='y2profile.tsv' )
                else:
                    mytest._exportProfile( profile,fname='x2profile.tsv' )
            except:
                pass

        fdelta = .1 # peak should differ by more than this fraction of max ps

        if mode == 'sensitivity':
            cs.sens_basewavelength = 0.
            cs.sens_ylim2 = -1

        # remove average component
        ywork = profile-np.average(profile)
        nx = len(ywork)
        x = np.array(range(nx))
        
        ffty = np.fft.rfft(ywork)
        ps = np.abs(ffty)**2
        integral = np.sum(ps[self.fft_skip:]) # calculate total content, but ignore zero component
        freqs = rfftfreq(nx)
            
        # find peaks
        xy_max,xy_min = wadwrapper_lib.peakdet(ps, delta=fdelta*np.max(ps[self.fft_skip:])) # returns indices for freqs

        # sort peaks from large to small
        xy_max = sorted(xy_max,key=operator.itemgetter(1))
        xy_max = xy_max[::-1]

        # find max component which is not zero-component and represents larger than given fraction of power
        base_ix = 0   # index of base-frequency
        start_ix = 0  # index of zero-peak
        for i,(xi,yi) in enumerate(xy_max):
            if xi<=self.fft_skip: # zero-peak must be located in skip
                start_ix = xi
            if xi>self.fft_skip and base_ix<1:
                base_ix = xi 
                if mode == 'sensitivity':
                    cs.sens_basewavelength = self.pixels2mm(cs, 1./freqs[xi])
                break

        # filter-out all non-background trend components to find extend of signal
        for i in range(len(ffty)):
            if i>(base_ix+start_ix)/2:
                ffty[i] = 0.

        # locate minimum (after zero-peak) of trend
        fy = np.fft.irfft(ffty)
        ix_max = np.argmax(fy)
        ix_min = ix_max+np.argmin(fy[ix_max:])
        if cs.verbose:
            print('max @',ix_max)
            print('min @',ix_min)

        if mode == 'sensitivity':
            cs.sens_ylim2 = ix_min
        
        if cs.verbose:
            plt.figure()
            plt.plot(ywork)
            plt.plot(fy)
            if mode == 'sensitivity':
                plt.title('sensitivity,filtered')
            else:
                plt.title('uniformity,filtered')
            plt.figure()
            plt.plot(ps)
            if mode == 'sensitivity':
                plt.title('sensitivity,powerspectrum')
            else:
                plt.title('uniformity,powerspectrum')
            cs.hasmadeplots = True

    def sensitivityAnalysis(self,cs):
        """
        Workflow:
        1. Calculate profile (vertically)
        2. Determine noise level
        3. Exponential fit to peaks in profile
        4. total sensitivity = penetration depth (intercept fit and noise)
        5. Repeat for subsets
        """
        error = True
    
        ## 1. Calculate profile (vertically)
        wid,hei = np.shape(cs.rect_image)
        offsetHOR= int(wid/4) #10 ; adjusted to remove influence of annotations through reverb data (enhancement line)

        if offsetHOR == 0:
            offsetHOR = self.offsetHOR
        crop = cropImage(cs.rect_image, [offsetHOR], [self.offsetVER])
        wid,hei = np.shape(crop)
        
        sensitivity = np.mean(crop,axis=0)
        time0 = time.time()
        self.fftAnalysis(cs,sensitivity,mode='sensitivity')
        cs.report.append(('sens_fft',time.time()-time0))

        nx = len(sensitivity)
        sensitivitysmoothed = mymath.movingaverage(sensitivity,self.smooth_sensitivity)
        if cs.testing:
            try:
                import QCUS_testing as mytest
                mytest._exportProfile( [ (s,ss) for s,ss in zip(sensitivity,sensitivitysmoothed) ],fname='yprofile.tsv' )
            except:
                pass
        
        ##2. Determine noise level (mean of last n values of sensitivity)
        cs.sens_noiseM = -1
        noiseRange = 2
        noise_av = np.mean(sensitivitysmoothed[-noiseRange:])
        noise_sd = np.std(sensitivitysmoothed[-noiseRange:])
        noise_snr = noise_av/noise_sd if noise_sd >0. else 0.
        noise_inc = False # sometimes snr starts with a local max
        for nr in range(noiseRange+2,int(nx/3)):
            av = np.mean(sensitivitysmoothed[-nr:])
            sd = np.std(sensitivitysmoothed[-nr:])
            snr = av/sd if sd >0. else 0.
            if snr>noise_snr or not noise_inc:
                if snr>noise_snr:
                    noise_inc = True 
                noise_av = av
                noise_sd = sd
                noise_snr = snr
            else:
                cs.sens_noiserange = nr
                break
        cs.sens_noiseM = noise_av

        top_val = np.max(sensitivitysmoothed)
        cut_val = cs.sens_noiseM+.1*(top_val-cs.sens_noiseM)
        for kk in reversed(range(nx)):
            if sensitivitysmoothed[kk]>cut_val:
                cs.sens_ylim = kk
                break
            
        # peak detection
        pmax = np.max(sensitivitysmoothed)
        pmin = np.min(sensitivitysmoothed)
        delta = (pmax-max(cs.sens_noiseM,pmin))*self.fdelta_sensitivity
        xy_max,xy_min = wadwrapper_lib.peakdet(sensitivitysmoothed, delta=delta,x=range(nx))
        
        # select only those peaks where max>noise; alternatively, select those peaks with min<noise
        if self.sens_hicut:
            xy_swap = []
            for xy in xy_max:
                if xy[1]>cs.sens_noiseM:
                    xy_swap.append(xy)
                else:
                    break
            xy_max = xy_swap
            if len(xy_max)>0:
                cs.sens_ylim = max(cs.sens_ylim,max(xy_max,key=operator.itemgetter(0))[0])

            xy_swap = []
            for xy in xy_min:
                if xy[0]<cs.sens_ylim:
                    xy_swap.append(xy)
                else:
                    xy_swap.append(xy) # last point
                    break
            xy_min = xy_swap
        else:
            xy_swap = []
            for xy in xy_min:
                if xy[1]<cs.sens_noiseM or len(xy_swap)<3:
                    xy_swap.append(xy)
                else:
                    break
            xy_min = xy_swap
            if len(xy_min)>1:
                cs.sens_ylim = max(xy_min,key=operator.itemgetter(0))[0]

            xy_swap = []
            for xy in xy_max:
                if xy[0]<cs.sens_ylim:
                    xy_swap.append(xy)
                else:
                    xy_swap.append(xy) # last point
                    break
            xy_max = xy_swap
        cs.sens_numtops = len(xy_max)
        cs.sens_numbots = len(xy_min)
        if cs.testing:
            try:
                import QCUS_testing as mytest
                mytest._exportProfile( xy_max, fname='xy_max.tsv' )
                mytest._exportProfile( xy_min, fname='xy_min.tsv' )
            except:
                pass
        
        if cs.verbose:
            x_max = np.array([xy[0] for xy in xy_max])
            y_max = np.array([xy[1] for xy in xy_max])
            x_min = np.array([xy[0] for xy in xy_min])
            y_min = np.array([xy[1] for xy in xy_min])

            x = np.array(range(nx))
            plt.figure()
            plt.plot(x,sensitivitysmoothed,'k-')
            plt.plot(x_min,y_min,'ro')
            plt.plot(x_max,y_max,'bo')
            cs.hasmadeplots = True


        error = False
        return error

    def Analyse(self,cs):
        """
        Analysis of reverberations in air.
        Workflow:
        1. Find reverberations as largest connected component != 0
        2. Straighten curved data if needed
        3. Uniformity of reverb analysis
        4. Sensitivity analysis
        """
        #return self.interpTest(cs)
        error = True
        
        time00 = time.time()
        cs.rect_image = None
        if cs.testing:
            try:
                import QCUS_testing as mytest
                cs.rect_image = mytest._importNDArray(fname='ndarray.npy')
    
            except Exception as e:
                print(e)

        if cs.rect_image is None:
            ## 1. CCA 
            # cluster connected components with pixelvalues>0
            error = self.isolateReverbrations(cs)
            if error:
                return error
    
            ## 2. Transform reverb data to rectangle is needed
            error = self.straightenCurve(cs)
            if error:
                return error
         
        time0 = time.time()
        error = self.sensitivityAnalysis(cs)
        cs.report.append(('sensitivity',time.time()-time0))
        if not error:
            time0 = time.time()
            error = self.reverbUniformity(cs)
            cs.report.append(('uniformity',time.time()-time0))

        cs.report.append(('total',time.time()-time00))

        print("#timing info#\npart seconds")
        for r,t in cs.report:
            print(r,t)
        return error
           
    def pixels2mm(self,cs,px):
        # translate pixels into mm

        dicomfields = [
            ["0018,6011, 0018,6024", "Physical Units X Direction"],
            ["0018,6011, 0018,602c", "Physical Delta X"],
        ] # Philips
    
        results = []
        units = self.readDICOMtag(cs,dicomfields[0][0])
        if units != 3: # 3 = cm
            return -1

        mm = 10.*self.readDICOMtag(cs,dicomfields[1][0]) # convert to mm
        return px*mm

    def imageID(self,cs,probeonly=False):
        if not probeonly and not cs.verbose:# and not self.guimode:
            return '0000'
        # make an identifier for this image
        if probeonly:
            di = self.DICOMInfo(cs,info='probe')
        else:
            di = self.DICOMInfo(cs,info='id')
        label = ''
        for k,v in di:
            label += '%s_'%v
        label = label[:-1]
        
        forbidden = '[,]\'" '
        label2 = ''
        for la in label:
            if la in forbidden:
                continue
            else:
                label2 += la
        label2 = label2.replace('UNUSED', '') # cleaning
        return label2

    def reportEntries(self,cs):
        # Convenience function to get all entries to report
        cs_items =[
            ('Curve_X',cs.curve_xyr[0]),
            ('Curve_Y',cs.curve_xyr[1]),
            ('Curve_R',cs.curve_xyr[2]),
            ('Curve_Residu',cs.curve_residu),
            ('Curve_OpenDeg',cs.curve_opendeg),
            ('Rect_width',np.shape(cs.rect_image)[0]),
            ('Sens_noiserange',cs.sens_noiserange),
            ('Sens_depth',cs.sens_ylim),
            ('Sens_fft_depth',cs.sens_ylim2),
            ('Sens_noise',cs.sens_noiseM),
            ('Sens_tops',cs.sens_numtops),
            ('Sens_bots',cs.sens_numbots),
            ('Sens_basewl_mm',cs.sens_basewavelength),
            ('Unif_bots',len(cs.unif_bots)),
            ('Unif_lowest',cs.unif_lowest),
            ('Unif_low',cs.unif_low),
            ('Unif_line',cs.unif_line),
        ]
        
        return cs_items
    
