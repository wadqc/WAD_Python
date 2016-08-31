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
input: dicomobject,pixeldata (the pixeldata is already cleaned; color is removed, leaving plain grayscale)
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
    20160831: bugfix no clusters left; Unified pywad1.0 and wad2.0; Transfer input parameters to xml; 
              Modifications in analysis paramters optimized for GE Voluson [PvH]; changed uniformity 
              to only analyse first bright reverb if possible
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
__version__ = '20160831'
__author__ = 'aschilham'

import copy
import numpy as np
import matplotlib.pyplot as plt
import operator
import scipy
import time
import os
import QCUS_math as mymath

# First try if we are running wad1.0, since in wad2 libs are installed systemwide
try: 
    # try local folder
    import wadwrapper_lib
except ImportError:
    # try pyWADlib from plugin.py.zip
    try: 
        from pyWADLib import wadwrapper_lib

    except ImportError: 
        # wad1.0 solutions failed, try wad2.0
        try: 
            # try system package wad_qc
            from wad_qc.modulelibs import wadwrapper_lib
        except ImportError: 
            # use parent wad_qc folder, and add it to search path
            import sys
            # add root folder of WAD_QC to search path for modules
            _modpath = os.path.dirname(os.path.abspath(__file__))
            while(not os.path.basename(_modpath) == 'Modules'):
                _new_modpath = os.path.dirname(_modpath)
                if _new_modpath == _modpath:
                    raise
                _modpath = _new_modpath
            sys.path.append(os.path.dirname(_modpath))
            from wad_qc.modulelibs import wadwrapper_lib

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
    unif_stdnorm  = 0 # Standard deviation over norm uniformity data
    unif_COV      = 0 # Coefficient of Variation over uniformity data
    unif_yrange   = [0,0] # Range in y that is analyzed for uniformity

    dipfarea_0_10   = 0 # fraction of image part of dip*depth for  0-10% width
    dipfarea_10_30  = 0 # fraction of image part of dip*depth for 10-30% width
    dipfarea_30_70  = 0 # fraction of image part of dip*depth for 30-70% width
    dipfarea_70_90  = 0 # fraction of image part of dip*depth for 70-90% width
    dipfarea_90_100 = 0 # fraction of image part of dip*depth for 90-100% width

    dipstrength_0_10   = 0 # number of dips*depth for  0-10% width
    dipstrength_10_30  = 0 # number of dips*depth for 10-30% width
    dipstrength_30_70  = 0 # number of dips*depth for 30-70% width
    dipstrength_70_90  = 0 # number of dips*depth for 70-90% width
    dipstrength_90_100 = 0 # number of dips*depth for 90-100% width
    
    lowfrac_0_10   = 0 # fraction of image part with norm.uniformity 5% or more below line average for  0-10% width
    lowfrac_10_30  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 10-30% width
    lowfrac_30_70  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 30-70% width
    lowfrac_70_90  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 70-90% width
    lowfrac_90_100 = 0 # fraction of image part with norm.uniformity 5% or more below line average for 90-100% width

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
    sens_bots = [] # store the y locations of the dips to use for uniformity analysis [px]
    
    # images
    image_fnames = [] # filenames of generated images

    # for gui
    cca_image = None # connected component analysis
    report = None # list ('description',value) of timing information

    #params
    uni_filter = 5 # running average of this width before peak detection
    uni_delta = 0.05 # A dip in normalized reverb pattern must be at least <delta> .05 = 5%
    sen_filter = 5 # running average of this width for sensitivity data
    sen_delta = 0.1 # A peak in sensitivity profile must be at least <fdelta>*(max-noise)
    ver_offset = 0 # default lines to exclude from top and bottom when making profiles (10)
    hor_offset = 0 # default rows to exclude from top and bottom when making profiles (10)
    fitcircle_frac = 1/3 # by default use only central 1/3 of circle for fitting, as deformations towards edge can occur. 
                            # use >1 for full fit. 1 is best for GE, 1/3 is best for Philips

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
        self.unif_stdnorm = 0
        self.unif_COV = 0

        self.dipfarea_0_10   = 0 # fraction of image part of dip*depth for  0-10% width
        self.dipfarea_10_30  = 0 # fraction of image part of dip*depth for 10-30% width
        self.dipfarea_30_70  = 0 # fraction of image part of dip*depth for 30-70% width
        self.dipfarea_70_90  = 0 # fraction of image part of dip*depth for 70-90% width
        self.dipfarea_90_100 = 0 # fraction of image part of dip*depth for 90-100% width
    
        self.dipstrength_0_10   = 0 # number of dips*depth for  0-10% width
        self.dipstrength_10_30  = 0 # number of dips*depth for 10-30% width
        self.dipstrength_30_70  = 0 # number of dips*depth for 30-70% width
        self.dipstrength_70_90  = 0 # number of dips*depth for 70-90% width
        self.dipstrength_90_100 = 0 # number of dips*depth for 90-100% width
        
        self.lowfrac_0_10   = 0 # fraction of image part with norm.uniformity 5% or more below line average for  0-10% width
        self.lowfrac_10_30  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 10-30% width
        self.lowfrac_30_70  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 30-70% width
        self.lowfrac_70_90  = 0 # fraction of image part with norm.uniformity 5% or more below line average for 70-90% width
        self.lowfrac_90_100 = 0 # fraction of image part with norm.uniformity 5% or more below line average for 90-100% width

        self.sens_ylim = 0
        self.sens_ylim2 = 0
        self.sens_noiseM = 0
        self.sens_numtops = 0
        self.sens_numbots = 0
        self.sens_noiserange = 2 
        self.sens_basewavelength = 0.
        self.sens_bots = []
        
        #input parameters
        self.uni_filter = 5
        self.uni_delta = 0.1
        self.sen_filter = 5
        self.sen_delta = 0.1
        self.ver_offset = 0
        self.hor_offset = 0
        self.fitcircle_frac = 1/3

        # helper stuff
        self.rev_miny = None 
        self.rev_maxy = None 
        self.rev_maxx = None 
        self.rev_minx = None 
        self.rev_midx = None 
        self.curve_angles = []
        self.curve_radii  = []
        self.rev_mask = None

class US_QC:
    qcversion = __version__
    fastmode = False
    modeLocalNorm = False

    # parameters:
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
            minsize = wid/10*hei/10
            while sum(cca.clusterSizes()> minsize)<2 and minsize>100:
                minsize = int(minsize/10)
            #cca.removeSmallClusters(wid/10*hei/10)
            cca.removeSmallClusters(minsize)
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

        if cs.verbose:
            scipy.misc.imsave('Reverbimage.jpg', cs.reverb_image)
        # bounds of reverb image
        time0 = time.time()
        if len(clus) == 0:
            error = True
        else:
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
            if cs.fitcircle_frac<1 and cs.fitcircle_frac>0:
                fff = 1.-cs.fitcircle_frac
                cf = mymath.CircleFit(circ_xy[int(fff*len(circ_xy)):])
            else:
                cf = mymath.CircleFit(circ_xy) # seems best for GE

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
            if cs.verbose:
                scipy.misc.imsave('Transfimage.jpg', cs.rect_image)
        else: # not a curved probe
            cs.rect_image = cs.reverb_image[cs.rev_minx[0]:cs.rev_maxx[0],cs.rev_miny[1]:cs.rev_maxy[1]]

        error = False
        return error

    def _uniformityAnalysis(self,cs, normuniformity, raw_uniformity):
        # Helper function for analysis of normalized uniformity profile:
        # normuniformity = normalized, smoothend uniformity
        # raw_uniformity = normalized uniformity
        error = True
        nx = len(normuniformity)
        time0 = time.time()
        self.fftAnalysis(cs,normuniformity,mode='uniformity')
        cs.report.append(('unif_fft',time.time()-time0))

        # peak detection with ugly hack to also allow peaks at index 0 and last
        #xy_max,xy_min = wadwrapper_lib.peakdet(normuniformity, delta=cs.uni_delta,x=range(nx))
        p0 = normuniformity[::-1]
        p2 = np.append(np.append(p0[:-1],normuniformity),p0[1:]) 
        x2 = np.array(range(-(nx-1),2*nx-1))
        xy_max,xy_min = wadwrapper_lib.peakdet(p2, delta=cs.uni_delta,x=x2)
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
                yran = [cs.ver_offset,ylim]

        if len(yran) == 0:
            yran = [cs.ver_offset,cs.sens_ylim] if cs.sens_ylim>0 else [cs.ver_offset]

        plt.figure()
        crop = cropImage(cs.rect_image, [cs.hor_offset], yran)
        plt.imshow(crop.transpose(),cmap=plt.gray())
        x = np.array(range(nx))
        plt.plot(x,-100*normuniformity,'y-',label='100*normalized')
        plt.plot(x_min,-100*y_min,'ro',label='accepted dips')
        plt.plot(x_max,-100*y_max,'bo',label='accepted peaks')
        plt.xlim([0,nx])

        if cs.verbose:
            np.savetxt('normuniformity.tsv', [ (s,ss) for s,ss in zip(raw_uniformity, normuniformity) ], delimiter='\t')
            np.savetxt('normuniformity_xy_max.tsv', xy_max, delimiter='\t')
            np.savetxt('normuniformity_xy_min.tsv', xy_min, delimiter='\t')

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
                delta = cs.uni_delta
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
                        delta -= .1*cs.uni_delta
                        if delta < .2*cs.uni_delta:
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
        cs.unif_stdnorm = np.std(normuniformity) # standard deviation of the normalised uniformity profile

        """
        uniformity_width = width of analysed profile [px]
        dip_width_i = width of dip_i as distance between zero-crossings
        frac_dip = sum_i (dip_width_i) / uniformity_width
        Ai = dip_width_i*dip_value_i
        calculate diparea in 10%-26%-28%-26%-10%
        WAD protocol: for L probe ignore outer 20% and action if dip deeper than 5 % of mean
        
        calculate fraction of region darker than 5% of line average

        """
        
        # calculate fraction of region part of dips weighted with depth of dip
        # number of dips in region weighted with depth of dip
        weights = [100/10, 100/20, 100/40, 100/20, 100/10] # 10%, 20%, 40%, 20%, 10%
        w_buckets = [0,0,0,0,0] # 10%, 20%, 40%, 20%, 10%
        p_buckets = [0,0,0,0,0] 
        for x,y in cs.unif_bots:
            if y<0: # skip positive dips
                width, w_contrib, p_contrib = self._dipanalysis(normuniformity, x)
                for i,(f,w,p) in enumerate(zip(weights, w_contrib, p_contrib)):
                    w_buckets[i]+= f*w/len(normuniformity)*-y #f*w/len(normuniformity)
                    p_buckets[i]+=-y*p
                    
        # calculate fraction of region darker than 5% of line average
        low = normuniformity<-0.05
        regions = [
            [int(0*nx),  int(.1*nx)],
            [int(.1*nx), int(.3*nx)],
            [int(.3*nx), int(.7*nx)],
            [int(.7*nx), int(.9*nx)],
            [int(.9*nx), nx],
        ]
        low_buckets = [0,0,0,0,0]

        for i, (x0,x1) in enumerate(regions):
            low_buckets[i] = np.average(low[x0:x1])
            
        if cs.verbose:
            print('  ','dip_frac','dips','low_frac')
            for w,p,low,avg,std,cov in zip(w_buckets, p_buckets,low_buckets):
                print('  ',w,p,low)

        # write to cs
        cs.dipfarea_0_10   = w_buckets[0]
        cs.dipfarea_10_30  = w_buckets[1]
        cs.dipfarea_30_70  = w_buckets[2]
        cs.dipfarea_70_90  = w_buckets[3]
        cs.dipfarea_90_100 = w_buckets[4]

        cs.dipstrength_0_10   = p_buckets[0]
        cs.dipstrength_10_30  = p_buckets[1]
        cs.dipstrength_30_70  = p_buckets[2]
        cs.dipstrength_70_90  = p_buckets[3]
        cs.dipstrength_90_100 = p_buckets[4]
        
        cs.lowfrac_0_10   = low_buckets[0]
        cs.lowfrac_10_30  = low_buckets[1]
        cs.lowfrac_30_70  = low_buckets[2]
        cs.lowfrac_70_90  = low_buckets[3]
        cs.lowfrac_90_100 = low_buckets[4]
        """
        meanvalue = np.mean(uniformity) 
        stdvalue = np.std(uniformity)
        CoefVar = stdvalue/meanvalue
        cs.unif_COV = CoefVar
        """
        
        error = False
        return error

    def _dipanalysis(self, normuniformity, y):
        # linear interpolate zero crossing on left of dip
        yleft = 0
        for yy in reversed(range(0,y)):
            if normuniformity[yy]>=0:
                yleft = yy+(0-normuniformity[yy])/(normuniformity[yy+1]-normuniformity[yy])
                break
        # linear interpolate zero crossing on right of dip
        yright = len(normuniformity)-1
        for yy in range(y,len(normuniformity)):
            if normuniformity[yy]>=0:
                y0 = yy-1
                yright = y0+(0-normuniformity[y0])/(normuniformity[y0+1]-normuniformity[y0])
                break
        width = yright-yleft
        
        # contributions
        profile_length = len(normuniformity)
        buckets = [ [0,profile_length*.1],
                    [profile_length*.1, profile_length*.3],
                    [profile_length*.3, profile_length*.7],
                    [profile_length*.7, profile_length*.9],
                    [profile_length*.9, profile_length-1],
                ]
        w_contrib = [0,0,0,0,0]
        p_contrib = [0,0,0,0,0]
        for i,edge in enumerate(buckets):
            if yleft<=edge[1] and yright>=edge[0]:
                w_contrib[i] += min(yright,edge[1])-max(yleft,edge[0])
            if y>=edge[0] and y<=edge[1]:
                p_contrib[i] += 1
        
        return width, w_contrib, p_contrib
    
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
        """
        As a fall back use 20% of the total y_range found from sensitivity analysis. 
        Prefer to analyse first full reverb band, as averaging over multiple bands will
        decrease the effect of a dip.
        """
        yran = []
        if self.fastmode:
            print('Uniformity: FASTMODE',cs.sens_basewavelength)
            if cs.sens_basewavelength >0:
                ylim = int(cs.sens_basewavelength/self.pixels2mm(cs,1.)+.5) 
                yran = [cs.ver_offset,ylim] # consider only first reverb length
            print(yran)
        if len(yran) == 0:
            yran = [cs.ver_offset,int(cs.sens_ylim/5)] if cs.sens_ylim>0 else [cs.ver_offset] # only first few reverbs needed
        if len(cs.sens_bots)>1:
            # find range between first good dip and second good dip.
            dist = []
            for i in range(1,len(cs.sens_bots)):
                dist.append(cs.sens_bots[i]-cs.sens_bots[i-1])
            dd = np.median(dist)/2 #at least start with top
            for i,yy in enumerate(cs.sens_bots):
                if yy<dd: continue
                yran = [yy, cs.sens_bots[i+1]]
                break

        cs.unif_yrange   = yran # Range in y that is analyzed for uniformity
                
        crop = cropImage(cs.rect_image, [cs.hor_offset], yran)

        uniformity = np.mean(crop,axis=1)
        
        #uniformitysmoothed = mymath.movingaverage(uniformity,cs.uni_filter)
        uniformitysmoothed = mymath.smooth(uniformity,cs.uni_filter, window='flat') # Smoothing without boundary effects

        ## line uniformity
        intemin = np.min(uniformity)
        intemax = np.max(uniformity)
        inteavg = np.mean(uniformity)
        cs.unif_line = np.max([intemax-inteavg,inteavg-intemin])/inteavg

        # Output of COV 
        meanvalue = np.mean(uniformity) 
        stdvalue = np.std(uniformity)
        CoefVar = stdvalue/meanvalue
        cs.unif_COV = CoefVar

        # normalized uniformity 
        if self.modeLocalNorm:
            frac = .1
            print("Uniformity, using Local Normalization",frac)
            normuniformity = mymath.localnormalization(uniformity, sigma=len(uniformity)*frac)
        else:
            normuniformity = ((uniformitysmoothed-meanvalue)/meanvalue)

        return self._uniformityAnalysis(cs, normuniformity, (uniformity-meanvalue)/meanvalue )

    def fftAnalysis(self,cs,profile,mode='sensitivity'):
        """
        Fourier analysis: find main frequencies and power therein
        """
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
        xy_max = sorted(xy_max,key=operator.itemgetter(1), reverse=True)

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
            offsetHOR = cs.hor_offset
        crop = cropImage(cs.rect_image, [offsetHOR], [cs.ver_offset])
        wid,hei = np.shape(crop)

        sensitivity = np.mean(crop,axis=0)
        time0 = time.time()
        self.fftAnalysis(cs,sensitivity,mode='sensitivity')
        cs.report.append(('sens_fft',time.time()-time0))

        nx = len(sensitivity)
        #sensitivitysmoothed = mymath.movingaverage(sensitivity,cs.sen_filter)
        sensitivitysmoothed = mymath.smooth(sensitivity,cs.sen_filter, window='flat') # Smoothing without boundary effects

        if cs.verbose:
            np.savetxt('sensitivity.tsv', [ (s,ss) for s,ss in zip(sensitivity,sensitivitysmoothed) ], delimiter='\t')

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
        delta = (pmax-max(cs.sens_noiseM,pmin))*cs.sen_delta
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
        cs.sens_bots = [x for x,y in xy_min] # store the y locations of the dips to use for uniformity analysis
        
        if cs.verbose:
            np.savetxt('sensitivity_xy_max.tsv', xy_max, delimiter='\t')
            np.savetxt('sensitivity_xy_min.tsv', xy_min, delimiter='\t')

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

        if cs.verbose:
            scipy.misc.imsave('Inputdata.jpg', cs.pixeldataIn)

        time00 = time.time()
        cs.rect_image = None

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
            ('Unif_stdnorm',cs.unif_stdnorm), 
            ('Unif_COV',cs.unif_COV),
            ('Unif_ymin',cs.unif_yrange[0]),
            ('Unif_ymax',cs.unif_yrange[1]),
            ('DipFracArea_0_10',cs.dipfarea_0_10),
            ('DipFracArea_10_30',cs.dipfarea_10_30),
            ('DipFracArea_30_70',cs.dipfarea_30_70),
            ('DipFracArea_70_90',cs.dipfarea_70_90),
            ('DipFracArea_90_100',cs.dipfarea_90_100),
            ('DipStrength_0_10',cs.dipstrength_0_10),
            ('DipStrength_10_30',cs.dipstrength_10_30),
            ('DipStrength_30_70',cs.dipstrength_30_70),
            ('DipStrength_70_90',cs.dipstrength_70_90),
            ('DipStrength_90_100',cs.dipstrength_90_100),
            ('LowFrac_0_10',cs.lowfrac_0_10),
            ('LowFrac_10_30',cs.lowfrac_10_30),
            ('LowFrac_30_70',cs.lowfrac_30_70),
            ('LowFrac_70_90',cs.lowfrac_70_90),
            ('LowFrac_90_100',cs.lowfrac_90_100),
        ]

        return cs_items
