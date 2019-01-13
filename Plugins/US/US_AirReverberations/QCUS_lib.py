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
4. Uniformity analysis of rect_image to find dips that suggest problems with the probe:
  a. Make a profile along the horizontal of the rect_image upto sensitivity_extend;
     FASTMODE: (experimental) use onle rect_image of first ring; should have all the info there already
  b. Detect dips that are deep enough, find the minimum in of the dips, and overal and the overal non-uniformity


Warning: THIS MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

TODO:
    o Maybe put hard threshold on peak value for uniformity (is normalized, so why not?!)
Changelog:
    20180731: fix error ValueError: assignment destination is read-only for work[cs.pixeldataIn ==0] = 1
    20170912: tried sensitivity profile analysis (per column, calculate depth relative to accepted depth; 
              then calculate COV skew kurt); does not help; removed again.
              tried combining DipFrac etc to 0-10+90-100, 10-30+70-90; does not help; removed again.
    20170909: remove deviation, add 10-90, 30-70 versions of kurtosis, COV, skew
    20170906: testing kurtosis, dev=max-min, and per region;
              adding alternative uni_range definitions; 
              error in calculating straightend image is build from several connected components.
    20170830: added boundingbox parameter to fix extend of reverbpattern (primarily to exclude know P problem philips); added auto_suffix
    20170510: extra parameters uni_start, clustermode; extra measurement skewness; 
              removed unused uniformity yrange stuff; removed unused fft stuff; removed localnormalization; removed fastmode;
              removed helper cropImage; added makeFigure for uniformity and rois; removed 'reverb' image;
              remove timing
    20170109: remove unprintable characters from DICOM info
    20161221: added level to reportvalues; added unif_depth and unif_sum; normalize on median;
              isolateReverbrations threshold is config param (PvH)
    20161220: removed class variables; removed testing stuff
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
__version__ = '20180731'
__author__ = 'aschilham, pvanhorsen'

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import operator
import scipy
import scipy.stats
import os

try:
    # wad2.0 runs each module stand alone
    import QCUS_math as mymath
except ImportError:
    from . import QCUS_math as mymath
    
from PIL import Image # image from pillow is needed
from PIL import ImageDraw # imagedraw from pillow is needed, not pil

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

class USStruct:
    def __init__ (self, dcmInfile, pixeldataIn, dicomMode):
        self.verbose = False

        # input image
        self.dcmInfile = dcmInfile
        self.pixeldataIn = pixeldataIn
        self.dicomMode = dicomMode
        self.resultslabel = None # allow override of label by config

        # for matlib plotting
        self.hasmadeplots = False

        # for gui
        self.report = [] # list ('description',value) of timing information
        self.cca_image = None # connected component analysis

        # images
        self.image_fnames = [] # filenames of generated images

        # reverbrations
        self.cluster_fminsize = 10*10*3 # ignore clusters of size smaller than imwidth*imheigth/minsizefactor (wid/10*hei/10)/3
        self.reverb_image = None        # isolated reverbrations (by largest connected component)
        self.rect_image = None          # straightend image for curvilinear probes
        self.signal_thresh = 0          # only pixelvalues above this number can be part of reverberations (set >0 is very noisy)
        self.cluster_mode = 'all_middle' # 'largest_only' default mode of dataroi selection

        self.curve_xyr = [0,0,0] # center and radius [x,y,r] of curvilinear probe
        self.curve_residu = 0    # residu of fitting of curvilinear probe
        self.curve_opendeg = 0   # opening angle in degrees


        self.unif_bots = []      # list of (pos,val) of kept valleys in normalized reverb profile
        self.unif_lowest = 0     # lowest value of valleys
        self.unif_low = 0        # lowest value 
        self.unif_line = 1       # Overal non-uniformity over whole, unnormalized pofile
        self.unif_stdnorm = 0    # Standard deviation over norm uniformity data
        self.unif_COV = 0        # Coefficient of Variation over uniformity data
        self.unif_yrange = [0,0] # Range in y that is analyzed for uniformity
        self.unif_sum = 0        # Sum of normalised uniformity values
        self.unif_skew = 0       # skewness of the distribution
        self.unif_kurtosis = 0   # kurtosis of the distribution

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

        self.unif_COV_10_90 = 0        # Coefficient of Variation over uniformity data 10-90%
        self.unif_skew_10_90 = 0       # skewness of the distribution
        self.unif_kurtosis_10_90 = 0   # kurtosis of the distribution
        self.unif_COV_30_70 = 0        # Coefficient of Variation over uniformity data 30-70%
        self.unif_skew_30_70 = 0       # skewness of the distribution
        self.unif_kurtosis_30_70 = 0   # kurtosis of the distribution

        # sensitivity
        self.sens_ylim = 0    # limit on y determined from peaks in y-dir 
        self.sens_noiseM = 0  # noise determined from (vertical) sensitivity profile
        self.sens_numtops = 0 # number of peaks found upto ylim
        self.sens_numbots = 0 # number of valleys found upto ylim
        self.sens_noiserange = 2      # number of lines used for calculation of noiselevel
        self.sens_bots = []   # store the y locations of the dips to use for uniformity analysis [px]
        
        #input parameters
        self.uni_filter = 5  # running average of this width before peak detection
        self.uni_delta = 0.1 # A dip in normalized reverb pattern must be at least <delta> .05 = 5%
        self.uni_depth = 5   # default depth in mm for uniformity analysis
        self.uni_start = 2   # default offset in mm to start of uniformity analysis
        self.sen_filter = 5  # running average of this width for sensitivity data
        self.sen_delta = 0.1 # A peak in sensitivity profile must be at least <fdelta>*(max-noise)
        self.ver_offset = 0  # default lines to exclude from top and bottom when making profiles (10)
        self.hor_offset = 0  # default rows to exclude from top and bottom when making profiles (10)
        self.fitcircle_frac = 1/3 # by default use only central 1/3 of circle for fitting, as deformations towards edge can occur. 
                                  # use >1 for full fit. 1 is best for GE, 1/3 is best for Philips

        self.uni_range_model = 'absolute' # restrict analysis from uni_start to uni_start+uni_depth (if uni_depth > sens_ylim, this will include many 0s)
        #self.uni_range_model = 'skip10pct'# restrict analysis from 0.1*sens_ylim to 0.9*sens_ylim
        #self.uni_range_model = 'skip20pct'# restrict analysis from 0.2*sens_ylim to 0.8*sens_ylim
        #self.uni_range_model = 'maxsenslimit' # restrict analysis from uni_start to min(sens_ylim, uni_start+uni_depth)

        # helper stuff
        self.rev_forcebbox = None # (xmin, xmax, ymin, ymax) in px to be used to restrict reverb pattern
        self.rev_miny = None # (x,y) of position of minimum y of reverb image
        self.rev_maxy = None # (x,y) of position of maximum y of reverb image
        self.rev_maxx = None # (x,y) of position of maximum x of reverb image
        self.rev_minx = None # (x,y) of position of minimum x of reverb image
        self.rev_midx = None # estimated x of middle of reverb image
        self.rev_mask = None # mask of reverb image
        self.curve_angles = [] # [min, max] radius
        self.curve_radii  = [] # [min, max] angle
        self.auto_suffix = True # add a probename as suffix to all results

class US_QC:
    def __init__(self, guimode=False):
        self.qcversion = __version__
        self.guimode = guimode

        # parameters:
        self.sens_hicut = True # restrict peaks to max > noise (True) or min < noise

    def readDICOMtag(self, cs, key, imslice=0): # slice=2 is image 3
        value = wadwrapper_lib.readDICOMtag(key, cs.dcmInfile, imslice)
        return value

    def DICOMInfo(self, cs, info='dicom'):
        # Different from ImageJ version; tags "0008","0104" and "0054","0220"
        #  appear to be part of sequences. This gives problems (cannot be found
        #  or returning whole sequence blocks)
        # Possibly this can be solved by using if(type(value) == type(dicom.sequence.Sequence()))
        #  but I don't see the relevance of these tags anymore, so set them to NO

        import string
        printable = set(string.printable)
        
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
                value = str(self.readDICOMtag(cs,key)).replace('&','')
                value = ''.join(list(filter(lambda x: x in printable, value)))
            except:
                value = ""
            results.append( (df[1],value) )

        return results

    #----------------------------------------------------------------------
    def isolateReverbrations(self, cs):
        """
        Find reverbrations part of image.
        Workflow:
        1. Restrict to bbox if provided
        2. Find reverberations as largest connected component != 0
        2. Return
        """
        error = True
        # cluster connected components with pixelvalues>0
        #work = (cs.pixeldataIn>0) * (cs.pixeldataIn<255) # must give a signal, but 255 often reserved for overlay
        work = cs.pixeldataIn>cs.signal_thresh

        # restrict to bbox if provided:
        if not cs.rev_forcebbox is None:
            xmin,xmax,ymin,ymax = cs.rev_forcebbox
            work[    :xmin,     :    ] = 0
            work[xmax:    ,     :    ] = 0
            work[    :    ,     :ymin] = 0
            work[    :    , ymax:    ] = 0
            
        cca = wadwrapper_lib.connectedComponents()
        cs.cca_image,nb_labels = cca.run(work)
        if cs.cluster_mode == 'largest_only': # model of PVH
            # select only largest cluster
            cluster_sizes = cca.clusterSizes()
            clus_val = np.argmax(cluster_sizes)
            cs.rev_mask = (cs.cca_image == clus_val)
            clus = cca.indicesOfCluster(clus_val)
        else: #'all_middle'
            # select all clusters present in vertical middle area of image, excluding top and 0
            wid,hei = np.shape(cs.cca_image)
            # first remove very small clusters (can be located around top edge!)
            minsize = wid*hei*cs.cluster_fminsize #(wid/10*hei/10)/3
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

        if cs.verbose:
            scipy.misc.imsave('Reverbimage.jpg', cs.reverb_image)

        # bounds of reverb image
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

            # transform reverb pattern to rectangle: interpolate at coords
            ang = np.linspace(cs.curve_angles[0],cs.curve_angles[1],cs.rev_maxx[0]-cs.rev_minx[0])
            #rad = np.linspace(cs.curve_radii[0],cs.curve_radii[1],cs.rev_maxy[1]-cs.rev_miny[1])
            rad = np.linspace(cs.curve_radii[0],cs.curve_radii[1],int(0.5+cs.curve_radii[1]-cs.curve_radii[0]))
            an,ra = scipy.meshgrid(ang,rad)
            xi = cs.curve_xyr[0]+ra*np.sin(an)
            yi = cs.curve_xyr[1]+ra*np.cos(an)
            coords = np.array( [xi,yi] )
            cs.rect_image = scipy.ndimage.map_coordinates(cs.reverb_image, coords).transpose()
            if cs.verbose:
                scipy.misc.imsave('Transfimage.jpg', cs.rect_image)
        else: # not a curved probe
            cs.rect_image = cs.reverb_image[cs.rev_minx[0]:cs.rev_maxx[0],cs.rev_miny[1]:cs.rev_maxy[1]]

        error = False
        return error

    def _uniformityAnalysis(self,cs, normuniformity, raw_uniformity):
        # Helper function for analysis of normalized uniformity profile: (find dips and dark areas)
        # normuniformity = normalized, smoothend uniformity
        # raw_uniformity = normalized uniformity
        error = True
        nx = len(normuniformity)

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

        # generate an image showing detected dips and peaks
        fig = self.makeFigure(cs, 'uniformity', 
                              {'x_min': x_min, 'x_max': x_max, 'y_min': y_min, 'y_max': y_max,
                               'normuniformity': normuniformity})
        if not self.guimode:
            fname = 'uniformity'
            suffix = self.imageID(cs,probeonly=True)
            if not suffix == '':
                fname += '_'+suffix
            fname = fname +'.jpg'
            fig.savefig(fname)
            cs.image_fnames.append(fname)
        else:
            cs.hasmadeplots = True
        

        if cs.verbose:
            np.savetxt('normuniformity.tsv', [ (s,ss) for s,ss in zip(raw_uniformity, normuniformity) ], delimiter='\t')
            np.savetxt('normuniformity_xy_max.tsv', xy_max, delimiter='\t')
            np.savetxt('normuniformity_xy_min.tsv', xy_min, delimiter='\t')

        if cs.verbose:
            plt.figure()
            x = np.array(range(nx))
            plt.plot(x,normuniformity,'k-')
            plt.plot(x_min,y_min,'ro')
            plt.plot(x_max,y_max,'bo')
            cs.hasmadeplots = True

        cs.unif_bots = xy_min
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
        weights = [100./10, 100./20, 100./40, 100./20, 100./10] # 10%, 20%, 40%, 20%, 10%
        w_buckets = [0,0,0,0,0] # 10%, 20%, 40%, 20%, 10%
        p_buckets = [0,0,0,0,0] 
        for x,y in cs.unif_bots:
            if y<0: # skip positive dips
                width, w_contrib, p_contrib = self._dipanalysis(normuniformity, x)
                for i,(f,w,p) in enumerate(zip(weights, w_contrib, p_contrib)):
                    w_buckets[i]+= f*w/len(normuniformity)*-y #f*w/len(normuniformity)
                    p_buckets[i]+=-y*p

        # calculate fraction of region darker than 2* stddev of line average
        low = normuniformity # low parameter, now used as temporary profile for low analysis per bucket.

        regions = [
            [int(0*nx),  int(.1*nx)],
            [int(.1*nx), int(.3*nx)],
            [int(.3*nx), int(.7*nx)],
            [int(.7*nx), int(.9*nx)],
            [int(.9*nx), nx],
        ]
        low_buckets = [0,0,0,0,0]
        cs.unif_sum = np.sum(low)
        if cs.unif_sum < 0:
            cs.unif_sum = np.abs(cs.unif_sum)
        else:
            cs.unif_sum = 0

        for i, (x0,x1) in enumerate(regions):
            if np.min(low[x0:x1])<0:
                low_buckets[i] = -100*np.min(low[x0:x1])# lowest value per bucket; multiply with -100 to yield percentage
            
        if cs.verbose:
            print('  ','dip_frac','dips','low_frac')
            for w,p,lw in zip(w_buckets, p_buckets,low_buckets):
                print('  ',w,p,lw)

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
        
        if cs.verbose:
            plt.figure()
            print("region: var std cov skew")
            lab = "all"
            reg = low
            hist,bins = np.histogram(reg, bins=100)
            #plt.plot(sorted(reg), label=lab)
            plt.plot(bins[1:], hist, label=lab)
            print("{}: {} {} {} {} {}".format(lab, (np.max(reg)-np.min(reg)), np.std(reg), scipy.stats.variation(reg), scipy.stats.skew(reg), scipy.stats.kurtosis(reg)))
            for reg, lab, num in [
                (low[regions[0][0]:regions[0][1]], "00_10", 10),
                (low[regions[1][0]:regions[1][1]], "10_30", 20),
                (low[regions[2][0]:regions[2][1]], "30_70", 40),
                (low[regions[3][0]:regions[3][1]], "70_90", 20),
                (low[regions[4][0]:regions[4][1]], "90_00", 10),]:  
                hist,bins = np.histogram(reg, bins=num)
                #plt.plot(sorted(reg), label=lab)
                plt.plot(bins[1:], hist, label=lab)
                print("{}: {} {} {} {} {}".format(lab, (np.max(reg)-np.min(reg)), np.std(reg), scipy.stats.variation(reg), scipy.stats.skew(reg), scipy.stats.kurtosis(reg)))

            plt.legend()
            cs.hasmadeplots = True

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
        
        return width, w_contrib, p_contrib # width of dip, width of dip in buckets, peak in which buckets
    
    def reverbUniformity(self,cs):
        """
        Analysis uniformity of reverb pattern (find dips and dark areas)
        Workflow:
        1. Define ROI of reverb (dependend on vertical profile)
        2. Average ROI horizontally and calculate COV
        3. Normalize profile as deviation wrt mean
        4. Find dips in profile
        """
        error = True

        ## 1. Define ROI of reverb (dependend on vertical profile)
        # Take a fixed depth in mm (as a max). Testing, if OK, then remove stuff above
        depth  = int(float(cs.uni_depth)/self.pixels2mm(cs,1.)+.5) # depth from mm to pixels
        offset = int(float(cs.uni_start)/self.pixels2mm(cs,1.)+.5) # offset from mm to pixels

        if cs.uni_range_model == 'absolute':
            # restrict analysis from uni_start to uni_start+uni_depth (if uni_depth > sens_ylim, this will include many 0s)
            yran = [offset, offset+depth]
        elif cs.uni_range_model == 'skip10pct':
            # restrict analysis from 0.1*sens_ylim to 0.9*sens_ylim
            skip = .1 # fraction
            yran = [ int(skip*cs.sens_ylim), int((1.-skip)*cs.sens_ylim) ]
        elif cs.uni_range_model == 'skip20pct':
            # restrict analysis from 0.2*sens_ylim to 0.8*sens_ylim
            skip = .2 # fraction
            yran = [ int(skip*cs.sens_ylim), int((1.-skip)*cs.sens_ylim) ]
        elif cs.uni_range_model == 'maxsenslimit':
            # restrict analysis from uni_start to min(sens_ylim, uni_start+uni_depth)
            yran = [ offset, min(offset+depth, cs.sens_ylim) ]
        else:
            raise ValueError("Unknown uni_range_model '{}'".format(cs.uni_range_model))

        # crop image
        im_wid, im_hei = np.shape(cs.rect_image)
        xran = [cs.hor_offset, -cs.hor_offset]
        if xran[1] == 0: xran[1] = im_wid-1
        if yran[1] == 0: yran[1] = im_hei-1

        crop = cs.rect_image[ xran[0]:xran[1], yran[0]:yran[1] ]

        # store for reference 
        cs.unif_yrange   = yran # Range in y that is analyzed for uniformity

        uniformity = np.mean(crop,axis=1)
        uniformitysmoothed = mymath.smooth(uniformity,cs.uni_filter, window='flat') # Smoothing without boundary effects
        if cs.verbose:
            np.savetxt('uniformity.tsv', [ (s,ss) for s,ss in zip(uniformity, uniformitysmoothed) ], delimiter='\t')

        ## line uniformity
        intemin = np.min(uniformity)
        intemax = np.max(uniformity)
        inteavg = np.mean(uniformity)
        cs.unif_line = np.max([intemax-inteavg,inteavg-intemin])/inteavg

        # Output of COV based on smoothed uniformity
        meanvalue = np.mean(uniformitysmoothed) 
        medianvalue = np.median(uniformitysmoothed)
        stdvalue = np.std(uniformitysmoothed)
        CoefVar = stdvalue/meanvalue
        cs.unif_COV = CoefVar*100 # in pct
        cs.unif_skew = scipy.stats.skew(uniformitysmoothed)
        cs.unif_kurtosis = scipy.stats.kurtosis(uniformitysmoothed)

        wid = len(uniformitysmoothed)
        part = uniformitysmoothed[ int(0.1*wid):int(0.9*wid) ]
        cs.unif_COV_10_90      = scipy.stats.variation(part)*100 # in pct
        cs.unif_skew_10_90     = scipy.stats.skew(part)
        cs.unif_kurtosis_10_90 = scipy.stats.kurtosis(part)

        part = uniformitysmoothed[ int(0.3*wid):int(0.7*wid) ]
        cs.unif_COV_30_70      = scipy.stats.variation(part)*100 # in pct
        cs.unif_skew_30_70     = scipy.stats.skew(part)
        cs.unif_kurtosis_30_70 = scipy.stats.kurtosis(part)

        # normalized uniformity 
        normuniformity = ((uniformitysmoothed-medianvalue)/medianvalue)

        return self._uniformityAnalysis(cs, normuniformity, (uniformity-medianvalue)/medianvalue )

    def sensitivityAnalysis(self,cs):
        """
        Workflow:
        1. Calculate profile (vertically): count rings, determine noise level
        2. Determine noise level
        3. Exponential fit to peaks in profile
        4. total sensitivity = penetration depth (intercept fit and noise)
        5. Repeat for subsets
        """
        error = True

        ## 1. Calculate profile (vertically)
        im_wid,im_hei = np.shape(cs.rect_image)
        offsetHOR = int(im_wid/4.) #10 ; adjusted to remove influence of annotations through reverb data (enhancement line)

        if offsetHOR == 0:
            offsetHOR = cs.hor_offset
        
        # crop image
        xran = [offsetHOR, im_wid-1-offsetHOR]
        cs.sens_xrange = xran
        yran = [cs.ver_offset, -cs.ver_offset]
        if xran[1] == 0: xran[1] = im_wid-1
        if yran[1] == 0: yran[1] = im_hei-1
        crop = cs.rect_image[ xran[0]:xran[1], yran[0]:yran[1] ]
        wid,hei = np.shape(crop)

        sensitivity = np.median(crop,axis=0) # more robust for defects
        nx = len(sensitivity)
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
            plt.plot([kk],[0],'c+')
            plt.plot([cs.sens_ylim],[0],'cx')
            plt.plot([0, nx-1], [cs.sens_noiseM, cs.sens_noiseM], 'c:')
            plt.plot([0, nx-1], [cut_val, cut_val], 'c-')

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

        error = self.sensitivityAnalysis(cs)
        if not error:
            error = self.reverbUniformity(cs)

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

    def imageID(self, cs, probeonly=False):
        if not cs.resultslabel is None:
            return cs.resultslabel
        
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
        cs.resultslabel = label2.replace('/','-')
        return label2

    def reportEntries(self,cs):
        # Convenience function to get all entries to report
        cs_items =[ # name, value, level
            ('Curve_X',cs.curve_xyr[0], 2),
            ('Curve_Y',cs.curve_xyr[1], 2),
            ('Curve_R',cs.curve_xyr[2], 2),
            ('Curve_Residu',cs.curve_residu, 2),
            ('Curve_OpenDeg',cs.curve_opendeg, 2),
            ('Rect_width',np.shape(cs.rect_image)[0], 2),
            ('Sens_noiserange',cs.sens_noiserange, 2),
            ('Sens_depth',self.pixels2mm(cs, cs.sens_ylim), 1),
            ('Sens_noise',cs.sens_noiseM, 2),
            ('Sens_tops',cs.sens_numtops, 2),
            ('Sens_bots',cs.sens_numbots, 2),
            ('Unif_bots',len(cs.unif_bots), 1),
            ('Unif_lowest',cs.unif_lowest, 2),
            ('Unif_low',cs.unif_low, 2),
            ('Unif_line',cs.unif_line, 2),
            ('Unif_stdnorm',cs.unif_stdnorm, 2), 
            ('Unif_sumlow',cs.unif_sum, 1),
            ('Unif_COV',cs.unif_COV, 1),
            ('Unif_skew',cs.unif_skew, 1),
            ('Unif_ymin',cs.unif_yrange[0], 2),
            ('Unif_ymax',cs.unif_yrange[1], 2),
            ('DipFracArea_0_10',cs.dipfarea_0_10, 2),
            ('DipFracArea_10_30',cs.dipfarea_10_30, 2),
            ('DipFracArea_30_70',cs.dipfarea_30_70, 2),
            ('DipFracArea_70_90',cs.dipfarea_70_90, 2),
            ('DipFracArea_90_100',cs.dipfarea_90_100, 2),
            ('DipStrength_0_10',cs.dipstrength_0_10, 2),
            ('DipStrength_10_30',cs.dipstrength_10_30, 2),
            ('DipStrength_30_70',cs.dipstrength_30_70, 2),
            ('DipStrength_70_90',cs.dipstrength_70_90, 2),
            ('DipStrength_90_100',cs.dipstrength_90_100, 2),
            ('LowFrac_0_10',cs.lowfrac_0_10, 1),
            ('LowFrac_10_30',cs.lowfrac_10_30, 1),
            ('LowFrac_30_70',cs.lowfrac_30_70, 1),
            ('LowFrac_70_90',cs.lowfrac_70_90, 1),
            ('LowFrac_90_100',cs.lowfrac_90_100, 1),

            ('Unif_kurtosis',cs.unif_kurtosis, 1),

            ('Unif_COV_10_90',cs.unif_COV_10_90, 1),
            ('Unif_skew_10_90',cs.unif_skew_10_90, 1),
            ('Unif_kurtosis_10_90',cs.unif_kurtosis_10_90, 1),
            ('Unif_COV_30_70',cs.unif_COV_30_70, 1),
            ('Unif_skew_30_70',cs.unif_skew_30_70, 1),
            ('Unif_kurtosis_30_70',cs.unif_kurtosis_30_70, 1),
        ]

        return cs_items


    def makeFigure(self, cs, what, xtra={}):
        """
        Use matplotlib to make an analysis image indicating all analyzed and interesting parts
        """
        if what == 'uniformity':
            # generate an image showing detected dips and peaks
            # crop image
            fig = plt.figure(facecolor='black') # black frame
            ax = plt.gca()
            ax.set_facecolor('black')
            
            im_wid, im_hei = np.shape(cs.rect_image)
            xran = [cs.hor_offset, -cs.hor_offset]
            #yran = [cs.ver_offset, cs.sens_ylim] if cs.sens_ylim>0 else [cs.ver_offset, -cs.ver_offset]
            yran = [cs.ver_offset, -cs.ver_offset]
            if xran[1] == 0: xran[1] = im_wid-1
            if yran[1] == 0: yran[1] = im_hei-1
            crop = cs.rect_image[ xran[0]:xran[1], yran[0]:yran[1] ]
            plt.imshow(crop.transpose(),cmap=plt.gray())
            plt.axis('off')
            
            # draw boxes around dip calc areas
            c_wid, c_hei = np.shape(crop)
            rectrois = [] # x0,y0, width, height, color, alpha
            parts = [ [0,.1], [.1,.3], [.3,.7], [.7,.9], [.9,1] ]
            for p0,p1 in parts:
                rectrois.append( (p0*c_wid, cs.unif_yrange[0],  # (x,y)
                                  c_wid*(p1-p0),          # width
                                  cs.unif_yrange[1]-cs.unif_yrange[0], # height
                                  'r', .5)
                                 )

            # add box around sensitivity part (only meaningfull on linear probe image)
            sens_xran = cs.sens_xrange
            if sens_xran[0] != xran[0]: # correct if part only
                sens_xran[0] -= xran[0]
                sens_xran[1] -= xran[0]
            rectrois.append( (sens_xran[0], yran[0], 
                              sens_xran[1]-sens_xran[0], cs.sens_ylim, 'c', .5 ) ) # correct! sens_ylim + yran is end
        

            for x0,y0,w,h,c,a in rectrois:
                ax.add_patch(Rectangle( 
                    (x0, y0),   # (x,y)
                    w, # width
                    h, # height
                    fill = None, color = c, alpha=a)
                )
                
            # add enhanced dips line on top
            normuniformity = xtra['normuniformity']
            nx = len(normuniformity)
            x_min = xtra['x_min']
            x_max = xtra['x_max']
            y_min = xtra['y_min']
            y_max = xtra['y_max']

            # ax1 = top left
            yoff = max(10, int(.1*(yran[1]-yran[0])))
            x = np.array(range(nx))
            plt.plot(x,-100*normuniformity-yoff,'y-',label='100*normalized')
            plt.plot(x_min,-100*y_min-yoff,'ro',label='accepted dips')
            plt.plot(x_max,-100*y_max-yoff,'co',label='accepted peaks')
            plt.xlim( [0, nx] )
            plt.tight_layout()
    
            return fig
            
    def saveAnnotatedImage(self, cs, fname, what='overview', xtra={}):
        """
        Make an jpg of the original image, indicating all analyzed and interesting parts
        """
        # make a palette, mapping intensities to greyscale
        pal = np.arange(0,256,1,dtype=np.uint8)[:,np.newaxis] * \
            np.ones((3,),dtype=np.uint8)[np.newaxis,:]
        # but reserve the first for red for markings
        pal[0] = [255,0,0]

        rectrois = []
        polyrois = []
        circlerois = []

        # convert to 8-bit palette mapped image with lowest palette value used = 1
        if what == 'overview':
            # first the base image
            work = np.array(cs.pixeldataIn)
            work[cs.pixeldataIn ==0] = 1
            im = scipy.misc.toimage(work.transpose(),low=1, pal=pal) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

            # add box around isolated reverb pattern
            minx = min([x for x,y in [cs.rev_minx, cs.rev_miny, cs.rev_maxx, cs.rev_maxy]])
            maxx = max([x for x,y in [cs.rev_minx, cs.rev_miny, cs.rev_maxx, cs.rev_maxy]])
            miny = min([y for x,y in [cs.rev_minx, cs.rev_miny, cs.rev_maxx, cs.rev_maxy]])
            maxy = max([y for x,y in [cs.rev_minx, cs.rev_miny, cs.rev_maxx, cs.rev_maxy]])

            # add box around reverb region
            rectrois.append( [(minx, miny),(maxx, maxy)] )

        # add extra rois if provided
        if 'circlerois' in xtra:
            for r in xtra['circlerois']:
                circlerois.append(r)
        if 'polyrois' in xtra:
            for r in xtra['polyrois']:
                polyrois.append(r)
        if 'rectrois' in xtra:
            for r in xtra['rectrois']:
                rectrois.append(r)

        # now draw all rois in reserved color
        draw = ImageDraw.Draw(im)
        for r in polyrois:
            #[ [ (x,y) ] ]
            roi =[]
            for x,y in r:
                roi.append( (int(x+.5),int(y+.5)))
            draw.polygon(roi,outline=0)

        for r in rectrois:
            # [ (x0,y0),(x1,y1) ]
            (x0,y0),(x1,y1) = r
            draw.rectangle(((x0,y0),(x1,y1)),outline=0)

        # now draw all cirlerois in reserved color
        for x,y,r in circlerois:
            # [ (x,y,r) ]
            draw.ellipse((x-r,y-r,x+r,y+r), outline=0)
        del draw

        # convert to RGB for JPG, cause JPG doesn't do PALETTE and PNG is much larger
        im = im.convert("RGB")

        imsi = im.size
        if max(imsi)>2048:
            ratio = 2048./max(imsi)
            im = im.resize( (int(imsi[0]*ratio+.5), int(imsi[1]*ratio+.5)),Image.ANTIALIAS)
        im.save(fname)

