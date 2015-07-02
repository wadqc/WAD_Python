"""
Changelog:
    20150121: Code cleanup + modification to make module compatible with Philips camera's (TdW)
    [history here]
    ...
"""

import scipy
from scipy import ndimage
from scipy import optimize
from skimage.measure import block_reduce
import dicom
import pylab as pl
import numpy as np 
from numpy import ma

import sys, getopt, os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

from pyWAD.plugindata import PluginData
from pyWAD.pluginresults import PluginResults

import itertools
import xml.etree.ElementTree as ET

import operator
from copy import deepcopy

"""NEMA Definitions:

UFOV: the area of the detector that is used for imaging gamma rays and x-rays.
It is defined by a dimensioned figure supplied by the manufacturer.

CVOF: The area defined by scaling all linear dimensions of the UFOV by a factor
of 75%

Differential uniformity: The amount of count density change per defined unit distance when the detector's incident gamma radiation is a homogeneous flux over the field of measurement

DU_x,y = = 100% * (MAX - MIN)/(MAX + MIN), where the largest difference between any two pixels within a set of 5 contiguous pixels in a row or column shall be 
calculated. The calculation shall be done for X and Y directions independently and expressed as:

Integral uniformity: A measure of the maximum count density variation over a defined large area of the scintillation detector for a uniform input gamma flux to the UFOV of the camera

IU = 100% * (MAX - MIN)/(MAX + MIN)

"""

def flux(A, R_0, center_x, center_y):
    ''' Function template for fitting 2D dome.
    '''
    A        = float(A)
    R_0      = float(R_0)
    center_x = float(center_x)
    center_y = float(center_y)
    return lambda x,y: A*R_0 /( R_0**2 + (center_x-x)**2 + (center_y-y)**2 )**(3/2.)


def fitflux(data):
    """ Function for fitting a theoretical dome-function to the measured (uncorrected) data.
    """
    # parameter initialization
    params = np.array([1, 1, data.shape[0]/2, data.shape[0]/2])
    # (data>0) : prevents fitting of the 0-values outside the crystal
    #errorfunction = lambda p: ravel( (flux(*p)(*indices(data.shape)) - data)*(data>0) )
    # strip outer 10 pixels for optimal fit (prevent issues with hot crystal corners)
    errorfunction = lambda p: np.ravel( (flux(*p)(*np.indices(data.shape)) - data)*ndimage.binary_erosion(data>0,iterations=10) )
    p, success = optimize.leastsq(errorfunction, params)
    return p
 

def dome_correction(imagearray):
   """ Function for performing dome-correction.
         Input  = uncorrected flood data
         Output = dome-corrected flood data
   """
   data = imagearray
   
   # fit function: flux() to the measured data
   params = fitflux(data)
   print 'A, R_0, center_x, center_y : ', params
   fitdata = flux(*params)(*np.indices(data.shape))

   # difference between data and fit (inspection of systematic errors)
   #diff=(data-fitdata)*(data>0)
   #diff.tofile('diff_%d.raw' % i)

   # divide out the geometry
   fitdata = data/fitdata
   
   # preserve the number of counts
   fitdata = data.sum()*fitdata/fitdata.sum()
   
   # write back the corrected data to the dicom object and convert float64->uint16
   data = np.rint(fitdata).astype(np.uint16) 

   return data


def set_threshold(imagearray):
     return np.mean(imagearray)

def testmask():
     tmp = np.zeros((100,100))
     tmp[20:40,30:50] = 1
     return ma.make_mask(tmp)


def unifcalc(vector):
    """ Function calculates the contrast between minimum and maximum pixelvalue in "vector".
        Returns fraction not percentage!
    """
    _min = float(np.ma.min(vector))
    _max = float(np.ma.max(vector))
    return (_max - _min) / (_max + _min)


def nema_smooth(inarray):
    """ Function smoothes inarray with a nine-point NEMA kernel.
        The total number of counts in "inarray" is preserved.
    """
    k = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]])

    outarray = ndimage.convolve(inarray, k)/k.sum()
    
    return outarray


def diff_data(inarray):
    """ Function calculates the Differential Uniformity in a
        horizontal and vertical range of five pixels with maximum contrast.

        Returns contrast percentages and upper-left pixel of the locations
        with highest contrast.
    """
    outputx = []
    outputy = []

    # convenient method to create pairs of 5 pixels in vert and horiz direction
    # by shifting the array 5 times in vert/horiz direction, and concatenating
    # the flattened shifted arrays. In order to correctly identify the first
    # pixel of the group of five pixels with the highest nonuniformity, the range
    # is running from -4..0
    for i in range(-4,1):
        tmpx = np.roll(inarray,i,0)
        tmpy = np.roll(inarray,i,1)
    
        outputx.append(tmpx.ravel())
        outputy.append(tmpy.ravel())


    difflistx = []
    difflisty = []

    # calculate the uniformity inside all possible pairs of five concurrent pixels
    # in vertical direction:
    for vector in np.ma.array(outputx).T:
        # ignore groups of five pixels (partly) outside the mask
        if np.ma.count_masked(vector) == 0:
            difflistx.append(100*unifcalc(vector))
            #difflistx.append(100*(ma.max(vector) - ma.min(vector))/(ma.max(vector) + ma.min(vector)))
        else:
            difflistx.append(0.0)

    img = np.reshape(np.array(difflistx),np.shape(inarray))
    vertmax=ndimage.maximum_position(img)

    # calculate the uniformity inside all possible pairs of five concurrent pixels
    # in horizontal direction:
    for vector in np.ma.array(outputy).T:
        # ignore groups of five pixels (partly) outside the mask
        if np.ma.count_masked(vector) == 0:
            difflisty.append(100*unifcalc(vector))
        else:
            difflisty.append(0.0)

    img = np.reshape(np.array(difflisty),np.shape(inarray))
    hormax=ndimage.maximum_position(img)

    highx = np.max(difflistx)
    # lowest DU in vert direction, larger than zero
    lowx =  np.min(np.array(difflistx)[np.nonzero(difflistx)])

    highy = np.max(difflisty)
    # lowest DU in horiz direction, larger than zero
    lowy =  np.min(np.array(difflisty)[np.nonzero(difflisty)])

    #return highx, highy, difflistx,difflisty
    return highx, highy, vertmax, hormax


def save_imgmap(inarray,vertmax,hormax,filename):
    """ Converts input-array to image, with superposed the 5-pixel row/column with
       highest DU, the minimal and maximum pixel coordinates.

        input: inarray  = ufov or cfov
        vertmax  = first coordinate of 5-pixel column (vertical)
        hormax   = first coordinate of 5-pixel row (horizontal)
        filename = png-filename of output
    """
    wi,he = np.shape(inarray)
    rgb = np.zeros((wi, he, 3), dtype=np.uint8)
    _max=np.max(inarray)
    _min=np.min(inarray)*0.95
    grayvalue = np.round(255/(_max-_min)*(ma.filled(inarray,fill_value=0)-_min))
    grayvalue[grayvalue<0]=0
    rgb[:,:, 0] = grayvalue
    rgb[:,:, 1] = grayvalue
    rgb[:,:, 2] = grayvalue
    
    rgb[vertmax[0]:vertmax[0]+5,vertmax[1],:] = (255,150,50)   # orange
    rgb[hormax[0],hormax[1]:hormax[1]+5,:] = (0,200,0)         # green
    minpos=ndimage.minimum_position(inarray)
    maxpos=ndimage.maximum_position(inarray) 
    rgb[minpos[0], minpos[1], :] = (0,0,255)                   # blue
    rgb[maxpos[0], maxpos[1], :] = (255,0,0)                   # red    

    rgb = np.zeros((wi, he, 3), dtype=np.uint8)
    _max=np.max(inarray)
    _min=np.min(inarray)*0.95
    grayvalue = np.round(255/(_max-_min)*(ma.filled(inarray,fill_value=0)-_min))
    grayvalue[grayvalue<0]=0
    rgb[:,:, 0] = grayvalue
    rgb[:,:, 1] = grayvalue
    rgb[:,:, 2] = grayvalue
    rgb = rgb.repeat(16, axis=0).repeat(16, axis=1)
    
    # d=line thickness of roi
    d=4    
    
    # 5 pixels in vertical direction, starting from 16*(vertmax[0],vertmax[1]) 
    # left edge
    rgb[16*vertmax[0]:16*(vertmax[0]+5),16*vertmax[1]:16*vertmax[1]+d,:] = (255,150,50)           # orange
    # right edge
    rgb[16*vertmax[0]:16*(vertmax[0]+5),16*(vertmax[1]+1)-d:16*(vertmax[1]+1),:] = (255,150,50)   # orange
    # top edge
    rgb[16*vertmax[0]:16*vertmax[0]+d,16*vertmax[1]:16*(vertmax[1]+1)-1,:] = (255,150,50)         # orange
    # bottom edge
    rgb[16*(vertmax[0]+5)-d:16*(vertmax[0]+5),16*vertmax[1]:16*(vertmax[1]+1)-1,:] = (255,150,50) # orange
 
    # 5 pixels in horizontal direction, starting from 16*(vertmax[0],vertmax[1])
    # left edge
    rgb[16*hormax[0]:16*(hormax[0]+1),16*hormax[1]:16*hormax[1]+d,:] = (0,200,0)                  # green
    # right edge
    rgb[16*hormax[0]:16*(hormax[0]+1),16*(hormax[1]+5)-d:16*(hormax[1]+5),:] = (0,200,0)          # green
    # top edge
    rgb[16*hormax[0]:16*hormax[0]+d,16*hormax[1]:16*(hormax[1]+5)-1,:] = (0,200,0)                # green
    # bottom edge
    rgb[16*(hormax[0]+1)-d:16*(hormax[0]+1),16*hormax[1]:16*(hormax[1]+5)-1,:] = (0,200,0)        # green    
 
    minpos=ndimage.minimum_position(inarray)
    maxpos=ndimage.maximum_position(inarray)
    
    # position of lowest pixel value
     # left edge
    rgb[16*minpos[0]:16*(minpos[0]+1),16*minpos[1]:16*minpos[1]+d,:] = (0,0,255)                  # blue
    # right edge
    rgb[16*minpos[0]:16*(minpos[0]+1),16*(minpos[1]+1)-d:16*(minpos[1]+1),:] = (0,0,255)          # blue
    # top edge
    rgb[16*minpos[0]:16*minpos[0]+d,16*minpos[1]:16*(minpos[1]+1)-1,:] = (0,0,255)                # blue
    # bottom edge
    rgb[16*(minpos[0]+1)-d:16*(minpos[0]+1),16*minpos[1]:16*(minpos[1]+1)-1,:] = (0,0,255)        # blue
    
    # position of highest pixel value    
    # left edge
    rgb[16*maxpos[0]:16*(maxpos[0]+1),16*maxpos[1]:16*maxpos[1]+d,:] = (255,0,0)                  # red
    # right edge
    rgb[16*maxpos[0]:16*(maxpos[0]+1),16*(maxpos[1]+1)-d:16*(maxpos[1]+1),:] = (255,0,0)          # red 
    # top edge
    rgb[16*maxpos[0]:16*maxpos[0]+d,16*maxpos[1]:16*(maxpos[1]+1)-1,:] = (255,0,0)                # red
    # bottom edge
    rgb[16*(maxpos[0]+1)-d:16*(maxpos[0]+1),16*maxpos[1]:16*(maxpos[1]+1)-1,:] = (255,0,0)        # red
    
    #rgb[16*minpos[0]:16*(minpos[0]+1), 16*minpos[1]:16*(minpos[1]+1), :] = (0,0,255)     # blue
    #rgb[16*maxpos[0]:16*(maxpos[0]+1), 16*maxpos[1]:16*(maxpos[1]+1), :] = (255,0,0)     # red 
    
    '''   
    rgb[vertmax[0]:vertmax[0]+5,vertmax[1],:] = (255,150,50)   # orange
    rgb[hormax[0],hormax[1]:hormax[1]+5,:] = (0,200,0)         # green
    minpos=ndimage.minimum_position(inarray)
    maxpos=ndimage.maximum_position(inarray) 
    rgb[minpos[0], minpos[1], :] = (0,0,255)                   # blue
    rgb[maxpos[0], maxpos[1], :] = (255,0,0)                   # red    
    '''
    
    #imshow(rgb,interpolation='None')
    
    # truncate image
    # UL, LR = bounding_box(rgb[:,:,0])
    # rgb = rgb[UL[0]:LR[0],UL[1]:LR[1]]

    pl.imsave(filename,ma.filled(rgb,fill_value=0))



def get_dims(squaremask):
    """ return number of unmasked pixels along horizontal and vertical profiles
        through the center of the image (=size of mask)
    """

    dimx,dimy = np.shape(squaremask)
    halfpointval = squaremask[dimx/2,dimy/2]
    
    horvec = squaremask[dimx/2,:]
    vervec = squaremask[:,dimy/2]

    return ma.count(horvec), ma.count(vervec)
        

# alternative for block_reduce (currently not used)
def downsample(arr, factor):
    assert isinstance(factor, int), type(factor)
    sx, sy = arr.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy/factor * (X/factor) + Y/factor
    result = ndimage.sum(arr, labels=regions, index=np.arange(regions.max() + 1))
    result.shape = (sx/factor, sy/factor)
    return result 


def bounding_box(inarray):
    """ Function returns the upper-left and lower-right coordinates of the bounding
        box enclosing the imagedata in "inarray". the UL and LR coordinates are
        inside the FOV.
    """
    temp=np.argwhere(inarray)
    UL_x, UL_y = temp.min(0)
    LR_x, LR_y = temp.max(0)
    return np.array(([UL_x,UL_y], [LR_x,LR_y]))


def create_cfov(ufov):
    """ Function calculates the CFOV from a given UFOV mask.
    """

    # determine UL and LR of bounding box
    UL, LR = bounding_box(ufov)

    # convert discrete pixel coordinates (i,j) to continuous coordinate grid (fi,fj)
    # (i,j)=(0,0) corresponds to center of UL pixel
    # (fi,fj)=(0,0) corresponds to UL corner of UL pixel
    #   fi=i+0.5
    #   fj=j+0.5
    # convert UL and LR of FOV to continuous coordinate grid (along pixel boundaries)
    fUL = UL - 0.5
    fLR = LR + 0.5
    fCenter = (fUL+fLR)/2
    fUL_CFOV = 0.25*fCenter + 0.75*fUL
    fLR_CFOV = 0.25*fCenter + 0.75*fLR
    
    ''' NEMA: "Any pixel that has at least 50% of its area inside the CFOV shall
               be included within the CFOV analysis."
    '''
    UL_CFOV = np.ceil(fUL_CFOV-0.5).astype('int')
    LR_CFOV = np.floor(fLR_CFOV-0.5).astype('int')
    
    cfov = deepcopy(ufov)
    
    # initialize mask with Trues (=masked)
    cfov.mask = True
    
    # use UL and LR of calculated CFOV to adjust the CFOV-mask
    cfov.mask[UL_CFOV[0]:LR_CFOV[0]+1,UL_CFOV[1]:LR_CFOV[1]+1] = False
    
    return cfov




def nema_data_preprocess(imagearray,resamplesize):
    """ Function performs preprocessing on the input image:
          - resample to 64x64
          - smoothing with nine-point NEMA kernel
          - calculation of UFOV and CFOV regions
        Returns: masked UFOV and CFOV arrays
    """    

    print "data preprocessing: "
    print "array size:",  np.shape(imagearray), "resamplesize: ", resamplesize
    # causes Fourier artifacts
    #imagearray = resample(resample(imagearray,resamplesize[0],axis=0),resamplesize[1],axis=1)
    if resamplesize[0]>0 and resamplesize[1]>0:
      imagearray = block_reduce(imagearray, block_size=(np.shape(imagearray)[0]/resamplesize[0],np.shape(imagearray)[1]/resamplesize[1]),func=np.sum)
    imagearray = imagearray.astype('float64')

    imagearray = nema_smooth(imagearray) 


    """ NEMA step 1: "First, any pixels at the edge of UFOV containing less
                      than 75% of the mean counts per pixel in the CFOV shall
                      be set to zero."
    """                      
    # first estimate of UFOV (use segmentation-threshold = mean value of entire image)
    threshold = set_threshold(imagearray)
    ufov = ma.masked_less(imagearray,threshold,copy=False)
    
    # use NEMA guidelines to determine UFOV
    cfov = create_cfov(ufov)
    
    # average of CFOV
    cfov_average = set_threshold(cfov)
    ufov = ma.masked_less(imagearray,0.75*cfov_average,copy=False)
    
    
    """ NEMA step 2: "Second, those pixels which now have at least one of their
        four directly abutted neighbors containing zero counts, will be also
        set to zero. The remaining non-zero pixels are the pixels to be included
        in the analysis for the UFOV.
    """
    ufov.mask=scipy.ndimage.binary_dilation(ufov.mask,iterations=1)
    # based on final UFOV, create a new CFOV
    cfov = create_cfov(ufov)

    # FIXME: inconsistent use of xy (in python y is the horizontal dimension)
    #ux, uy = get_dims(ufov)

    ufov.fill_value=0
    cfov.fill_value=0

    return ufov, cfov



def calculate_nema_uniformity (imagearray, resamplesize, results, domecorrection=False):
    """ Wrapper function for flood calculation according to NEMA recommendations.
        Input:
          imagearray     : NxN numpy input array
          resamplesize   : downsample size (MxM), typically (64,64)
          results        : instance of PluginData-class (container for generated results)
          domecorrection : Perform dome correction? [True, False]

        Dome correction can be used for intrinsic uniformity measurements (e.g. with
        Siemens camera's) where the distance between point-source and detector is
        smaller than 5 times the maximum FOV dimension.
    """

    if domecorrection == True:
        print 'Performing dome-correction...'
        imagearray = dome_correction(imagearray)

    IUufov = 0
    IUcfov = 0
    DUxufov = 0
    DUyufov = 0 
    DUxcfov = 0
    DUycfov = 0
    
    imshape = np.shape(imagearray)
    
    try:
         ufov, cfov = nema_data_preprocess(imagearray,resamplesize)
    except:
         print "warning: could not preprocess ufov, cfov"
         ufov, cfov = np.ones((resamplesize))

    ufov.fill_value=0
    cfov.fill_value=0

    #unifcalc = lambda arr: 100*(ma.max(arr) - ma.min(arr))/(ma.max(arr) + ma.min(arr))
    unifxy_min = lambda arr: ndimage.minimum_position(arr) 
    unifxy_max = lambda arr: ndimage.maximum_position(arr) 

    IUufov = 100*unifcalc(ufov)
    IUufov_min = unifxy_min(ufov)
    IUufov_max = unifxy_max(ufov)
    IUcfov = 100*unifcalc(cfov)
    IUcfov_min = unifxy_min(cfov)
    IUcfov_max = unifxy_max(cfov) 


    DUxufov_val,DUyufov_val, DUxufov_coord, DUyufov_coord = diff_data(ufov)
    DUxcfov_val,DUycfov_val, DUxcfov_coord, DUycfov_coord = diff_data(cfov)

    output = DUxufov_val, DUyufov_val, DUxufov_coord, DUyufov_coord, DUxcfov_val, DUycfov_val, DUxcfov_coord, DUycfov_coord, IUufov, IUcfov, ufov, cfov

    return output
