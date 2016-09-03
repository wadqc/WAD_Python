# wadwrapper_lib.py
from __future__ import print_function

"""
This provides a number of DICOM handling and math routines
which are often used in plugins developed at UMCU.
This includes:
 * reading a 2D image
 * reading of a folder of 2D slices as a 3D image (excluding non-image files)
 * reading an Enhanced DICOM object as a 3D image
 * scaling data values accoriding to slope and offset tags
 * accessing (nested) DICOM tags in 2D, 3D, and Enhanced DICOM

"""
# Changelog:
# 20160902: python3 compatible
# 20160826: for RGB US default to use only red channel
# 20150420: Added scipy connectedComponents
# 20150203: Added docstrings; added raise exceptions; clean-up of code

US_RGB_USE_RED = True # force using only the RED channel in RGB, or remove all data where there is a difference between R,G,B

import dicom
try:
    import pydicom_series as dcmseries
except ImportError:
    from pyWADLib import pydicom_series as dcmseries
    
try: # moved to errors in 'bleeding edge' pydicom
    import dicom.errors as dicomExceptions
except ImportError: 
    import dicom.filereader as dicomExceptions
    
import sys
import os
import numpy as np
import scipy.ndimage as scind

stMode2D = "2D"
stMode3D = "3D"
stModeEnhanced = "Enhanced"
stModeBTO = "BTO"

### DICOM
# Helper functions
def _readSingleImage(filename,headers_only):
    """
    Internal function to read a single file as a DICOM image, optionally skipping pixel data.
    """

    dcmInfile = None
    try:
        dcmInfile = dicom.read_file(filename, stop_before_pixels=headers_only) # Read only header!
    except dicomExceptions.InvalidDicomError:
        try:
            dcmInfile = dicom.read_file(filename,force=True, stop_before_pixels=headers_only)
        except:
            pass
    return dcmInfile

def getDICOMMode(dcmInfile):
    """
    Determine the structure of the DICOM image: 2D, 3D, Enhanced, BTO
    """
    dcmMode = stMode3D
    try:
        num = len(dcmInfile._datasets)
    except AttributeError:
        dcmMode = stMode2D
        try:
            sopclass = dcmInfile.SOPClassUID
            if sopclass == 'Enhanced MR Image Storage': #'1.2.840.10008.5.1.4.1.1.66'
                dcmMode = stModeEnhanced
            elif sopclass == 'Breast Tomosynthesis Image Storage' or sopclass == '1.2.840.10008.5.1.4.1.1.13.1.3':
                dcmMode = stModeBTO

        except AttributeError:
            pass
    return dcmMode

def testIfEnhancedDICOM(filename):
    """
    Test if DICOM file to be read is an EnhancedDICOM object (but not BTO)
    """
    dcmInfile = _readSingleImage(filename,headers_only=True)
    if dcmInfile is None:
        # Not a dicom file
        return False

    #Skip the file if it is of SOPClassUID 1.2.840.10008.5.1.4.1.1.66
    try:
        sopclass = dcmInfile.SOPClassUID
    except AttributeError:
        # some other kind of dicom file and no reason to exclude it
        return False
    if sopclass == 'Enhanced MR Image Storage': #'1.2.840.10008.5.1.4.1.1.66'
        print('Image is Enhanced MR')
        return True
    else:
        if sopclass == 'Breast Tomosynthesis Image Storage' or sopclass == '1.2.840.10008.5.1.4.1.1.13.1.3':
            print('Image is BTO')
        return False

def _readGroupSequenceTag(sequencetags,keytag):
    """
    Return value of tag in a sequence of tags of a EnhancedDICOM header object.
    Sequence should be either a SharedFunctionalGroupsSequence or a PerFrameFunctionalGroupsSequence
    """
    value = ''
    for seq in sequencetags.keys():
        subsetA = sequencetags[(seq.group,seq.elem)]
        if subsetA.VR != 'SQ' or subsetA.VM == 0:
            continue
        subsetA = sequencetags[(seq.group,seq.elem)][0]
        for elemA in subsetA.keys():
            if (elemA.group,elemA.elem) == (keytag.group, keytag.elem):
                value = subsetA[(elemA.group,elemA.elem)].value
    return value

def readDICOMtag(key,dcmInfile,imslice=0): # slice=2 is image 3
    """
    Returns the value of the dicomtag 'key' for image slice 'imslice'
    imslice is ignored for 2D images.
    A key is defined as e.g. "0018,11A0", but could also be a nested tag "0018,11A0,0008,1001"
    Returns empty string if tag is not defined in the dicomfile.
    """
    
    # transcribe "0018,11A0" as a tag entry
    keytag = dicom.tag.Tag(key.split(',')[0],key.split(',')[1])
    
    if len(key.split(','))>2: # key has more than 2 groups of digits, so presumable a nested tag
        return readNestedDICOMtag(key,dcmInfile,imslice)
    try:
        dicomMode = getDICOMMode(dcmInfile)
        if dicomMode == stMode2D:
            value =  dcmInfile[keytag].value
        elif dicomMode == stMode3D:
            value =  dcmInfile._datasets[imslice][keytag].value
        else: # EnhancedDicom or BTO
            try:
                value =  dcmInfile[keytag].value
            except:
                value = ''
            
            if value == '': # not found, look for it in sharedtags
                value = _readGroupSequenceTag(dcmInfile.SharedFunctionalGroupsSequence[0],keytag)

            if value == '': # not found, look for it in perframetags
                value = _readGroupSequenceTag(dcmInfile.PerFrameFunctionalGroupsSequence[imslice],keytag)

            if value == '': # not found, look for it in perframetags
                print("[ww_readDICOMtag] key",key,"not found")

    except:
        print("[ww_readDICOMtag] Exception for key",key)
        value = ""
    return value

def readNestedDICOMtag(key,dcmInfile,imslice): # slice=2 is image 3
    """
    Returns the value of the nested dicomtag 'key' for image slice 'imslice'
    imslice is ignored for 2D images.
    A nested tag is defined as e.g "0018,11A0,0008,1001" meaning tag "0008,1001" in subset tag "0018,11A0"
    Returns empty string if tag is not defined in the dicomfile.
    """
    # transcribe "0018,11A0" as a tag entry
    lv0 = dicom.tag.Tag(key.split(',')[0],key.split(',')[1])
    lv1 = dicom.tag.Tag(key.split(',')[2],key.split(',')[3])

    dicomMode = getDICOMMode(dcmInfile)
    try:
        if dicomMode == stMode2D:
            value =  dcmInfile[lv0][0][lv1].value
        elif dicomMode == stMode3D:
            value =  dcmInfile._datasets[imslice][lv0][0][lv1].value
        else: # EnhancedDicom or BTO
            #sharedTags = dcmInfile.SharedFunctionalGroupsSequence[0]
            try:
                value =  dcmInfile[lv0][0][lv1].value
            except:
                value = ''
            if value == '':
                sharedTags = dcmInfile.SharedFunctionalGroupsSequence[0]
                for seq in sharedTags.dir():
                    subsetA = sharedTags.data_element(seq)[0]
                    for elemA in subsetA.keys():
                        if (elemA.group,elemA.elem) == (lv0.group, lv0.elem):
                            subsetB = subsetA[(elemA.group,elemA.elem)]
                            if subsetB.VM>0 and subsetB.VR == 'SQ':
                                for elemB in subsetB[0].keys():
                                    if (elemB.group,elemB.elem) == (lv1.group, lv1.elem):
                                        value = elemB.value


    except:
        print("[ww_readNestedDICOMtag] Exception for key",key,dicomMode)
        value = ""
    return value

def removeBogusDICOMfiles(instancedict):
    """
    When a folder of files is offered, prepare a list of files that are valid DICOM image files;
    this removes some report files that would break reading a folder of slices as a 3D image
    """
    results = []
    for fname in instancedict:
        skipMe = False
        dcmInfile = _readSingleImage(fname,headers_only=True)
        if dcmInfile is None:
            # Not a dicom file
            continue

        #AS: Skip the file if it is of SOPClassUID 1.2.840.10008.5.1.4.1.1.66
        try:
            sopclass = dcmInfile.SOPClassUID
        except AttributeError:
            # some other kind of dicom file and no reason to exclude it
            results.append(fname)
            continue
        if sopclass == 'Raw Data Storage': #'1.2.840.10008.5.1.4.1.1.66'
            print('Skipping RAW file %s' %fname)
            continue
        else:
            results.append(fname)

    return results

def prepareEnhancedInput(filename,headers_only,logTag="[prepareEnhancedInput] "):
    """
    Reads filename as an EnhancedDICOM object. If not headers_only , scaling the pixel values according to RescaleIntercept and RescaleIntercept of each frame
    and transposing for pyqtgraph format.
    Raises ValueError if file cannot be opened as a DICOM object
    Returns raw dcmfile, scaled and transposed pixeldata (or None), type of DICOM object.
    """
    dcmInfile = None
    pixeldataIn = None

    # Repair Mapping for MR; not correct for all scanner
    keymapping = [
                     ("0040,9096,0040,9224","0028,1052"), # Real World Value Intercept -> Rescale Intercept
                     ("0040,9096,0040,9225","0028,1053"), # Real World Value Slope -> Rescale Slope
    ]

    dcmInfile = _readSingleImage(filename,headers_only=headers_only)
    if dcmInfile is None:
        # Not a dicom file
        raise ValueError("{} ERROR! {} is not a valid Enhanced DICOM object".format(logTag, filename))

    if not headers_only:
        pixeldataIn = dcmInfile.pixel_array.astype(int)
        perframe = dcmInfile.PerFrameFunctionalGroupsSequence # Sequence of perframeTags
        for i in range(len(perframe)):
            intercept = perframe[i].PixelValueTransformationSequence[0].RescaleIntercept
            slope = perframe[i].PixelValueTransformationSequence[0].RescaleSlope
            pixeldataIn[i] = intercept + slope*pixeldataIn[i]
        pixeldataIn = np.transpose(pixeldataIn,(0,2,1))
    return dcmInfile,pixeldataIn,getDICOMMode(dcmInfile)

def prepareInput(instancedict,headers_only,logTag="[prepareInput] "):
    """
    Reads inputfile as an EnhancedDICOM object; if an EnhancedDICOM object is detected, prepareEnhancedInput is called.
    Checks if the input is as expected: number of slices etc.
    If not headers_only, scaling the pixel values according to RescaleIntercept and RescaleIntercept of each frame
    and transposing for pyqtgraph format.
    Raises ValueError if file cannot be opened as a DICOM object
    Returns raw dcmfile, scaled and transposed pixeldata (or None), type of DICOM object.
    """
    # compile a list of valid files to read
    instancedict = removeBogusDICOMfiles(instancedict)

    dcmInfile = None
    pixeldataIn = None

    # Repair Mapping for MR; not correct for all scanner
    keymapping = [
                     ("0040,9096,0040,9224","0028,1052"), # Real World Value Intercept -> Rescale Intercept
                     ("0040,9096,0040,9225","0028,1053"), # Real World Value Slope -> Rescale Slope
    ]

    # Check if input data is as expected: MR needs 3D data
    if len(instancedict) == 1:
        ModeEnhanced = False
        ModeEnhanced = testIfEnhancedDICOM(instancedict[0])
        if ModeEnhanced:
            return prepareEnhancedInput(instancedict[0],headers_only=headers_only) # scaled and transposed

        if dcmInfile is None:
            filename = instancedict[0]
            dcmInfile = _readSingleImage(filename,headers_only=headers_only)
            if dcmInfile is None:
                # Not a dicom file
                raise ValueError("{} ERROR! {} is not a valid non-Enhanced DICOM object".format(logTag, filename))

            modality = dcmInfile.Modality
            
            if not headers_only:
                # need scaling for single slice, already done for series in pydicom_series, but not for some MR files
                if modality == 'CT' or modality == 'MR':
                    if not "RescaleIntercept" in dcmInfile: # in wrong place define for some MR files
                        dcmInfile.RescaleIntercept = readDICOMtag(keymapping[0][0],dcmInfile,0)
                        dcmInfile.RescaleSlope = readDICOMtag(keymapping[1][0],dcmInfile,0)

                    pixeldataIn = np.int16(np.transpose(dcmInfile.pixel_array.astype(int),(1,0)))
                    slope = dcmInfile.RescaleSlope
                    intercept = dcmInfile.RescaleIntercept
                    pixeldataIn = intercept + slope*pixeldataIn
                elif modality == 'MG' and getDICOMMode(dcmInfile) == stModeBTO:
                    print('!WARNING! MG BTO dataset! DICOM info is NOT properly adjusted, no scaling applied yet!')
                    pixeldataIn = np.transpose(dcmInfile.pixel_array,(0,2,1))
                elif modality == 'MG' or modality == 'CR' or modality == 'DX':
                    pixeldataIn = dcmInfile.pixel_array.transpose()
                elif modality == 'RF': # fixme! 2D and 3D
                    pixeldataIn = dcmInfile.pixel_array.transpose()
                elif modality == 'US':
                    rgbmode = (dcmInfile.SamplesPerPixel == 3)
                    if not rgbmode:
                        if len(np.shape(dcmInfile.pixel_array)) == 2:
                            pixeldataIn = dcmInfile.pixel_array.transpose()
                        elif len(np.shape(dcmInfile.pixel_array)):
                            pixeldataIn = np.transpose(dcmInfile.pixel_array,(0,2,1))
                            pixeldataIn = pixeldataIn[0]
                    else:
                        # AS: this fix is needed in pydicom < 1.0; maybe solved in later versions?
                        try:
                            nofframes = dcmInfile.NumberOfFrames
                        except AttributeError:
                            nofframes = 1
                        if dcmInfile.PlanarConfiguration==0:
                            pixel_array = dcmInfile.pixel_array.reshape(nofframes, dcmInfile.Rows, dcmInfile.Columns, dcmInfile.SamplesPerPixel)
                        else:
                            pixel_array = dcmInfile.pixel_array.reshape(dcmInfile.SamplesPerPixel, nofframes, dcmInfile.Rows, dcmInfile.Columns)

                            # force using only the RED channel in RGB.
                            if US_RGB_USE_RED == True:
                                if len(np.shape(pixel_array)) == 3: #2d rgb
                                    pixeldataIn = pixel_array[:,:,0].transpose()
                                elif len(np.shape(pixel_array)) == 4: #3d rgb
                                    pixeldataIn = (pixel_array[-1,:,:,0]).transpose()
                            else:
                                # remove all data where there is a difference between R,G,B
                                if len(np.shape(pixel_array)) == 3: #2d rgb
                                    pixeldataIn = pixel_array[:,:,0].transpose()
                                    pixeldataInR = (pixel_array[:,:,0]).transpose()
                                    pixeldataInG = (pixel_array[:,:,1]).transpose()
                                    pixeldataInB = (pixel_array[:,:,2]).transpose()
                                elif len(np.shape(pixel_array)) == 4: #3d rgb
                                    pixeldataIn = (pixel_array[-1,:,:,0]).transpose()
                                    pixeldataInR = (pixel_array[-1,:,:,0]).transpose()
                                    pixeldataInG = (pixel_array[-1,:,:,1]).transpose()
                                    pixeldataInB = (pixel_array[-1,:,:,2]).transpose()
                                # remove rgb info
                                for y in range(dcmInfile.Rows):
                                    for x in range(dcmInfile.Columns):
                                        r = pixeldataInR[x,y]
                                        g = pixeldataInG[x,y]
                                        b = pixeldataInB[x,y]
                                        ma = max(r,g,b)
                                        mi = min(r,g,b)
                                        if ma != mi:
                                            pixeldataIn[x,y] = 0

    else:
        path = os.path.dirname(instancedict[0])
        for ip in instancedict:
            if os.path.dirname(ip) != path:
                raise ValueError("{} ERROR! multiple would-be dicom files scattered over multiple dirs!".format(logTag))

        try:
            dcmInfile = dcmseries.read_files(path, True, readPixelData=False)[0]
            if not headers_only: # NOTE: Rescaling is already done pydicom_series, but maybe not for stupid MR
                pixeldataIn = np.transpose(dcmInfile.get_pixel_array(),(0,2,1))

        except Exception as e:
            raise ValueError("{} ERROR! {} is not a valid non-Enhanced DICOM series".format(logTag, path))

    return dcmInfile,pixeldataIn,getDICOMMode(dcmInfile)


### math
class connectedComponents():
    maskIn = None
    cca    = None
    nb_labels = None
    cluster_sizes = []
    
    def __init__(self):
        self.maskIn = None
        self.cca = None
        self.nb_labels = 0
        self.cluster_sizes = []

    def run(self,pBool):
        self.maskIn = pBool
        self.cca,self.nb_labels = scind.label(self.maskIn)
        self.cluster_sizes = []
        return self.cca,self.nb_labels
    
    def removeSmallClusters(self,minsize):
        mask_size = self.clusterSizes()<minsize
        self.cca[mask_size[self.cca]] = 0

        labels = np.unique(self.cca)
        self.nb_labels = len(labels)
        self.cca = np.searchsorted(labels, self.cca)
        self.cluster_sizes = []
        
    def clusterSizes(self):
        if len(self.cluster_sizes) == 0:
            self.cluster_sizes = scind.sum(self.maskIn, self.cca, range(self.nb_labels + 1))
        return self.cluster_sizes

    def indicesOfCluster(self,val):
        clus = np.where(self.cca == val)
        clus = [(x,y) for x,y in zip(clus[0],clus[1])]
        return clus

### peak finding
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
    #return array(maxtab), array(mintab) 3 converts indices to float!
    return maxtab, mintab

### thresholding
def threshold_adaptive(image, block_size, method='gaussian', offset=0,
                       mode='reflect', param=None):
    """
    from skitimage 0.8
    Applies an adaptive threshold to an array.

    Also known as local or dynamic thresholding where the threshold value is
    the weighted mean for the local neighborhood of a pixel subtracted by a
    constant. Alternatively the threshold can be determined dynamically by a a
    given function using the 'generic' method.

    Parameters
    ==========
    image : (N, M) ndarray
        Input image.
    block_size : int
        Uneven size of pixel neighborhood which is used to calculate the
        threshold value (e.g. 3, 5, 7, ..., 21, ...).
    method : {'generic', 'gaussian', 'mean', 'median'}, optional
        Method used to determine adaptive threshold for local neighbourhood in
        weighted mean image.

        * 'generic': use custom function (see `param` parameter)
        * 'gaussian': apply gaussian filter (see `param` parameter for custom\
                      sigma value)
        * 'mean': apply arithmetic mean filter
        * 'median': apply median rank filter

        By default the 'gaussian' method is used.
    offset : float, optional
        Constant subtracted from weighted mean of neighborhood to calculate
        the local threshold value. Default offset is 0.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
        Default is 'reflect'.
    param : {int, function}, optional
        Either specify sigma for 'gaussian' method or function object for
        'generic' method. This functions takes the flat array of local
        neighbourhood as a single argument and returns the calculated
        threshold for the centre pixel.

    Returns
    =======
    threshold : (N, M) ndarray
        Thresholded binary image

    References
    ==========
    .. [1] http://docs.opencv.org/modules/imgproc/doc/miscellaneous_transformations.html?highlight=threshold#adaptivethreshold

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> binary_image1 = threshold_adaptive(image, 15, 'mean')
    >>> func = lambda arr: arr.mean()
    >>> binary_image2 = threshold_adaptive(image, 15, 'generic', param=func)
    """
    thresh_image = np.zeros(image.shape, 'double')
    if method == 'generic':
        scind.generic_filter(image, param, block_size,
                             output=thresh_image, mode=mode)
    elif method == 'gaussian':
        if param is None:
            # automatically determine sigma which covers > 99% of distribution
            sigma = (block_size - 1) / 6.0
        else:
            sigma = param
        scind.gaussian_filter(image, sigma, output=thresh_image,
                              mode=mode)
    elif method == 'mean':
        mask = 1. / block_size * np.ones((block_size,))
        # separation of filters to speedup convolution
        scind.convolve1d(image, mask, axis=0, output=thresh_image,
                         mode=mode)
        scind.convolve1d(thresh_image, mask, axis=1,
                         output=thresh_image, mode=mode)
    elif method == 'median':
        scind.median_filter(image, block_size, output=thresh_image,
                            mode=mode)

    return image > (thresh_image - offset)

def __IJIsoData(data):
    """
    This is the original ImageJ IsoData implementation
    """
    count0 = data[0]
    data[0] = 0  # set to zero so erased areas aren't included
    countMax = data[-1]
    data[-1] = 0

    maxValue = len(data)- 1
    amin = 0
    while data[amin]==0 and amin<maxValue:
        amin += 1
    amax = maxValue
    while data[amax]==0 and amax>0:
        amax -= 1
    if amin>=amax:
        data[0]        = count0
        data[maxValue] = countMax;
        level = len(data)/2
        return level

    movingIndex = amin
    cond = True
    while cond:
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0
        for i in range(amin,movingIndex+1):
            sum1 += i*data[i]
            sum2 += data[i]

        for i in range(movingIndex+1,amax+1):
            sum3 += i*data[i]
            sum4 += data[i]

        result = (sum1/sum2 + sum3/sum4)/2.0
        movingIndex += 1
        cond = ((movingIndex+1)<=result and movingIndex<amax-1)

    data[0]        = count0
    data[maxValue] = countMax;
    level = int(round(result))
    return level
    
def threshold_isodata2(data):
    """
    Ripped from ImageJ "defaultIsoData"
    """
    maxCount = np.max(data)
    mode = np.argmax(data)

    data2 = np.copy(data)
    maxCount2 = 0
    for i,v in enumerate(data2):
        if v>maxCount2 and i != mode:
            maxCount2 = v

    if maxCount>maxCount2*2 and maxCount2 != 0:
        data2[mode] = int(maxCount2 * 1.5)

    return __IJIsoData(data2)
    
###
def extract(Z, shape, position, fill=np.NaN):
    """ Extract a sub-array from Z using given shape and centered on position.
      If some part of the sub-array is out of Z bounds, result is padded
         with fill value.

         **Parameters**
             `Z` : array_like
                Input array.

            `shape` : tuple
                Shape of the output array

            `position` : tuple
                Position within Z

            `fill` : scalar
                Fill value

         **Returns**
             `out` : array_like
                 Z slice with given shape and center

         **Examples**

         >>> Z = numpy.arange(0,16).reshape((4,4))
         >>> extract(Z, shape=(3,3), position=(0,0))
         [[ NaN  NaN  NaN]
          [ NaN   0.   1.]
          [ NaN   4.   5.]]

         Schema::

         
              +-----------+
              | 0   0   0 |  =  extract (Z, shape=(3,3), position=(0,0))
              |   +-------|---------+  
              | 0 | 0   1 | 2   3   |  =  Z
              |   |       |         |
              | 0 | 4   5 | 6   7   |
              +-----------+         |
                  | 8    9   10  11 |
                  |                 |
                  | 12   13  14  15 |
                  +-----------------+
             

         >>> Z = numpy.arange(0,16).reshape((4,4))
         >>> extract(Z, shape=(3,3), position=(3,3))
         [[ 10.  11.  NaN]
          [ 14.  15.  NaN]
          [ NaN  NaN  NaN]]

         Schema::

             +---------------+
             | 0   1   2   3 | = Z
             |               |
             | 4   5   6   7 |
             |       +-----------+
             | 8   9 |10  11 | 0 | = extract (Z, shape=(3,3),position=(3,3))
             |       |       |   |
             | 12 13 |14  15 | 0 |
             +---------------+   |
                     | 0   0   0 |
                     +-----------+
    """
    #    assert(len(position) == len(Z.shape))
    #    if len(shape) < len(Z.shape):
    #        shape = shape + Z.shape[len(Z.shape)-len(shape):]

    R = np.ones(shape, dtype=Z.dtype)*fill
    P  = np.array(list(position)).astype(int)
    Rs = np.array(list(R.shape)).astype(int)
    Zs = np.array(list(Z.shape)).astype(int)

    R_start = np.zeros((len(shape),)).astype(int)
    R_stop  = np.array(list(shape)).astype(int)
    Z_start = (P-Rs//2)
    Z_stop  = (P+Rs//2)+Rs%2

    R_start = (R_start - np.minimum(Z_start,0)).tolist()
    Z_start = (np.maximum(Z_start,0)).tolist()
    R_stop = (R_stop - np.maximum(Z_stop-Zs,0)).tolist()
    Z_stop = (np.minimum(Z_stop,Zs)).tolist()

    r = [slice(start,stop) for start,stop in zip(R_start,R_stop)]
    z = [slice(start,stop) for start,stop in zip(Z_start,Z_stop)]

    R[r] = Z[z]

    return R

