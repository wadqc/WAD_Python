__author__ = 'aschilha'
import scipy.ndimage as scind
import numpy as np

def FiniteDifference1D(pSrc,BC="BC_MIRROR",order=0):

    pDest = np.zeros(pSrc.shape[0],dtype=float)
    if(order == 0):  # Gaussian just blurring, here nothing
        for i in range(0,pSrc.shape[0]):
            pDest[i] = pSrc[i]

        return pDest

    length = pSrc.shape[0]
    if(length <2):
        print "FiniteDifference1D: SKIP len<2:",length


    if(order == 1):
        for i in range(0,pSrc.shape[0]-1):
            pDest[i] = pSrc[i+1]-pSrc[i]

        if(BC == "BC_ZERO"):
            pDest[-1] = 0 - pSrc[-1]
        elif(BC == "BC_CONT"):
            pDest[-1] = 0.
        else:
            pDest[-1] = pSrc[-2]-pSrc[-1]

        return pDest

    print "FiniteDifference1D: ERROR: order",order,"not implemented yet!"
    return None

def linearInterExtrapolate(xarr,yarr,xpos):
    """
    linearly interpolates if xpos within bounds of xarr,
    else linearly extrapolates
    """
    result = 0.
    xref_id = -1
    num = len(xarr)
    if(xarr[0]<xarr[-1]): # xarr ordered from high to low
        if(xarr[0]>xpos): # needs extrapolation
            xref_id = 0
        if(xarr[-1]<=xpos): # needs extrapolation
            xref_id = num-2
        if(xref_id<0):
            for ix in range(0,num-1):
                if(xarr[ix]<=xpos and xarr[ix+1]>xpos):
                    xref_id = ix
                    break
    else: # xarr ordered form low to high
        xref_id = -1
        if(xarr[0]<xpos): # needs extrapolation
            xref_id = 0
        if(xarr[-1]>=xpos): # needs extrapolation
            xref_id = num-2
        if(xref_id<0):
            for ix in reversed(range(0,num-1)):
                if(xarr[ix]>xpos and xarr[ix+1]<=xpos):
                    xref_id = ix
                    break
    if(xref_id<0):
         print "***ERROR: Cannot linearly interpolate at ",xpos
         return result
    result = (xpos-xarr[xref_id])/(xarr[xref_id+1]-xarr[xref_id])*(yarr[xref_id+1]-yarr[xref_id])+yarr[xref_id]
    return result

def AreaUnderCurve(freqs,merits):
    area = 0.
    if(len(freqs) == 0 or len(merits) == 0):
        return area
    for k in range(0,len(freqs)-1):
        y1 = merits[k]
        y2 = merits[k+1]
        dx = freqs[k+1]-freqs[k]
        if(y2<y1):
            y2 = merits[k]
            y1 = merits[k+1]
        area += dx*(y1+(y2-y1)/2.)
    area /= freqs[-1]
    return area

def MTF10pct(freqs,mtf):
    if(len(freqs) == 0 or len(mtf) == 0):
        return 0
    return linearInterExtrapolate(mtf, freqs, 0.1)

