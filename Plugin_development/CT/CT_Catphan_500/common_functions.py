'''
This module contains common functions used to analyze the catphan phantom.
Original code: https://pypi.python.org/pypi/pylinac
'''

import numpy as np
from numpy import sqrt
import os, os.path as osp
from scipy.ndimage import center_of_mass as cm
from skimage.feature import peak_local_max


def find_slices_of_interest(data,params):
    '''
    This routine expects the data to be ordered z-y-x such that a specific z-slice can be addressed using array[id,:,:]
    
    The workflow is as follows:
    1) walk through the stack of slices and for each slice:
     - determine if it is the positioning slice
     - determine if it is the center of a module 

    2) Given the scan direction and location of the positioning slice we identify the HU, SR and UN modules


    Parameters that can be set in the params block of the config-xml file:

    <scandirection> 
    <linewidth> 
    <profilelen>
    <detthreshold> 
    <roidim> 
    <refroidim>
    <refroioffset>

    '''

    scandirection=params.find('scandirection').text
    linewidth=int(params.find('linewidth').text)
    profilelen=int(params.find('profilelen').text)
    detthreshold=int(params.find('detthreshold').text)
    roidim=int(params.find('roidim').text)
    refroidim=int(params.find('refroidim').text)
    refroioffset=int(params.find('refroioffset').text)

    com = cm(np.sum(data,0))
    print ('com',com)
    number_of_slices = np.shape(data)[0]
    possiblerefslice = []
    modslice = []

    
    for sliceid in range(number_of_slices):
        tmpslice = data[sliceid,:,:]
        lineprofile_E = tmpslice[com[0]-linewidth:com[0]+linewidth,com[1]:com[1]+profilelen]
        lineprofile_W = tmpslice[com[0]-linewidth:com[0]+linewidth,com[1]-profilelen:com[1]]
        lineprofile_S = tmpslice[com[0]:com[0]+profilelen,com[1]-linewidth:com[1]+linewidth]
        lineprofile_N = tmpslice[com[0]-profilelen:com[0],com[1]-linewidth:com[1]+linewidth]
        
        lines = [lineprofile_N,lineprofile_S,lineprofile_E,lineprofile_W]
        


        if np.min([np.max(line) for line in lines]) > detthreshold:
            tmpcommean = np.mean(tmpslice[com[0]-roidim:com[0]+roidim,com[1]-roidim:com[1]+roidim]) 
            if (0 < tmpcommean) and (tmpcommean < 100):
                possiblerefslice.append(sliceid)

        modslice.append(np.max(tmpslice[com[0]-refroioffset-refroidim:com[0]-refroioffset+refroidim,com[1]-refroidim:com[1]+refroidim]))




    print modslice
    print possiblerefslice
    outdict = {'HU':None,'SR':None,'UN':None}

    if len(possiblerefslice)==1:
        refslice =  possiblerefslice[0]
        
        pkcrds = peak_local_max(np.array(modslice), min_distance=2)
 
        print ('pkcrds',pkcrds)
        idx = findindex(refslice,list(pkcrds))
        

        if scandirection == 'down':
            outdict['HU'] = int(pkcrds[idx][0])
            outdict['SR'] = int(pkcrds[idx-1][0])
            outdict['UN'] = int(pkcrds[idx-2][0])
        elif scandirection == 'up':
            outdict['HU'] = int(pkcrds[idx][0])
            outdict['SR'] = int(pkcrds[idx+1][0])
            outdict['UN'] = int(pkcrds[idx+2][0])

    return outdict

def findindex(number,mylist):
    if number in mylist:
        return mylist.index(number)
    elif number-1 in mylist:
        return mylist.index(number-1)
    elif number+1 in mylist:
        return mylist.index(number+1)
    else:
        return None


def point_to_2point_line_dist(point, line_point):
    """
    Determine the minimum distance from a point to a line defined by two points.
    Based on Wikipedia article: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line, paragraph: Line defined by two points

    :param point: (y,x)
    :param line_point1: (y,x)
    :return:
    """
    x0 = point[1]
    y0 = point[0]
    x1 = line_point[0,1]
    y1 = line_point[0,0]
    x2 = line_point[1,1]
    y2 = line_point[1,0]

    Dx = x2 - x1
    Dy = y2 - y1
    numerator = np.abs(Dy*x0 - Dx*y0 - x1*y2 + x2*y1)
    denom = np.sqrt(Dx**2 + Dy**2)
    distance = numerator/denom
    return distance

def point_line_dist(p, seg, testSegmentEnds=False):
    """
    Minimum Distance between a Point and a Line
    Written by Paul Bourke,    October 1988
    http://astronomy.swin.edu.au/~pbourke/geometry/pointline/

    input:
    p: point, y,x
    seg: y1,x1,y2,x2

    """

    y3, x3 = p
    y1, x1, y2, x2 = seg[0, 0], seg[0, 1], seg[1, 0], seg[1, 1]

    dx21 = (x2 - x1)
    dy21 = (y2 - y1)

    lensq21 = dx21 * dx21 + dy21 * dy21
    if lensq21 == 0:
        #20080821 raise ValueError, "zero length line segment"
        dy = y3 - y1
        dx = x3 - x1
        return sqrt(dx * dx + dy * dy)  # return point to point distance

    u = (x3 - x1) * dx21 + (y3 - y1) * dy21
    u = u / float(lensq21)

    x = x1 + u * dx21
    y = y1 + u * dy21

    if testSegmentEnds:
        if u < 0:
            x, y = x1, y1
        elif u > 1:
            x, y = x2, y2

    dx30 = x3 - x
    dy30 = y3 - y

    return sqrt(dx30 * dx30 + dy30 * dy30)

def point_line_dist_multiline(p, segs, minormax='max'):
    """
    smallest/biggest distance of a point to a sequence of line segments
    """
    if minormax == 'min':
        return min([point_line_dist(p, seg) for seg in segs])
    elif minormax == 'max':
        return max([point_line_dist(p, seg) for seg in segs])

def point2edge_min(image, point):
    """
    Calculates minimum distance from user point to image edges
    point = (y,x)
    """
    rows, cols = size(image)
    disttoedge = np.zeros(4)
    disttoedge[0] = rows - point[0]
    disttoedge[1] = cols - point[1]
    disttoedge[2] = point[0]
    disttoedge[3] = point[1]
    return min(disttoedge)

def size(matrix):
    """
    Matlab equivalent of size; returns the size of the matrix in [rows, columns]
    """
    rows = np.size(matrix, 0)
    cols = np.size(matrix, 1)
    return rows, cols

def invert(matrix):
    """
    Return the imcomplement of the matrix/image. Equivalent to Matlab's imcomplement function.
    """
    newmatrix = -matrix + np.max(matrix) + np.min(matrix)
    return newmatrix

def dist_2points(point1, point2):
    """
    Find the distance from point1 to point2
    """
    #TODO: make this multi-dimensional
    dist = np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)
    return dist

