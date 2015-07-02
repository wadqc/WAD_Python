# PyWAD is an open-source set of plugins for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# This package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN) and Arnold Schilham (UMCU) 
#
#
# Date: 
# Version: 0.2
# Authors: Koen v. Gils, DD
# Changelog:
#
#
# Description of this plugin:
#
# 
#

#####

__version__ = 10062015

import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy import ndimage
import math

import scipy.misc
import dicom
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from skimage.feature import corner_peaks
from skimage import data, img_as_float


def xselector(coordinates,location,interval):
    subsetx=[]
    subsety=[]
    for coordinate in coordinates:
        if location-interval<coordinate[0]<location+interval:
            subsetx.append(coordinate[0])
            subsety.append(coordinate[1])
    return subsetx, subsety



def findcenter(listy,size):
    diff=size
    center=0
    counter=0

    for listitem in listy:
        if abs(listitem-size/2)<diff:
            diff=abs(listitem-size/2)
            center=listitem
            location=counter
        counter+=1
    return location

def fluoroscopydistortion(data,results, **kwargs):
    '''
    Function receives an image of a grid and calculates S-distortion and Pincushion distortion.
    '''
    relevantfile = data.getAllInstances()[0]


     
    #Read DCM file, threshold to separate grid from background
    rawdata=relevantfile
    pixeldata=rawdata.pixel_array
    array_np = np.asarray(pixeldata)
    size=array_np.shape[0]
    if size!=512:
        print 'Resizing image array...'
        array_np=scipy.misc.imresize(array_np,(512,512))
        size=array_np.shape[0]
    forexport=array_np
    mean_pixel=np.mean(pixeldata)
    high_values_indices = array_np > mean_pixel
    low_values_indices = array_np<=mean_pixel
    zero_values_indices = array_np==0
    array_np[high_values_indices]=2
    array_np[low_values_indices] = 1
    array_np[zero_values_indices] = 0
    resultxml = []

    #Define an area for peak detection, resize image matrix if other than 512x512
    area=10

        
    #Plot resized and threholded image
    plt.imshow(array_np)
    plt.show()
    #print (array_np==1).sum()


    #Detect intersections by adjacent cell values
    dimensionrange=range(size)
    xarray=[0]
    yarray=[0]

    for i in dimensionrange:
        for j in dimensionrange:
            
            if array_np[i,j]==1:
                if array_np[i,j-3]==1 and array_np[i,j-2]==1 and array_np[i,j-1]==1 and array_np[i,j+1]==1 and array_np[i,j+2]==1 and array_np[i,j+3]==1 and array_np[i-3,j]==1 and array_np[i-2,j]==1 and array_np[i-1,j]==1 and array_np[i+1,j]==1 and array_np[i+2,j]==1 and array_np[i+3,j]==1:
                    xarray.append(i)
                    yarray.append(j)


    dotsarray=np.zeros((size,size))
    dotsarray[xarray,yarray]=1
    #plt.imshow(dotsarray)
    #plt.show()

    #Detect peaks and plot
    coordinates=corner_peaks(dotsarray, min_distance=area)
    plt.scatter(coordinates[:,1],coordinates[:,0])
    plt.show()

    #Find center row in grid, sort peaks from left to right, find central peak and determine grid rotation
    (subsetx,subsety)=xselector(coordinates,size/2,15)
    (subsety,subsetx)=[list(x) for x in zip(*sorted(zip(subsety, subsetx), key=lambda pair: pair[0]))]
    location=findcenter(subsety,size)
    rotation=math.degrees(math.atan((float(subsetx[location-2])-float(subsetx[location+2]))/(float(subsety[location-2])-float(subsety[location+2]))))
    print rotation

    #Rotate image to correct for grid rotation and plot corrected peaks
    rotateddotsarray=scipy.misc.imrotate(dotsarray, rotation, interp='bilinear')
    rotatedcoordinates=corner_peaks(rotateddotsarray, min_distance=area)

    plt.scatter(rotatedcoordinates[:,1],rotatedcoordinates[:,0])
    plt.show()


    #Find center row in rotated grid
    (subsetx,subsety)=xselector(rotatedcoordinates,size/2,15)
    (subsety,subsetx)=[list(x) for x in zip(*sorted(zip(subsety, subsetx), key=lambda pair: pair[0]))]

    #Determine position of central line to aid in finding consecutive lines in grid
    centralline=np.mean(subsetx)


    #find central peak of rotated grid
    location=findcenter(subsety,size)
    ycenter=subsetx[location]
    xcenter=subsety[location]
    print np.diff(subsety)

    #Determine unit size of grid in center of image
    #unitsize= np.amin(np.nonzero(np.diff(subsety)))
    unitsize= (subsety[location+1]-subsety[location-1])/2
    print unitsize


    #Determine the squares x squares area of grid investigated
    #squares=6
    squares=int(math.floor((0.60*512)/unitsize))/2
    print squares

    #Determine distortion in x direction to aid in finding consecutive lines in grid
    xdistortion=(float(subsety[location+squares])-float(subsety[location-squares]))/(2*squares*unitsize)



    #Determine grid points in lower part of image
    (subsetx,subsety)=xselector(rotatedcoordinates,centralline-squares*unitsize*xdistortion,unitsize/2)
    (subsety,subsetx)=[list(x) for x in zip(*sorted(zip(subsety, subsetx), key=lambda pair: pair[0]))]
    location=findcenter(subsety,size)

    #Determine corners of selected square area
    x1=subsety[location-squares]
    y1=subsetx[location-squares]
    print (x1,y1)
    x3=subsety[location+squares]
    y3=subsetx[location+squares]
    s1=subsetx[location]



    (subsetx,subsety)=xselector(rotatedcoordinates,centralline+squares*unitsize*xdistortion,unitsize/2)
    (subsety,subsetx)=[list(x) for x in zip(*sorted(zip(subsety, subsetx), key=lambda pair: pair[0]))]
    print subsetx,subsety
    location=findcenter(subsety,size)

    x2=subsety[location+squares]
    y2=subsetx[location+squares]
    x4=subsety[location-squares]
    y4=subsetx[location-squares]
    s2=subsetx[location]

 
    pincushion1=float(math.sqrt((x2-x1)**2+(y2-y1)**2)/math.sqrt(2*(2*squares*unitsize)**2))
    pincushion2=float(math.sqrt((x3-x4)**2+(y3-y4)**2)/math.sqrt(2*(2*squares*unitsize)**2))

    pincushionratio=float(pincushion2/pincushion1)

    sartifacttop1=float(s1-y1)
    sartifacttop2=float(s1-y3)
    sartifactbottom1=float(s2-y2)
    sartifactbottom2=float(s2-y4)




    print pincushion1,pincushion2,pincushionratio
    results.addFloat('Pincushion 1', pincushion1,level=1)
    results.addFloat('Pincushion 2', pincushion2,level=1)
    results.addFloat('Pincushion ratio',pincushionratio,level=1)
    results.addFloat('S distortion top 1',sartifacttop1,level=1)
    results.addFloat('S distortion top 2',sartifacttop2,level=1)
    results.addFloat('S distortion bottom 1',sartifactbottom1,level=1)
    results.addFloat('S distortion bottom 2',sartifactbottom2,level=1)

#------------------------- 


    object_naam = 'export.png'


    size=pixeldata.shape[0]
    if size!=512:
        print 'Resizing image array...'
        pixeldata=scipy.misc.imresize(pixeldata,(512,512))

    rotated=scipy.misc.imrotate(pixeldata, rotation, interp='bilinear')
    rotated[y1-3:y1+3,x1-3:x1+3]=255
    rotated[y2-3:y2+3,x2-3:x2+3]=255
    rotated[y3-3:y3+3,x3-3:x3+3]=255
    rotated[y4-3:y4+3,x4-3:x4+3]=255
    rotated[y4-3:y4+3,x4-3:x4+3]=255
    rotated[ycenter-3:ycenter+3,xcenter-3:xcenter+3]=255
    plt.imshow(rotated)
    plt.show()


    scipy.misc.imsave(object_naam,rotated)
    results.addObject('Image',object_naam,level=1)
    return
