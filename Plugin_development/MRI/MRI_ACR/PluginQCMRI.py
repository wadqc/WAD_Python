# -*- coding: utf-8 -*-
"""
Created on Fri Dec 05 10:11:57 2014

@author: k.henken,dd
"""

##################### PLUGIN VOOR QC MRI MET ACR FANTOOM #####################

__version__ = 01062015

# IMPORT FUNCTIONS
import os
import numpy as np
import matplotlib
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scipy
import math
import copy
import skimage
from skimage import data, filters, color
from skimage.transform import hough_circle
from skimage.feature import peak_local_max
from skimage.draw import circle_perimeter
from skimage.util import img_as_ubyte

from pyWADLib import findmax 


# IMPORT IMAGES




# FUNCTIONS


def mid_phantom(image_ACR, imarray):
    # Detect edges in image
    edges = skimage.feature.canny(imarray, sigma=3, low_threshold=200, high_threshold=1000)
    hough_radii = np.array([190/2/image_ACR.PixelSpacing[1]])

    print np.shape(edges)
    print hough_radii

    hough_res = hough_circle(edges, hough_radii)
    
    # Detect contours and middle of phantom
    centers = []
    radii = []
    
    for radius, h in zip(hough_radii, hough_res):
        peaks = peak_local_max(h, num_peaks=1)
        centers.extend(peaks)
        radii.extend([radius, radius])
    
    center_x, center_y = centers[0]
    radius = radii[1] # Niet nodig?
    radius = np.int32(radius) # Niet nodig?
    cy, cx = circle_perimeter(center_y, center_x, radius) # Niet nodig?
    return center_x, center_y, radii

def ROI_circ(image_ACR, imarray, a, b, radius):                 #a,b=center circle x,y

    n = image_ACR.Rows                      #size of matrix
    r = radius/image_ACR.PixelSpacing[0]    #radius



    y,x = np.ogrid[-a:n-a, -b:n-b]


    ROI = copy.deepcopy(imarray)
    ROI[x*x + y*y >= r*r] = 0


    ROIm = np.average(ROI, weights=(ROI>0)) # mean PV in ROI
    

    return ROI, ROIm
    
def ROI_ellipse(image_ACR, imarray, a, b, radius, orientation):
    n = image_ACR.Columns
    r = radius
    y,x = np.ogrid[-a:n-a, -b:n-b]
    if orientation == "horizontal":
        mask = x*x/16 + y*y >= r*r
    elif orientation == "vertical":
        mask = x*x + y*y/16 >= r*r
    ROI = copy.deepcopy(imarray)
    ROI[mask] = 0
    ROIm = np.average(ROI, weights=(ROI>0))
    return ROI, ROIm
    
# SLICE THICKNESS ACCURACY
def ACR_STA(header, pixarray, result, label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = header,pixarray
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)  
    # select areas that contain bars
    x_offset = np.int(np.round([4/image_ACR.PixelSpacing[0]]))
    center_x = center_x-1
    imarray_STA = imarray[center_x-x_offset:center_x+x_offset, center_y-100:center_y+100]
    imarray_STA1 = imarray[center_x-x_offset:center_x, center_y-100:center_y+100]
    imarray_STA2 = imarray[center_x:center_x+x_offset, center_y-100:center_y+100]   
    # Thresholding
    threshold = 0.95 * np.mean(imarray_STA)
    imarray_STA11 = imarray_STA1>threshold
    imarray_STA22 = imarray_STA2>threshold
    # Calculate length bars
    len1, len2 = [], []
    for idx1 in range(len(imarray_STA11)):
        len1.append(np.sum(imarray_STA11[idx1,:]))
    for idx1 in range(len(imarray_STA22)):
        len2.append(np.sum(imarray_STA22[idx1,:]))    
    u = np.mean(len1)
    s = np.std(len1)
    len1 = [e for e in len1 if (u - 2 * s < e < u + 2 * s)]
    length1 = np.mean(len1) * image_ACR.PixelSpacing[1]
    u = np.mean(len2)
    s = np.std(len2)
    len2 = [e for e in len2 if (u - 2 * s < e < u + 2 * s)]
    length2 = np.mean(len2) * image_ACR.PixelSpacing[1]
    # Slice thickness
    slice_thickness = 0.2 * (length1*length2)/(length1+length2)
    # Show results
    txt1 = 'Slice thickness is ' + repr(round(slice_thickness,2)) + 'mm (Should be 4.3-5.7 mm)'
    print txt1

    result.addFloat('STA %s'%label,round(slice_thickness,2))
    

#IMAGE INTENSISY UNIFORMITY
def ACR_IIU(header,pixarray,result,label=None):

    # Read image and store pixel values in array
    image_ACR, imarray = header, pixarray
    # Determine midpoint of phantom

    center_x, center_y, radii = mid_phantom(image_ACR, imarray)

    # Determine ROI
    ROI1, ROI1m = ROI_circ(image_ACR, imarray, center_x, center_y, 78)
    print 'hello world'
    # ROI black in white
    temp = np.reshape(ROI1,(1,np.product(ROI1.shape))) #gray values in een rij
    values = np.sort(temp) #gray values van laag naar hoog
    order = np.argsort(temp) #originele positie van gesorteerde gray values
    temp2  = np.where(values>0) #eerste nonzero gray value
    loc = temp2[1][0] #locatie eerst nonzero gray value

    n = np.round(np.pi*(78/image_ACR.PixelSpacing[0])**2/500,0) #aantal benodigde "zwarte punten"
    pixels =  order[0,loc:loc+n] #posities van n "zwarte punten"
    # x,y-positions "zwarte punten"
    pixel_pos = []
    for idx1 in range(len(pixels)):
        row_temp = np.floor(pixels[idx1]/imarray.shape[0])
        col_temp = pixels[idx1]-row_temp*imarray.shape[0]
        pixel_pos.append([row_temp,col_temp])
    pixel_pos = np.array(pixel_pos, dtype=np.int32)  
    #  calculate mean PV/gray value in 1cm2 area around relevant pixels
    PV_mean = []
    for idx2 in range(len(pixel_pos)):
        a, b = pixel_pos[idx2,0], pixel_pos[idx2,1]
        n = image_ACR.Rows
        r = 2/image_ACR.PixelSpacing[0]
        y,x = np.ogrid[-a:n-a, -b:n-b]
        mask = x*x + y*y >= r*r
        ROI = copy.deepcopy(imarray)
        ROI[mask] = 0
        temp = np.average(ROI, weights=(ROI>0))
        PV_mean.append(temp)
    # Determine lowest PV

    PV_low = np.min(PV_mean)
    # ROI white in black
    temp = np.reshape(ROI1,(1,np.product(ROI1.shape))) #gray values in een rij
    order = np.argsort(temp) #originele positie van gesorteerde gray values
    r = 78/image_ACR.PixelSpacing[0] #radius
    n = np.round(np.pi*r*r/500,0) #aantal benodigde "zwarte punten"
    pixels =  order[0,order.shape[1]-n:order.shape[1]] #posities van n "zwarte punten"
    # x,y-positions "witte punten"
    pixel_pos = []
    for idx1 in range(len(pixels)):
        row_temp = np.floor(pixels[idx1]/imarray.shape[0])
        col_temp = pixels[idx1]-row_temp*imarray.shape[0]
        pixel_pos.append([row_temp,col_temp])
    pixel_pos = np.array(pixel_pos, dtype=np.int32)
    #  calculate mean PV/gray value in 1cm2 area around relevant pixels
    PV_mean = []
    for idx2 in range(len(pixel_pos)):
        # circle around point    
        a, b = pixel_pos[idx2,0], pixel_pos[idx2,1]
        n = image_ACR.Rows
        r = 2/image_ACR.PixelSpacing[0]
        y,x = np.ogrid[-a:n-a, -b:n-b]
        mask = x*x + y*y >= r*r
        # set PV's in circle around point to zero
        ROI = copy.deepcopy(imarray)
        ROI[mask] = 0
        # calculate average PV in cicle
        temp = np.average(ROI, weights=(ROI>0))
        PV_mean.append(temp)

    # Determine highest PV
    PV_high = np.max(PV_mean)
    # PIU uitrekenen en printen
    PIU = 100*(1-(PV_high-PV_low)/(PV_high+PV_low))
    txt1 = 'Percent integral uniformity (should be >82% for 3T) is ' + repr(round(PIU,0)) + '%'
    print txt1

    result.addFloat('IIU %s'%label,round(PIU,2))

# GHOSTING
def ACR_Ghosting(header,pixarray,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = header,pixarray
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)
    # Determine mean PV in middle of phantom
    ROI, ROI1m = ROI_circ(image_ACR, imarray, center_x, center_y, 79)
    # Determine mean PV in ROIs around phantom
    ROI, ROI2m = ROI_ellipse(image_ACR, imarray, (center_x-np.int(radii[1]))/2, center_y, (center_x-np.int(radii[1]))/2-2, "horizontal")
    b = center_y+np.int(radii[1])+(image_ACR.Columns-center_y-np.int(radii[1]))/2
    ROI, ROI3m = ROI_ellipse(image_ACR, imarray, center_x, b, image_ACR.Columns-b-2, "vertical")
    ROI, ROI4m = ROI_ellipse(image_ACR, imarray, (center_x+np.int(radii[1]))+(image_ACR.Rows-center_x-np.int(radii[0]))/2, center_y, image_ACR.Rows-((center_x+np.int(radii[1]))+(image_ACR.Rows-center_x-np.int(radii[0]))/2)-2, "horizontal")
    ROI, ROI5m = ROI_ellipse(image_ACR, imarray, center_x, (center_y-np.int(radii[1]))/2, ((center_y-np.int(radii[1]))/2-2), "vertical")
    # Calculate ghosting ratio (GR)
    GR = ((ROI2m+ROI4m)-(ROI1m+ROI3m))/2/ROI1m
    # Print outcome
    txt1 = 'Ghosting ratio (should be <0.025) is ' + repr(round(GR,2))
    print txt1
    
    
    result.addFloat('Ghosting %s'%label,round(GR,2))

# SLICE POSITION ACCURACY
def ACR_SPA(header,pixarray,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = header, pixarray
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)
    # Select ROI
    
    print center_x,center_y

    xa = max(0,center_x-195) #quick fix, weet niet of dit goed is.
    xb = max(0,center_x-90)
    ya = center_y-10
    yb = center_y+10
    imarray_SPA = imarray[xa:xb,ya:yb]
    print imarray_SPA

    # Calculate length of bars
    # Set threshold
    threshold = np.mean(imarray_SPA)

    # Count pixels lower than threshold per vertical cross section 
    length = []
    temp = 0
    for idx1 in range(len(list(imarray_SPA[0,:]))):

        profile = imarray_SPA[:,idx1]
        temp = 0
        for idx2 in range(len(profile)):
            if profile[idx2] < threshold:
                temp = temp+1           
        length.append(temp)
    # Calculate difference length bars in mm
    l_diff_pix=np.mean(length[0:(len(length)/2)])-np.mean(length[(len(length)/2):len(length)])
    l_diff = np.abs(l_diff_pix) * np.float16(image_ACR.PixelSpacing[0])
    # Print outcome
    temp_txt = 'Verschil in lengte van staafjes is ' + repr(round(l_diff,1)) + ' mm'
    print(temp_txt)

    result.wad.addFloat('SPA %s'%label,round(l_diff,1))
    
# HIGH CONTRAST SPATIAL RESOLUTION
def ACR_HCSR(header,pixarray,result,label=None):
    '''
    This routine is broken.
    '''
    # Read image and store pixel values in array
    image_ACR, imarray = header,pixarray
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)

    HCSR = []       
    for idx11 in [1,2,3]:
        ########## Select ROI ##########
        ROI = copy.deepcopy(imarray)
        dxy = 22/image_ACR.PixelSpacing[0]

        if idx11 == 1:
            # ROI1
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y-25/image_ACR.PixelSpacing[0]
        if idx11 == 2:
            # ROI2
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y
        if idx11 == 3:
            # ROI3
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y+21/image_ACR.PixelSpacing[0]
    
        ROI = ROI[x:x+dxy,y:y+dxy]
    
        # plt.imshow(ROI)
    
        # Detect local maxima / corners of grid
        T = np.max(ROI)/2
        #max1 = findmax.find_max(ROI,neighborhood_size=20,threshold=T)
        


        #ROI[max1[1],max1[0]] = 2000
        #plt.imshow(ROI)
    
        # Thresholding
        threshold = np.mean(imarray)*0.8
        for idxx in range(ROI.shape[0]):
            for idxy in range(ROI.shape[1]):
                if ROI[idxx,idxy] < threshold:
                    ROI[idxx,idxy] = 0
                else:
                    ROI[idxx,idxy] = 10
                    #plt.imshow(ROI)
    
        # Count spaces in between dots per row
        dots_row = []
        # pos_tst = []

        for idx1 in range(4):

            profile =  ROI[max1[1][2]-idx1*(max1[1][2]-max1[1][0])/3,max1[0][2]-5:max1[0][2]+15]
            # pos_tst.append(max1[1][2]-idx1*(max1[1][2]-max1[1][0])/3)
            print 'sdfas',profile
            num_dots = 0
            for idx in range(len(profile)):
                if profile[idx]==0 and profile[idx-1]==10:
                    num_dots = num_dots+1
            dots_row.append(num_dots)

        # Count spaces in between dots per column
        dots_col = []
        pos_tst = []
        for idx1 in range(4):
            profile =  ROI[max1[1][4]-14:max1[1][4]+5,max1[0][4]+idx1*(max1[0][5]-max1[0][4])/3]
            pos_tst.append(max1[0][4]+idx1*(max1[0][5]-max1[0][4])/3)
            num_dots = 0
            for idx in range(len(profile)):
                if profile[idx]==0 and profile[idx-1]==10:
                    num_dots = num_dots+1
            dots_col.append(num_dots)
        
        # Store findings
        dots_col = dots_col[dots_col>3]    
        if dots_col >3:
            HCSR.append('Pass')
        else:
            HCSR.append('Fail')
        dots_row = dots_row[dots_row>3]    
        if dots_row >3:
            HCSR.append('Pass')
        else:
            HCSR.append('Fail')
            
    txt1 = 'Resolution 1.1mm - column/row: ' + HCSR[0] + '/' + HCSR[1]
    print txt1
    txt2 = 'Resolution 1.0mm - column/row: ' + HCSR[2] + '/' + HCSR[3]
    print txt2
    txt3 = 'Resolution 0.9mm - column/row: ' + HCSR[4] + '/' + HCSR[5]
    print txt3

    outputdict = { 'Res 1.1 column':HCSR[0],
        'Res 1.1 row':HCSR[1],
        'Res 1.0 column':HCSR[2],
        'Res 1.0 row':HCSR[3],
        'Res 0.9 column':HCSR[4],
        'Res 0.9 row':HCSR[5]
    }


    tmpout = []

    for key in outputdict:
        result.addChar(key+'%s'%label,outputdict[key])




def MRI_ACR_main(data,results, **kwargs):

     acr_test_dict = {
        'ACR_STA':ACR_STA,
         'ACR_IIU':ACR_IIU,
         'ACR_HCSR':ACR_HCSR,
         'ACR_SPA':ACR_SPA,
         'ACR_Ghosting':ACR_Ghosting
    }

     params = kwargs.get('params', None)
     print(data.series_filelist)
     if len(data.getAllInstances())>0:
         tmpfile = data.getAllInstances()[0]
         for child in params:
            print child.text
            if child.text in acr_test_dict.keys():
                try:
                     acr_test_dict[child.text](tmpfile,tmpfile.pixel_array,results)
                except:
                     pass



