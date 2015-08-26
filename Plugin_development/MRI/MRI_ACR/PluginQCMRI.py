# -*- coding: utf-8 -*-
'''
This plugin is an implementation of the ACR quality control protocol for MRI as specified in the following documents published by the ACR:

http://www.acr.org/~/media/ACR/Documents/Accreditation/MRI/LargePhantomInstructions.pdf
http://www.acr.org/~/media/ACR/Documents/Accreditation/MRI/LargePhantomGuidance.pdf

This code contains the following routines that can be called:
ACR_GA_l: Geometric accuracy - measure length of phantom
ACR_GA_wh: Geometric accuracy - measure width of phantom
ACR_GA_whd: Geometric accuracy - measure length of diagonal
ACR_STA: Slice Thickness Accuracy
ACR_IIU: Image intensity uniformity
ACR_HCSR: High contrast spatial resolution
ACR_SPA: Slice position accuracy
ACR_Ghosting: Ghosting

Low contrast object detectability (LCOD) and Percent-signal ghosting (PSG) are currently not implemented.

#Changelog:
26-08-15 - Updated original code to a working wadplugin version (DD)
05-12-14 - First version

'''

##################### PLUGIN VOOR QC MRI MET ACR FANTOOM #####################

__version__ = '01062015'
__author__ = 'k.henken'



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
    edges = filter.canny(imarray, sigma=3, low_threshold=200, high_threshold=1000)

    hough_radii = np.array([190/2/image_ACR.PixelSpacing[1]])
    print type(edges)
    print type(hough_radii)
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
#    imarray[center_x, center_y] = 10000
#    imarray[cx,cy] = 10000
#    plt.imshow(imarray)
    return center_x, center_y, radii


def ROI_circ(image_ACR, imarray, a, b, radius):                 #a,b=center circle x,y

    n = image_ACR.Rows                      #size of matrix
    r = radius/image_ACR.PixelSpacing[0]    #radius
    y,x = np.ogrid[-a:n-a, -b:n-b]
    mask = x*x + y*y >= r*r
    ROI = copy.deepcopy(imarray)
    ROI[mask] = 0
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
    

# GEOMETRIC ACCURACY
def ACR_GA_l(image,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image)
    ########## Measure length of phantom ##########    
     # set threshold
    threshold = 0.12*np.mean(imarray)
    
    pvsumloc = []
    pvsum = sum(imarray.T)
    for idx in range(len(imarray[:,0])):
        if pvsum[idx] > threshold*image_ACR.Rows:
            pvsumloc.append(idx)
    
    # calculate heigth (last - first column)
    height = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[0]
    
    txt_temp = 'length is '+repr(round(height,1))+' mm (should be 146-150 mm)'
    print(txt_temp)

def ACR_GA_wh(image,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image) 
    ########## Measure width of phantom ##########
    # set threshold
    threshold = 0.12*np.mean(imarray)
    
    # columns with a summed pixel value > threshold
    pvsumloc = []
    pvsum = sum(imarray)
    for idx in range(len(imarray[0,:])):
        if pvsum[idx] > threshold*image_ACR.Columns:
            pvsumloc.append(idx)
    
    # calculate width (last - first column)
    width = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[1]
    
    ########## Measure height of phantom ##########
    # columns with a summed pixel value > threshold
    pvsumloc = []
    pvsum = sum(imarray.T)
    for idx in range(len(imarray[:,0])):
        if pvsum[idx] > threshold*image_ACR.Rows:
            pvsumloc.append(idx)
    
    # calculate heigth (last - first column)
    height = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[0]
    
    txt_temp = 'width is '+repr(round(width,1))+' mm (should be 188-192 mm)'
    print(txt_temp)
    txt_temp = 'height is '+repr(round(height,1))+' mm (should be 188-192 mm)'
    print(txt_temp)

def ACR_GA_whd(image,result):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image)
        ########## Measure width of phantom ##########
    # set threshold
    threshold = 0.12*np.mean(imarray)
    
    # columns with a summed pixel value > threshold
    pvsumloc = []
    pvsum = sum(imarray)
    for idx in range(len(imarray[0,:])):
        if pvsum[idx] > threshold*image_ACR.Columns:
            pvsumloc.append(idx)
    
    # calculate width (last - first column)
    width = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[1]
    
    ########## Measure height of phantom ##########
    # columns with a summed pixel value > threshold
    pvsumloc = []
    pvsum = sum(imarray.T)
    for idx in range(len(imarray[:,0])):
        if pvsum[idx] > threshold*image_ACR.Rows:
            pvsumloc.append(idx)
    
    # calculate heigth (last - first column)
    height = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[0]
    
    txt_temp = 'width is '+repr(round(width,1))+' mm (should be 188-192 mm)'
    print(txt_temp)
    txt_temp = 'height is '+repr(round(height,1))+' mm (should be 188-192 mm)'
    print(txt_temp)

    ########## Rotate image to measure diagonal ##########
    imarrayrot = scipy.misc.imrotate(imarray,45)
    plt.imshow(imarrayrot)
    
    ########## Measure diagonal 1 of phantom ##########
    # set threshold
    threshold = 0.12*np.mean(imarrayrot)
    
    # columns with a summed pixel value > threshold
    pvsumloc = []
    #pvsum = sum(imarrayrot.T)
    for idx in range(len(imarrayrot[:,0])):
        pvsum = sum(imarrayrot[idx,:])
    #    print pvsum
        if pvsum > threshold*image_ACR.Rows:
            pvsumloc.append(idx)
    
    # calculate heigth (last - first column)
    diag1 = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[0]
    
    ########## Measure diagonal 2 of phantom ##########
    # columns with a summed pixel value > threshold
    pvsumloc = []
    #pvsum = sum(imarrayrot)
    for idx in range(len(imarrayrot[0,:])):
        pvsum = sum(imarrayrot[:,idx])
        if pvsum > threshold*image_ACR.Columns:
            pvsumloc.append(idx)
    
    # calculate heigth (last - first column)
    diag2 = (np.max(pvsumloc)-np.min(pvsumloc)) * image_ACR.PixelSpacing[0]
    
    txt_temp = 'diagonal 1 is '+repr(round(diag1,1))+' mm (should be 188-192 mm)'
    print(txt_temp)
    txt_temp = 'diagonal 2 is '+repr(round(diag2,1))+' mm (should be 188-192 mm)'
    print(txt_temp)


# SLICE THICKNESS ACCURACY
def ACR_STA(header, pixarray, result, label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image)
#    plt.figure(1)
#    plt.imshow(imarray)
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)  
    # select areas that contain bars
    x_offset = np.int(np.round([4/image_ACR.PixelSpacing[0]]))
    center_x = center_x-1
    imarray_STA = imarray[center_x-x_offset:center_x+x_offset, center_y-100:center_y+100]
    imarray_STA1 = imarray[center_x-x_offset:center_x, center_y-100:center_y+100]
    imarray_STA2 = imarray[center_x:center_x+x_offset, center_y-100:center_y+100]   
    # Thresholding
#    threshold = 0.95 * np.mean(imarray_STA)
    threshold = 0.9 * np.mean(imarray_STA)
    imarray_STA11 = imarray_STA1>threshold
    imarray_STA22 = imarray_STA2>threshold
#    plt.figure(2)        
#    plt.imshow(imarray_STA11)
#    plt.figure(3)        
#    plt.imshow(imarray_STA22)
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
    image_ACR, imarray = load_image(image)
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)
    # Determine ROI
    ROI1, ROI1m = ROI_circ(image_ACR, imarray, center_x, center_y, 78)
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
    image_ACR, imarray = load_image(image)
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
#    print ROI1m
#    print ROI2m
#    print ROI3m
#    print ROI4m
#    print ROI5m

    # Calculate ghosting ratio (GR)
    GR = np.abs(((ROI2m+ROI4m)-(ROI3m+ROI5m))/(2*ROI1m)) #--> hier klopt iets niet!!!!!!!!!!!!!!!!!1
    # Print outcome
    txt1 = 'Ghosting ratio (should be <0.025) is ' + repr(round(GR,3))
    print txt1
    

    
    result.addFloat('Ghosting %s'%label,round(GR,2))

# SLICE POSITION ACCURACY
def ACR_SPA(header,pixarray,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image)
#    plt.figure(1)
#    plt.imshow(imarray)
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)
    # Select ROI
#    imarray_SPA = imarray[center_x-195:center_x-90, center_y-10:center_y+10]
    imarray_SPA = imarray[center_x-(95/image_ACR.PixelSpacing[0]):center_x-(44/image_ACR.PixelSpacing[0]), center_y-(4.8/image_ACR.PixelSpacing[0]):center_y+(4.8/image_ACR.PixelSpacing[0])]
#    plt.figure(2)
#    plt.imshow(imarray_SPA)
#    print image_ACR.PixelSpacing[0] 
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
    temp_txt = 'Verschil in lengte van staafjes is ' + repr(round(l_diff,1)) + ' mm (should be <5 mm)'
    print(temp_txt)

    result.wad.addFloat('SPA %s'%label,round(l_diff,1))
    
# HIGH CONTRAST SPATIAL RESOLUTION
def ACR_HCSR(header,pixarray,result,label=None):
    # Read image and store pixel values in array
    image_ACR, imarray = load_image(image)   
    plt.figure(1)    
    plt.imshow(imarray)
    
    # Determine midpoint of phantom
    center_x, center_y, radii = mid_phantom(image_ACR, imarray)
      
    HCSR = []       
    for idx11 in [1,2,3]:
#    for idx11 in [1]:
        ########## Select ROI ##########
        ROI = copy.deepcopy(imarray)
#        plt.imshow(ROI)
        dxy = 22/image_ACR.PixelSpacing[0]
        
        if idx11 == 1:
            # ROI1
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y-25/image_ACR.PixelSpacing[0]
#            ROI[x,y] = 10000
#            ROI[x+dxy,y+dxy] = 10000
#            plt.imshow(ROI)
        if idx11 == 2:
            # ROI2
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y
        if idx11 == 3:
            # ROI3
            x = center_x+28/image_ACR.PixelSpacing[0]
            y = center_y+21/image_ACR.PixelSpacing[0]
    
#        print image_ACR.PixelSpacing[0]       
#        print image_ACR.PixelSpacing[1]
        ROI = ROI[x:x+dxy,y:y+dxy]
#        plt.imshow(ROI)
    
        # Detect local maxima / corners of grid
#        T = np.max(ROI)/2
#        T = 0.5*np.max(ROI)
##        max1 = find_max(ROI,neighborhood_size=20,threshold=T)
#        max1 = find_max(ROI,neighborhood_size=20,threshold=T)
##        print max1
##        ROI[max1] = 10000
#        mid_tst = scipy.ndimage.measurements.center_of_mass(ROI)
#        print mid_tst
#        ROI[mid_tst] = 10000
#        plt.imshow(ROI)
         
        #ROI[max1[1],max1[0]] = 2000
        #plt.imshow(ROI)
    
        # Thresholding
        threshold = np.mean(imarray)*0.7
        #        threshold = np.mean(imarray)*0.8
        for idxx in range(ROI.shape[0]):
            for idxy in range(ROI.shape[1]):
                if ROI[idxx,idxy] < threshold:
                    ROI[idxx,idxy] = 0
                else:
                    ROI[idxx,idxy] = 1
                    
        # Detect midpoint of grid
        row_val = []
        col_val = []
        for idxx in range(ROI.shape[0]):
            temp = ROI[idxx].sum()
            row_val.append(temp)
        midx = np.argmax(row_val)
        for idxy in range(ROI.shape[1]):
            temp = ROI[:,idxy].sum()
            col_val.append(temp)
        midy = np.argmax(col_val)
        
#        print midx
#        print midy
#        ROI[midx, midy] = 2
        
#        plt.imshow(ROI)
        # Count spaces in between dots per row
        dots_row = []
        # pos_tst = []
        for idx1 in range(4):
#            profile =  ROI[max1[1][2]-idx1*(max1[1][2]-max1[1][0])/3,max1[0][2]-5:max1[0][2]+15]
            profile = ROI[midx-idx1*2.1/image_ACR.PixelSpacing[0],midy-8.1/image_ACR.PixelSpacing[1]:midy+1]     
#            ROI[midx-idx1*2.1/image_ACR.PixelSpacing[0],midy-6.8/image_ACR.PixelSpacing[1]:midy] = 2
            # pos_tst.append(max1[1][2]-idx1*(max1[1][2]-max1[1][0])/3)
#            print profile
            num_dots = 0
            for idx in range(len(profile)):
                if profile[idx]==0 and profile[idx-1]==1:
                    num_dots = num_dots+1
            dots_row.append(num_dots)
#        print dots_row
#        plt.imshow(ROI)
        
        
        # Count spaces in between dots per column
        dots_col = []
#        pos_tst = []
        for idx1 in range(4):
#            profile =  ROI[max1[1][4]-14:max1[1][4]+5,max1[0][4]+idx1*(max1[0][5]-max1[0][4])/3]
            profile =  ROI[midx:midx+9.2/image_ACR.PixelSpacing[1],midy+idx1*2.2/image_ACR.PixelSpacing[0]]
#            ROI[midx:midx+9.2/image_ACR.PixelSpacing[1],midy+idx1*2.2/image_ACR.PixelSpacing[0]] = 2
#            pos_tst.append(max1[0][4]+idx1*(max1[0][5]-max1[0][4])/3)
#            print profile
            num_dots = 0
            for idx in range(len(profile)):
                if profile[idx]==0 and profile[idx-1]==1:
                    num_dots = num_dots+1
            dots_col.append(num_dots)
#        print dots_col
#        plt.figure(2)        
#        plt.imshow(ROI)
        
        # Store findings
        dots_col = [dot for dot in dots_col if dot>3]
        
#        print 'hiero        '
#        print dots_col
#        print len(dots_col)
        if len(dots_col) >0:
            HCSR.append('Pass')
        else:
            HCSR.append('Fail')
        dots_row = [dot for dot in dots_row if dot>3]
#        print dots_row
        if len(dots_row) >1:
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
         'ACR_GA_l':ACR_GA_l,
         'ACR_GA_wh':ACR_GA_wh,
         'ACR_GA_whd':ACR_GA_whd,
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



