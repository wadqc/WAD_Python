# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#

__version__ = '01062015'
__author__ = 'DD'


#
# Description of this plugin:
# This plugin was specifically developed for the analysis of a home made uniformity phantom in the St. Antonius hospital. The phantom consists of a square slab of perspex. The middle of the phantom consists of an aluminum ring with an aluminum ball inside to check if the phantom has been correctly positioned. 
#The plugin looks for the center of the phantom and draws four ROIs relative to the location of the center. For each ROI a mean, std and snr are calculated. 
# Finally the results for each ROI are combined and averaged and a thumbnail of the phantom acquisition is stored.
#


import os
from scipy import ndimage
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import dicom
from dicom import tag
import numpy as np 
from numpy import ma
import sys, getopt

import itertools
import xml.etree.ElementTree as ET

import scipy
import scipy.misc

import operator
import os
from pyWADLib import findmax


def calculate_beam_uniformity (imagearray,roidim=100,neighborhood_size=400):
     ''' Function receives image array and optionally ROI dimension (#pixels). The optional neighborhood_size parameter determines the scale of the distances between local maxima.

    Four ROIs are defined and for each ROI a mean, std and snr is calculated from the image array. The result is returned as a dictionary.
     '''

     thr = imagearray.min()/4.0     
     output = imagearray
     x,y = findmax.find_max(imagearray,neighborhood_size, threshold=thr)

     width = roidim
     height = roidim

     xdim ,ydim = np.shape(imagearray)

     middle = [elem for elem in zip(x,y) if xdim/2 - neighborhood_size < elem[0] < xdim/2 + neighborhood_size and ydim/2 - neighborhood_size < elem[1] < ydim/2 + neighborhood_size]
     
     x0 = int(middle[0][0]/2)
     x1 = int(xdim + middle[0][0])/2
     y0 = int(middle[0][1]/2)
     y1 = int(ydim + middle[0][1])/2
     
     widthx = int(xdim/10)
     widthy = int(xdim/10)

     deltax = int(xdim/10)
     deltay = int(ydim/10)

     roi1 = np.zeros(np.shape(imagearray))
     roi1[x0 - widthx : x0 +widthx, y0 - widthy:y0 + widthy] = 1
     roi1 = ma.make_mask(roi1)

     roi2 = np.zeros(np.shape(imagearray))
     roi2[x1 - widthx : x1 + widthx, y0 - widthy:y0 + widthy]=1
     roi2 = ma.make_mask(roi2)

     roi3 = np.zeros(np.shape(imagearray))
     roi3[x0 - widthx : x0 + widthx, y1 - widthy:y1 + widthy]=1
     roi3 = ma.make_mask(roi3)

     roi4 = np.zeros(np.shape(imagearray))
     roi4[x1 - widthx : x1 + widthx, y1 - widthy:y1 + widthy]=1
     roi4 = ma.make_mask(roi4)



     results = {}

     tmp1 = ma.array(imagearray,mask=1-roi1)
     print ma.count(tmp1)
     results['roi1']={'mean':np.mean(tmp1),'std':np.std(tmp1)}


     tmp2 = ma.array(imagearray,mask=1-roi2)
     results['roi2']={'mean':np.mean(tmp2),'std':np.std(tmp2)}


     tmp3 = ma.array(imagearray,mask=1-roi3)
     results['roi3']={'mean':np.mean(tmp3),'std':np.std(tmp3)}


     tmp4 = ma.array(imagearray,mask=1-roi4)
     results['roi4']={'mean':np.mean(tmp4),'std':np.std(tmp4)}



     
     avgmean = (results['roi1']['mean']+results['roi2']['mean']+results['roi3']['mean']+results['roi4']['mean'])/4.0

     avgstd = (results['roi1']['std']+results['roi2']['std']+results['roi3']['std']+results['roi4']['std'])/4.0

     if avgstd > 0:
          tmpsnr = avgmean/avgstd
     else:
          tmpsnr = -999

     results['avg']={'mean':avgmean,'std':avgstd,'snr':tmpsnr}
     results['middle'] = middle
     results['image'] = imagearray

     return results


def print_beam_output(dcmfile,results):
     '''Function receives a dicom file and results class object. The fuction calls the function calculate_beam_uniformity to calculate mean,std and snr and subsequently writes out the received results into the results class object.
     '''
     
     try:
          detectorname = str(dcmfile.SeriesDescription)
     except:
         detectorname = 'UnknownDetector'

     calc_output = calculate_beam_uniformity(dcmfile.pixel_array)

     results.addFloat('Gemiddelde ROI SNR %s'%detectorname,calc_output['avg']['snr'],level=1)
     results.addFloat('Gemiddelde ROI mean %s'%detectorname,calc_output['avg']['mean'],level=1)
     results.addFloat('Gemiddelde stdev %s'%detectorname,calc_output['avg']['std'],level=1)
     results.addFloat('Gemiddelde ROI1 %s'%detectorname,calc_output['roi1']['mean'],level=2)
     results.addFloat('Gemiddelde ROI2 %s'%detectorname,calc_output['roi2']['mean'],level=2)
     results.addFloat('Gemiddelde ROI3 %s'%detectorname,calc_output['roi3']['mean'],level=2)    
     results.addFloat('Gemiddelde ROI4 %s'%detectorname,calc_output['roi4']['mean'],level=2)

     results.addFloat('stdev ROI1 %s'%detectorname,calc_output['roi1']['std'],level=2 )
     results.addFloat('stdev ROI2 %s'%detectorname,calc_output['roi2']['std'],level=2 )
     results.addFloat('stdev ROI3 %s'%detectorname,calc_output['roi3']['std'],level=2 )    
     results.addFloat('stdev ROI4 %s'%detectorname,calc_output['roi4']['std'],level=2 )

     print('CWD:',os.getcwd())
     fig = plt.figure()
     ax = fig.add_subplot(111)
     plt.title("QC bucky")
     img = calc_output['image']
     middlex, middley = zip(*calc_output['middle'])
     plt.imshow(img)
     plt.scatter(middlex,middley)
     plt.savefig('%s.png'%detectorname.replace(" ",""))
     results.addObject('QCbeeld','%s.png'%detectorname.replace(" ",""))
     print('CWD',os.getcwd())     




def QC_bucky_run(data, results, **kwargs):
     '''Function extracts instances from data object. From the params section in the config XML file the filter for relevant bucky files is determined. For each relevant file the print_beam_output function is called.
     '''

     paramdict = kwargs.get('params', None) #read out all the parameter tags from the config_xml

     select_instance_lst = [] #create a list of all buckydata to be processed
     for child in paramdict:
          print(child, type(child), child.tag, type(child.tag))
          if 'bucky' in child.tag:
               select_instance_lst.append(child.attrib) 

     for crit in select_instance_lst:
          tmpfile = data.getInstanceByTags(crit)
          for elem in tmpfile:
               try:
                    print_beam_output(elem,results)
               except:
                    print("Warning, failed %s "%crit)

if __name__ == "__main__":
    from pyWAD import PluginTools
    
    series_lst = [[],[]] #to test from the command line add test files from a series here
    result_dir = "./"
    
    data = PluginData()
    results = PluginResults()
    QC_bucky_run()
    
    for result in tools.results:
        print result
