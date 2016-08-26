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

# PyWAD is an open-source set of plugins for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Date: 01062015
# Version: 1.0
# Authors: DD
# Changelog:
#
#
# Description of this plugin:
# 
#
'''
A very crude and simple implementation to estimate the MTF of an Edge phantom.
'''


import os,sys
import numpy as np
import scipy 

import dicom
import pylab as plt
from scipy import ndimage
from skimage import filter
from skimage.feature import corner_harris, corner_peaks
from skimage import measure
from itertools import izip, tee



def filter_min_diff(inlist):
    x_min_diff = 10
    outlist = []
    if len(inlist) > 0:
        outlist.append(inlist[0])    
    return outlist

def diff_list(rowlist):
    return np.diff(rowlist)




def MTF_main(data, results, **kwargs):
    '''
    Function estimates MTF from an exposure of an edge phantom.
    Workflow:
    1. From data object obtain first instance
    2. Remove an arbitrarily chosen offset, which is hardcoded (20 pixels). This needs improvement.
    3. Define what the minimal stepsize between adjacent pixels is
    4. For each row (hardcoded range!, needs improvement) determine where the step is located
    5. Sum the edge curves and finally take the Fourier transform to calculate the MTF
    '''    
    
    relevantfile = data.getAllInstances[0]
    print "relevantfile:", relevantfile

    dicomobject = relevantfile
    pixel_map = dicomobject.pixel_array

    xdim,ydim = np.shape(pixel_map)
    
    offset = 20
    fmin = np.min(pixel_map[offset:xdim - offset,offset:ydim-offset])
    fmax = np.max(pixel_map[offset:xdim - offset,offset:ydim-offset])

    min_stepsize = float((fmax - fmin)/8.0)
    nr_edge_px = 10
    edge_loc = []

    testline = np.zeros(np.shape(pixel_map)[0]-1 - 2*offset)
    
    for row in range(1400,2400):#ydim):
       curr_row = np.array(pixel_map[offset:xdim-offset ,row],dtype=float)
       tmpdiff = np.diff(curr_row)


       ## we determine the x and y coordinates of the steps
       tmp_edge_loc = [ [n,row] for n,v in enumerate(tmpdiff) if abs(v) > min_stepsize]

       ## there could be multiple steps in a row and we only want the first ones
       tmp_edge_loc = filter_min_diff(tmp_edge_loc)

       for elem in tmp_edge_loc:
           edge_loc.append(elem)

       testline = np.add(testline,np.diff(curr_row))

        
    sum_edge_curves = np.zeros(2*nr_edge_px - 2)
    
    norm = 1
    for edge in edge_loc:
        norm+=1
        tmplist =  np.zeros(2*nr_edge_px - 1,dtype=float)
        for m in range(2*nr_edge_px - 1):
            tmplist[m] = pixel_map[edge[0] + offset - nr_edge_px + m,edge[1]]
        
        diff_tmplist = np.diff(tmplist)
        sum_edge_curves = np.add(sum_edge_curves, diff_tmplist)
            
        

    sum_edge_curves = sum_edge_curves/norm
    sum_edge_curves_fft = np.fft.fft(sum_edge_curves)
    sum_edge_curves_fft = sum_edge_curves_fft[range(len(sum_edge_curves)/2)]/sum_edge_curves_fft[0]

    frq = (2*nr_edge_px)*np.arange(len(sum_edge_curves))/len(sum_edge_curves) # two sides frequency range
    frq = frq[range(len(sum_edge_curves)/2)] # one side frequency range

    object_naam = 'MTFim.png'
    img = pixel_map
    plt.savefig(object_naam_pad,img)
    results.addObject('MTF image',object_naam,level='1')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title("MTF")
    print "plot abs(sum edge curves)"
    plt.plot(frq,np.abs(sum_edge_curves_fft))
    plt.xlabel('lp/mm')
    plt.ylabel('Normalized modulation factor')
    plt.xlim([0,10])
    plt.grid()
    
    object_naam = 'MTFcurve.png'
    plt.savefig(object_naam,bbox_inches=0)
    results.addObject(    'MTF curve',object_naam,level='1')




    return resultxml

    









if __name__ == "__main__":
    sys.exit(main())
