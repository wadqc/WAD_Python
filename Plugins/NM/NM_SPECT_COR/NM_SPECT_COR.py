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

# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Description of this plugin:
#
# This program reads projection data and 
# determines the COR shift.
#
__version__ = '01062015'
__author__ = 'DD'




import sys,os
import dicom, getopt
from dicom import tag
import numpy as np
from numpy import random as rnd

import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

#import pylab as plt
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.pyplot

import sys


from scipy import optimize
from scipy.optimize import curve_fit

from scipy.stats import norm
import matplotlib.mlab as mlab
from pyWADLib import findmax
import os


def gaussian(x,a,mu,sigma):
    ''' Returns the gaussian function for B=m,stdev,max,offset '''
    return a*np.exp(-((x-mu)**2/(2*sigma**2)))



def sinogramfit(x, a, b, c,d):
    return a*np.sin(b*x+c) + d

class Datatree(file):
 
    
    def __init__(self, file):
        self.dcmfile = file
        self.pixelmap = file.pixel_array
        self.zdim, self.xdim, self.ydim = np.shape(self.pixelmap)

        ## get the locations of the point sources from the first projection
        
        self.bandwidth = 20
        self.sources = self.findmax(frame_id=0,neighborhood_size=self.bandwidth)
        
        ## the "y" component contains the number and location of sources
        ## The function show sinogram shows the sinogram for these sources
        


        

    def fitted_gauss(self,ydata):
        
        x = np.linspace(0,self.xdim, self.xdim)
        y = ydata
        
        p0 = [1,100,10]

        popt, pcov = curve_fit(gaussian,x,y,p0)
        return popt



        
        
    def show_sinogram(self):

        

        
        for source in self.sources[1]:
#            print "source =", source

            tmp_sinogram_data = np.sum(self.pixelmap[:,source-self.bandwidth/2:source+self.bandwidth/2,:],1)
            

            ## First detector

            det1data = range(0,self.zdim/2)

            ## Second detector

            det2data = range(self.zdim/2,self.zdim)

            detlist = [det1data,det2data]
            detdict = {}

            for detector in detlist:
                source_sinograms = {}

                tmpx = []
                detidx = detlist.index(detector) + 1
                print "DETECTOR", detidx
                print "\n"
                for angle in detector:
                    tmpx.append(self.fitted_gauss(tmp_sinogram_data[angle,:])[1])


                source_sinograms[source] = tmpx
                detdict[detidx] = source_sinograms
                
        return detdict

            


    def show_histogram(self, frame_id, pixelmap=None):

        if pixelmap == None:
            pixelmap = self.pixelmap

        pylab.hist(pixelmap[frame_id])
        return

    def findmax(self,frame_id,neighborhood_size):
        x,y = findmax.find_max(self.pixelmap[frame_id],neighborhood_size)
        return x,y


        

def spect_cor_main(data,results, **kwargs):

    print "Starting SPECT cor analysis"
    relevantfile = data.getAllInstances()[0]
    pixel_map = relevantfile.pixel_array
    pixelsize = float(relevantfile[tag.Tag("0028","0030")].value[0])
  
    data = Datatree(relevantfile)
    detdict = data.show_sinogram()

    corlist = []
    print "detdict.keys()",detdict.keys()
    for det in detdict.keys():
        print det, detdict[det]
        source_sinograms = detdict[det]
        for key in source_sinograms.keys():
            print "key=", key, source_sinograms[key]

            x = np.linspace(0,2*np.pi, data.zdim/2)

            y = source_sinograms[key]
        
            p0 = [60.,1/360.,0.,data.zdim/2.]

            popt, pcov = curve_fit(sinogramfit,x,y,p0)

            tmp_cor =  ((popt[3]+1)-((data.xdim+1)/2.0))*pixelsize
            corlist.append(tmp_cor)

            yfit =  sinogramfit(x,*popt)
           
            fig = plt.figure(str(det))
            ax = fig.add_subplot(111)
            plt.title("detector "+str(det))
            ax.plot(y,'bo')
            ax.plot(yfit)
            plt.text(0, 100, "COR$_{avg}$ = %.2f"%np.mean(corlist), size=25, rotation=0.,ha="left", va="top", bbox = dict(boxstyle="square",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),))                 

            object_naam = 'Sinogram_det_%s.png'%str(det)   
            plt.savefig(object_naam)

            results.addObject('Sinogram det %s'%str(det),object_naam,level=1)



        results.addFloat('COR det %s'%str(det),np.mean(corlist),level=1 )

if __name__ == "__main__":
    sys.exit(main())

