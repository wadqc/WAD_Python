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
# This package includes plugins for the automated analysis of QC images for various imaging modalities.
#
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN) and Arnold Schilham (UMCU)
#
#
# Date: 21/05/2014
# Version: 0.1
# Authors: D. Dickerscheid
# Changelog:
#
#
#
'''
 This plugin receives the EARL SUV qc data and makes an entry in the 
 The EARL_SUV tool needs to be run to set the correct values.
'''

import sys,os
import dicom, getopt
import numpy as np
from numpy import random as rnd
import sys

if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def earl_iq_main(data,results, **kwargs):

    results.addChar('status','unprocessed')
    
    floatlist = ['Sphere activity','Calibration time sphere','Background activity','Background act time','Sphere rem act', 'Sphere rem act time', 'Bg rem act', 'Bg rem act time','Scantime','Act Spheres (scantime)',
'Act BG (scantime)','RC-1','RC-2','RC-3','RC-4','RC-5','RC-6']

    for elem in floatlist:
        results.addFloat(elem,'0')

    imlist = ['zprofile','halfslice','rc mean curve','rc max curve']

    for im in imlist:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.title("EARL %s"%str(im))
        img = np.zeros((128,128))
        object_naam = '%s.png'%str(im)
        plt.savefig(object_naam)
        results.addObject('EARL %s'%im,object_naam)









