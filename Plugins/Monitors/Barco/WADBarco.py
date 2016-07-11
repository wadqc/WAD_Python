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
# This package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN) and Arnold Schilham (UMCU) 
#
#
# Date: 
# Version: 0.2
# Authors: K. van Gils, D. Dickerscheid
# Changelog:
#
#
# Description of this plugin:
#
# 
#


__version__=01062015

import os, sys
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
#import pylab as plt
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import csv
import StringIO
import dicom
from dicom import tag
from datetime import datetime

def convertboolean(inputboolean):
    if inputboolean:
        output='1'
    else:
        output='0'
    return output


def barco_compliance(data, results, **kwargs):
    
    dcmlist = data.getAllInstances()
    if len(dcmlist) > 1:
        print "only one instance supported!"
        sys.exit()

    dicomobject = dcmlist[0]
    print dicomobject

    csvfile=dicomobject[tag.Tag("0071","0071")].value
    csvfileread=StringIO.StringIO(csvfile)
    data = csv.reader(csvfileread, delimiter='\t')

    contrastlijst=[]
    contrastlijsttarget=[]
    contrasterrorlijst=[]
    counter=0
    countertarget=0
    resultxml = []
    for row in data:

            info=row[0].split(' ')
            datum=''.join(info[0:3])
            scherm=info[6].split('-')
            scherm=scherm[2]
            try:
                print info[10]
                if info[9]=="contrast" and info[11]=="error":
                    data=info[-1].strip('%')
                    contrasterrorlijst.append(float(data))
            except IndexError:
                pass
            
            try:
                if info[9]=="Ambient" and info[10]=="light:":
                    AmbientLightvalue=info[11]
                    results.addFloat('Omgevingslicht',AmbientLightvalue,level=1)
            except IndexError:
                pass

            try:
                if info[9]=="Error" and info[10]=="tolerance:":
                    Barcotolerance=info[11]
                    results.addFloat('Barco contrast tolerantiewaarde',Barcotolerance,level=1)
            except IndexError:
                pass


            try:     
                if info[9]=="Measured" and info[10]=="Contrast":
                    counter=17
            except IndexError:
                pass

            try:              

                 if info[9]=="contrastResponse:" and counter>0:
                    contrastlijst.append(float(info[-1]))
                    counter-=1
            except IndexError:
                 pass

            try:
                 if info[9]=="Target" and info[10]=="Contrast":
                    countertarget=17
                    print 'hoi'
            except IndexError:
                 pass
            try:
                if info[9]=="contrastResponse:" and countertarget>0:
                    contrastlijsttarget.append(float(info[-1]))
                    countertarget-=1
            except IndexError:
                pass

            try:
                if info[13]=="Status:":
                    if info[14]=="true":
                        DICOMluminance=True
                    else:
                        DICOMluminance=False
            except IndexError:
                pass

            try:
                if info[9]=="Result:":
                    if info[10]=="OK":
                        Result=True
                    else:
                        Result=False
            except IndexError:
                pass


    

    tolerance=10
    error=True

    for contrasterror in contrasterrorlijst:
        if contrasterror>=tolerance:
            error=False


    print contrastlijst

    print contrastlijsttarget

    print contrasterrorlijst

    results.addBool('Contrast eigen tolerantie',error,level=1)

    print DICOMluminance
    results.addBool('Barco Luminance',DICOMluminance,level=1)

    print Result
    results.addBool('Barco Result',Result,level=1)

    helpcounter=0
    for row in contrastlijst:
        results.addFloat('Contrast '+str(helpcounter),row,level=1)
        helpcounter+=1



#-------------------------    
    object_naam = 'GSDFcurve.png' 

    x=range(17)

    toleranties=[y*tolerance/100 for y in contrastlijsttarget]

    tolerantiehoog = []
    tolerantielaag = []

    for index, item in enumerate(contrastlijsttarget):
        tolerantiehoog.append(contrastlijsttarget[index] + toleranties[index])
        tolerantielaag.append(contrastlijsttarget[index] - toleranties[index])


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,contrastlijst,'b.',x,tolerantiehoog,'r--',x,tolerantielaag,'r--',x,contrastlijsttarget,'b:')
    plt.suptitle('Display ' + scherm+' '+datum)
    plt.yscale('log')
    plt.ylabel('dL/L')
    plt.grid(True)
    plt.xlabel('Measurement point')
    plt.savefig(object_naam)

    results.addObject('GSDF curve',object_naam,level=1)


    #object_naam = 'csvfile.csv'
    #results.addObject('GSDF curve data',object_naam,level=2)




