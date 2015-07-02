# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Changelog:
#
#
# Description of this plugin:
# This plugin calculates the NEMA differential and integral uniformity of an image 
#


__version__='01062015'
__author__ = 'DD, tdw'



import sys,os
import dicom, getopt
from dicom import tag
import numpy as np
from numpy import random as rnd

import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

from pyWAD.plugindata import PluginData
from pyWAD.pluginresults import PluginResults
import NEMA_unif_lib as nemalib

import os


class uniformity_nm:
    def __init__(self, data, results, **kwargs):

        self.params = kwargs.get('params', None)

        p = {}
        for param in self.params:
            p[param.tag] = (param.text,param.attrib)
        
        self.unif1 = p.get('unif1')
        self.unif2 = p.get('unif2')
        try:
            self.dome_corr = p.get('perform_dome_correction')[0]
        except:
            self.dome_corr = False

        self.dome_corr = self.dome_corr in ['yes','ja','1','true']

        print "Perform dome correction = ", self.dome_corr

        print "unif and lin matching criteria:"
        print self.unif1
        print self.unif2


        # data comes either as:
        #   - two series of each 1 instance (1 instance for each detector)
        #   - one series, one instance (multiframe: frame0=detector1, frame1=detector2)
        #   - one series, two instances (1 instance for each detector)
        #
        # FIXME: better handling of all three possibilities (incl naming of detector labels)
        serieslist = data.getAllSeries()

 
        for series in serieslist:
            print "Current series: ", series[0].SeriesDescription

            #instances = series[1]

            for instance in series:
                print instance.SeriesDescription
  
                dicomobject = instance
                tmpseriesdescription =  dicomobject[tag.Tag("0008","103E")].value

                try:
                    if self.unif1[1]['seriesdescription'] == tmpseriesdescription:
                        print "UNIF1 series Match!"
                        tmp = self.unif1[1]['instancetag'].split(',') 
                        tmpinstancetag = tag.Tag(tmp[0],tmp[1])
                        print "instance tag: ", tmp
                        print "instance tag val: ", dicomobject[tmpinstancetag].value
                        print '--',self.unif1[0]


                        if self.unif1[0] == str(dicomobject[tmpinstancetag].value):
                            print "UNIF1 instance Match!"

                            try:
                                if int(dicomobject[tag.Tag("0028","0008")].value) > 1:
                                    print "multiframe! "
                
                                    bits = int(dicomobject[tag.Tag("0028","0100")].value)
                                    frames = int(dicomobject[tag.Tag("0028","0008")].value)
                                    rows = int(dicomobject[tag.Tag("0028","0010")].value)
                                    columns = int(dicomobject[tag.Tag("0028","0011")].value)

                                    count = rows*columns*frames

                                    if bits == 8:
                                        dtype = np.uint8
                                    elif bits == 16:
                                        dtype = np.uint16
                                    elif bits == 32:
                                        dtype = np.uint32

                                    tmppixelmap = np.reshape(np.fromstring(dicomobject.PixelData,count=count,dtype=dtype),(frames,rows,columns))

                                    for n in range(frames):
                                        print "Frame: ", n                                    
                                        tmpresult = self.Uniformity_main(tmppixelmap[n],results,'det%s'%str(n+1),self.dome_corr)

                                else:
                                    frames = 1
                                    tmppixelmap = dicomobject.pixel_array
                                    tmpresult = self.Uniformity_main(tmppixelmap,results,'det1',self.dome_corr)

                            except:
                                pass

                except:
                    pass

                try:   
                    if self.unif2[1]['seriesdescription'] == tmpseriesdescription:
                        print "UNIF2 series Match!"
                        tmp = self.unif2[1]['instancetag'].split(',')
                        tmpinstancetag = tag.Tag(tmp[0],tmp[1])
                        if self.unif2[0] == str(dicomobject[tmpinstancetag].value):
                            print "UNIF2 instance Match!"

                            try:
                                if int(dicomobject[tag.Tag("0028","0008")].value) > 1:
                                    print "multiframe! "
                
                                    bits = int(dicomobject[tag.Tag("0028","0100")].value)
                                    frames = int(dicomobject[tag.Tag("0028","0008")].value)
                                    rows = int(dicomobject[tag.Tag("0028","0010")].value)
                                    columns = int(dicomobject[tag.Tag("0028","0011")].value)

                                    count = rows*columns*frames

                                    if bits == 8:
                                        dtype = np.uint8
                                    elif bits == 16:
                                        dtype = np.uint16
                                    elif bits == 32:
                                        dtype = np.uint32

                                    tmppixelmap = np.reshape(np.fromstring(dicomobject.PixelData,count=count,dtype=dtype),(frames,rows,columns))

                                    for n in range(frames):
                                        print "Frame: ", n
                                        tmpresult = self.Uniformity_main(tmppixelmap[n],results,'det%s'%str(n+1),self.dome_corr)
                                        #self.resultxml.append(tmpresult)

                                else:
                                    frames = 1
                                    tmppixelmap = dicomobject.pixel_array
                                    tmpresult = self.Uniformity_main(tmppixelmap,results,'det2',self.dome_corr)
                                    #self.resultxml.append(tmpresult)

                            except:
                                pass

                except:
                    pass


                print "End of class init"



    def check_dicomdict_studies(self,dicomdict):
        if len(dicomdict) > 1:
            print "dicomdict contains more than 1 study! aborting!"
            sys.exit()


    def savewadfig(self,impath,pixelmap):

        try:
            plt.imsave(impath,pixelmap)

        except:
            print "could not save image"

        return


    def Uniformity_main(self, pixelmap, results, label, dome_corr):

        print 'Starting Uniformity_main'
        detector = label
        print 'Calculating results'
        print 'pixelmap size: ', np.shape(pixelmap)
        
        output = nemalib.calculate_nema_uniformity(pixelmap, (64,64), results, dome_corr)
        
        '''
          output[0] = DUx in UFOV
          output[1] = DUy in UFOV
          output[2] = coordinates max(DUx) in UFOV
          output[3] = coordinates max(DUy) in UFOV

          output[4] = DUx in CFOV
          output[5] = DUy in CFOV
          output[6] = coordinates max(DUx) in CFOV
          output[7] = coordinates max(DUy) in CFOV

          output[8] = IU in UFOV
          output[9] = IU in CFOV
          output[10] = ufov pixelmap (numpy masked array)
          output[11] = cfov pixelmap (numpy masked array)
        '''

        results.addFloat('%s_DU_x (UFOV)'%detector, output[0], 1, 'Differential Uniformity', '')
        results.addFloat('%s_DU_y (UFOV)'%detector, output[1], 1, 'Differential Uniformity', '')
        results.addFloat('%s_DU_x (CFOV)'%detector, output[4], 1, 'Differential Uniformity', '')
        results.addFloat('%s_DU_y (CFOV)'%detector, output[5], 1, 'Differential Uniformity', '')

        results.addFloat('%s_Integral Uniformity (UFOV)'%detector, output[8], 1, 'Integral Uniformity', '')
        results.addFloat('%s_Integral Uniformity (CFOV)'%detector, output[9], 1, 'Integral Uniformity', '')
        
        print "writing image results"

        try:
           filename = 'ufov_%s.png'%detector
           nemalib.save_imgmap(output[10],output[2],output[3],filename)
           results.addObject('UFOV image %s'%detector,filename)
        except:
           print "Could not save %s"%filename

        try:
           filename = 'cfov_%s.png'%detector
           nemalib.save_imgmap(output[11],output[6],output[7],filename)
           results.addObject('CFOV image %s'%detector,filename)
        except:
           print "Could not save %s"%filename

    
    
    
   


def nm_uniformity(data, results, **kwargs):

    print '--'*20
    print 'nm_uniformity function'

    uniformity_nm(data, results, **kwargs)

    print 'done calculating uniformity'

    print "\n\n"
    print "---->"
    print "Results from processing module:"

