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
#
# PyWAD is open-source software and consists of a set of modules written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes modules for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Changelog:
#   20160802: sync with wad2.0
#
#
from __future__ import print_function

__version__ = '20160802'
__author__ = 'aschilham'

import os
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import dicom
import QCCT_lib
import QCCT_constants as lit

try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib
import scipy.misc

def logTag():
    return "[QCCT_wadwrapper] "

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

# helper functions
def _getScannerDefinition(params):
    # Use the params in the config file to construct an Scanner object
    try:
        # a name for identification
        scannername = params.find('scannername').text

        # three materials in Philips Performance Phantom Head
        headHU_air   = float(params.find('headHU_air').text)
        headHU_water = float(params.find('headHU_water').text)
        headHU_pvc   = float(params.find('headHU_pvc').text)

        # inner and outer diameter (in mm) of PVC skull (container) of head phantom
        headdiammm_in    = float(params.find('headdiammm_in').text)
        headdiammm_out   = float(params.find('headdiammm_out').text)

        # three materials in Philips Performance Phantom Body
        bodyHU_aculon   = float(params.find('bodyHU_aculon').text)
        bodyHU_teflon   = float(params.find('bodyHU_teflon').text)
        try:
            bodyHU_water    = float(params.find('bodyHU_water').text)
        except:
            bodyHU_water    = float(params.find('bodyHU_air').text)
            
    except AttributeError as e:
        raise ValueError(logTag()+" missing scanner definition parameter!"+str(e))

    return QCCT_lib.Scanner(scannername, 
                            [headHU_air,headHU_water,headHU_pvc],
                            headdiammm_in, headdiammm_out,
                            [bodyHU_aculon,bodyHU_teflon,bodyHU_water])
    
##### Real functions
def ctqc_series(data,results,params):
    """
    QCCT_UMCU Checks: extension of Philips QuickIQ (also for older scanners without that option), for both Head and Body if provided
      Uniformity
      HU values
      Noise
      Linearity 

    Workflow:
        1. Read image or sequence
        2. Run test
        3. Build xml output
    """
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data.series_filelist[0],headers_only=False,logTag=logTag())
    qclib = QCCT_lib.CT_QC()
    cs = QCCT_lib.CTStruct(dcmInfile=dcmInfile,pixeldataIn=pixeldataIn,dicomMode=dicomMode)
    cs.forceScanner = _getScannerDefinition(params)
    cs.verbose = False

    ## id scanner
    error = qclib.DetermineCTID(cs)
    if(error == True or cs.guessScanner.name == lit.stUnknown):
        raise ValueError("{} ERROR! Cannot determine CT ID".format(logTag))

    error = qclib.HeadOrBody(cs)
    if(error == True or cs.anatomy == lit.stUnknown):
        raise ValueError("{} ERROR! Cannot determine Anatomy".format(logTag))

    idname = ""
    if cs.anatomy == lit.stHead:
        idname = "_H"
    else:
        idname = "_B"
    # only uncomment if same config used for all scanners: idname += cs.guessScanner.name

    ## 2. Run tests
    error = qclib.AnalyseCT(cs)
    if error:
        raise ValueError("{} ERROR! Error in AnalyseCT".format(logTag))

    ## Struct now contains all the results and we can write these to the
    ## WAD IQ database
    includedlist = [
        'roiavg',
        'roisd',
        'snr_hol',
        'unif',
        'linearity',
        'maxdev',
        'shiftxypx',
    ]
    excludedlist = [
        'verbose',
        'dcmInfile',
        'pixeldataIn',
        'dicomMode',
        'hasmadeplots',
        'guessScanner',
        'anatomy',
        'roiavg',
        'roisd',
        'snr_hol',
        'unif',
        'linearity',
        'maxdev',
        'shiftxypx',
        'valid'
        'skull_avg',
        'unif_slice',
        'unif_rois'
    ]
    for elem in cs.__dict__:
        if elem in includedlist:
            newkeys = []
            newvals = []
            try:
                elemval =  cs.__dict__[elem]
                if 'roiavg' in elem: # array of avgs
                    newkeys.append('MeanCenter')
                    newvals.append(elemval[0])
                    newkeys.append('MeanAir')
                    newvals.append(elemval[3])
                elif 'shiftxypx' in elem:
                    newkeys.append('shiftxpx')
                    newvals.append(elemval[0])
                    newkeys.append('shiftypx')
                    newvals.append(elemval[1])
                else:
                    newkeys.append(str(elem))
                    newvals.append(elemval)
            except:
                print(logTag()+"error for",elem)
                elemval = -1.
            for key,val in zip(newkeys,newvals):
                results.addFloat(key+str(idname), val, quantity=str(key))

    ## Build thumbnail
    filename = 'test'+idname+'.jpg' # Use jpg if a thumbnail is desired
    if cs.dicomMode == wadwrapper_lib.stMode2D:
        scipy.misc.imsave(filename,cs.pixeldataIn.transpose()) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!
    else:
        scipy.misc.imsave(filename,(cs.pixeldataIn[cs.unif_slice]).transpose()) # MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!
    results.addObject('CTslice'+idname,filename)

def ctheader_series(data,results,params):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Run tests
        2. Build xml output
    """

    info = 'dicom'
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data.series_filelist[0],headers_only=True,logTag=logTag())
    qcctlib = QCCT_lib.CT_QC()
    cs = QCCT_lib.CTStruct(dcmInfile=dcmInfile,pixeldataIn=pixeldataIn,dicomMode=dicomMode)
    cs.verbose = False

    error = qcctlib.HeadOrBody(cs)
    if(error == True or cs.anatomy == lit.stUnknown):
        raise ValueError("{} ERROR! Cannot determine Anatomy".format(logTag))

    result_dict = {}
    idname = ""
    if cs.anatomy == lit.stHead:
        idname = "_H"
    else:
        idname = "_B"
    # only uncomment if same config used for all scanners: idname += "_"+cs.guessScanner.name

    ## 1. Run tests
    dicominfo = qcctlib.DICOMInfo(cs,info)

    ## 2. Add results to 'result' object
    # plugionversion is newly added in for this plugin since pywad2
    results.addChar('pluginversion'+idname, str(qcctlib.qcversion)) # do not specify level, use default from config
    results.addChar('Anatomy', str(cs.anatomy)) # do not specify level, use default from config
    for di in dicominfo:
        results.addChar(di[0]+idname, str(di[1])[:min(len(str(di[1])),128)]) # do not specify level, use default from config


