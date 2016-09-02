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
#   20160825: added extra config parameters (PvH)
#   20160802: sync with pywad1.0
from __future__ import print_function

__version__ = '20160825'
__author__ = 'aschilham'

import os
import numpy as np
if not 'MPLCONFIGDIR' in os.environ:
    # using a fixed folder is preferable to a tempdir, because tempdirs are not automatically removed
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import QCUS_lib
import ocr_lib
try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib

def logTag():
    return "[QCUS_wadwrapper] "

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

##### Series wrappers
def setup_series(inputfile,params,headers_only):
    """
    Shared routine to set runtime parameters and build structure
    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
    """
    # 1. Set runtime parameters
    # 2. Check data format
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(inputfile,headers_only=headers_only,logTag=logTag())

    # 3. Build and populate qcstructure
    qclib = QCUS_lib.US_QC(guimode=False)
    cs = QCUS_lib.USStruct(dcmInfile,pixeldataIn,dicomMode)
    cs.verbose = False # do not produce detailed logging
    cs.uni_filter = int(params.find("uni_filter").text)
    cs.uni_delta = float(params.find("uni_delta").text)
    cs.uni_low = float(params.find("uni_low").text)
    cs.sen_filter = int(params.find("sen_filter").text)
    cs.sen_delta = float(params.find("sen_delta").text)
    cs.ver_offset = int(params.find("ver_offset").text)
    cs.hor_offset = int(params.find("hor_offset").text)
    cs.fitcircle_frac = float(params.find("fitcircle_frac").text)
    cs.cluster_fminsize = float(params.find("cluster_fminsize").text)
    try:
        cs.verbose = params.find('verbose').text.lower() in ['1', 'true', 'y', 'yes']
    except:
        pass
    return qclib,cs

def qc_series(data, results, params):
    """
    US Reverberations in Air analysis:
        Check the uniformity of the reverberation patterns

    Params needs to define:
       nothing yet
       
    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
        6. Build artefact picture thumbnail
    """

    # 1.-3.
    inputfile = data.series_filelist[0]  # give me a [filename]
    qclib,cs = setup_series(inputfile, params, headers_only=False)

    # 4. Run tests
    error = qclib.Analyse(cs)

    # 5. Build xml output
    idname = '_'+qclib.imageID(cs,probeonly=True)
    labvals = qclib.reportEntries(cs)
    #    image_fnames = [] # filenames of generated images

    for key,val in labvals:
        results.addFloat(key+str(idname), val, quantity=str(key))

    # 6. also store images as results
    for fn in cs.image_fnames:
        results.addObject(os.path.splitext(fn)[0],fn)

    # also run ocr_series; needed as part of qc because of the boxes it generates
    ocr_rois, error, msg = ocr_series(data, results, params)
    xtra= {'rectrois': ocr_rois }
    fname = 'overview.jpg'
    qclib.saveAnnotatedImage(cs, fname, what='overview',xtra=xtra)
    results.addObject(os.path.splitext(fname)[0],fname)

    fname = 'reverb.jpg'
    qclib.saveAnnotatedImage(cs, fname, what='reverb')
    results.addObject(os.path.splitext(fname)[0],fname)
    if error:
        raise ValueError('%s Cannot read OCR box for %s'%(logTag(),msg))

def ocr_series(data, results, params):
    """
    Use pyOCR which for OCR
    returns rect rois for plotting in overview
    """
    inputfile = data.series_filelist[0]  # give me a [filename]
    dcmInfile, pixeldataIn, dicomMode = wadwrapper_lib.prepareInput(inputfile, headers_only=False, logTag=logTag())

    rectrois = []
    error = False
    msg = ''

    # solve ocr params
    regions = {}
    for param in params:
        #'OCR_TissueIndex.xywh' = 'x;y;w;h'
        #'OCR_TissueIndex.prefix' = 'prefix'
        #'OCR_TissueIndex.suffix' = 'suffix'
        if param.tag.startswith('OCR_'):
            split = param.tag.find('.') # ':' is illegal in xml names
            name = param.tag[:split]
            stuff = param.tag[split+1:]
            if not name in regions:
                regions[name] = {'prefix':'', 'suffix':''}
            if stuff == 'xywh':
                regions[name]['xywh'] = [int(p) for p in param.text.split(';')]
            elif stuff == 'prefix':
                regions[name]['prefix'] = param.text
            elif stuff == 'suffix':
                regions[name]['suffix'] = param.text
            elif stuff == 'type':
                regions[name]['type'] = param.text

    for name, region in regions.items():
        rectrois.append([ (region['xywh'][0],region['xywh'][1]), 
                          (region['xywh'][0]+region['xywh'][2],region['xywh'][1]+region['xywh'][3])])

        txt, part = ocr_lib.OCR(pixeldataIn, region['xywh'])
        if region['type'] == 'object':
            import scipy
            im = scipy.misc.toimage(part) 
            fn = '%s.jpg'%name
            im.save(fn)
            results.addObject(name, fn)
            
        else:
            try:
                value = ocr_lib.txt2type(txt, region['type'], region['prefix'],region['suffix'])
                if region['type'] == 'float':
                    results.addFloat(name, value)
                elif region['type'] == 'string':
                    results.addChar(name, value)
                elif region['type'] == 'bool':
                    results.addBool(name, value)
            except:
                error = True
                msg += name + ' '
    return rectrois, error, msg

def header_series(data,results,params):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
    """
    # 1.-3.
    inputfile = data.series_filelist[0]  # give me a [filename]
    qclib,cs = setup_series(inputfile, params, headers_only=True)

    # 4. Run tests
    dicominfo = qclib.DICOMInfo(cs)

    ## find probename
    idname = '_'+qclib.imageID(cs,probeonly=True)
        
    ## 2. Add results to 'result' object
    results.addChar('pluginversion'+idname, str(qclib.qcversion)) # do not specify level, use default from config
    for di in dicominfo:
        results.addChar(di[0]+idname, str(di[1])[:min(len(str(di[1])),128)]) # do not specify level, use default from config

