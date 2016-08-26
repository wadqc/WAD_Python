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

__version__ = '20150826'
__author__ = 'aschilham'


import sys
import numpy as np
import os
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import CDMam_lib
try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib

def logTag():
    return "[CDMam_wadwrapper] "

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!


##### Series wrappers
def cdmamsetup_series(inputfile,params,headers_only):
    """
    Shared routine to set runtime parameters and build structure
    Workflow:
        1. Set runtime parameters
        2. Check data format
        3. Build and populate qcstructure
    """
    # 1. Set runtime parameters
    if not headers_only: # then we can skip the parameter check
        try:
            phantomversion = params.find("phantomversion").text
        except AttributeError:
            raise ValueError(logTag()+" missing phantomversion parameter!")
    
        try:
            modeCDCOM = (params.find("modeCDCOM").text == 'True')
        except AttributeError:
            raise ValueError(logTag()+" missing cdcommode parameter!")
    else:
        phantomversion = '3.2' # dummy for headers
        modeCDCOM = False

    # 2. Check data format
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(inputfile,headers_only=headers_only,logTag=logTag())

    # 3. Build and populate qcstructure
    qclib = CDMam_lib.CDMam()
    cs = CDMam_lib.CDMamStruct(dcmInfile,pixeldataIn,phantomversion)
    cs.imnum = 0 # needs to be set to an increasing number if multiple files are used
    cs.verbose = False # do not produce detailed logging
    cs.imageFileName = inputfile[0] # only for CDCOM.exe
    qclib.determineScannerID(cs)
    return qclib,cs,modeCDCOM

def identifyName(cs):
    """
    Come up with a simple identification tag
    """
    idname = "_XX"
    if "RHO" in cs.filtername:
        idname = "_RH"
    elif "MOL" in cs.filtername:
        idname = "_MO"
    elif "ALU" in cs.filtername:
        idname = "_AL" # only one image "_SUMPRES"
    elif "SILV" in cs.filtername:
        idname = "_AG"

    if 'SUM' in cs.energypresentation: # sum, high, low
        idname += 'su'
    elif 'HIGH' in cs.energypresentation: # sum, high, low
        idname += 'hi'
    elif 'LOW' in cs.energypresentation: # sum, high, low
        idname += 'lo'
    if 'PROC' in cs.energypresentation: # presentation or processing
        idname += 'proc'
    elif 'PRES' in cs.energypresentation: # presentation or processing
        idname += 'pres'
    
    return idname

def cdmamqc_list(data,results,params):
    """
    CDMAM analysis for a list of images (as should be!) see cdmamqc_series for more information
    """
    
    finalscore = (np.zeros(np.shape(cs.phantom.groundtruth),dtype=float)).tolist()
    # 1.-3.
    imnum = 0
    for fname in data.series_filelist[0]:
        inputfile = [fname]
        qclib,cs,modeCDCOM = cdmamsetup_series(inputfile, params, headers_only=False)
        cs.imnum = imnum

        if not modeCDCOM:
            score = qclib.CDMamSingle(cs)
        else:
            score = qclib.CDCOMSingle(cs)
        for y in range(len(score)):
            for x in range(len(score[0])):
                finalscore[y][x] += score[y][x]
        
        imnum +=1
        
    # calculate average score
    for y in range(len(finalscore)):
        for x in range(len(finalscore[0])):
            finalscore[y][x] /= len(imagelist)

    qclib.thresholdThickness(cs,finalscore)

def cdmamqc_series(data, results, params):
    """
    CDMAM analysis:
        Use a model observer to detect location of gold disc, for different thicknesses and diameters of the disc.
        Predict the threshold thicknesses for each diameter for a human observer
        Calculate a comprehensive image quality factor from the predicted human observer performance.
        NOTE: CDMAM analysis should be performed on 10-16 slightly shifted images to get a meaningful number.
              Here the analysis is done on one image, but the interface for multiple images already exist as cdmamqc_list.
    Params needs to define:
      1. phantomversion = phantomversion (which version of the phantom is used? '3.2' or '3.4')
      2. modeCDCOM = 'False' or 'True' (use external executable cdcom.exe ('True') or not ('False'))

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
    qclib,cs,modeCDCOM = cdmamsetup_series(inputfile, params, headers_only=False)

    # 4. Run tests
    if not modeCDCOM:
        score = qclib.CDMamSingle(cs)
    else:
        score = qclib.CDCOMSingle(cs)
    qclib.thresholdThickness(cs,score)

    # 5. Build xml output
    # Struct now contains all the results and we can write these to the WAD IQ database
    includedlist = [
        'diam_mm', # list of diameters in mm
        'limit_um', # list of predicted human detectable thickness thresholds in um
        'iqf', # image quality factor
        'threshold_fnames',
        'fit_fnames'
    ]
    idname = identifyName(cs)
    for di,lim in zip(cs.diam_mm,cs.limit_um):
        results.addFloat('limit_um_'+str('%0.2f'%di).zfill(4)+idname, lim, quantity=str('threshold_limit'))
    results.addFloat('iqf'+idname, cs.iqf, quantity=str('iqf'))

    # 6. also store images as results
    for fn in cs.threshold_fnames:
        results.addObject(os.path.splitext(fn)[0],fn)
    for fn in cs.fit_fnames:
        results.addObject(os.path.splitext(fn)[0],fn)

def cdmamheader_series(data,results,params):
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
    qclib,cs,modeCDCOM = cdmamsetup_series(inputfile, params, headers_only=True)

    # 4. Run tests
    dicominfo = qclib.DICOMInfo(cs,'qc')

    ## find filtername
    idname = identifyName(cs)
        
    ## 2. Add results to 'result' object
    results.addChar('pluginversion'+idname, str(qclib.qcversion)) # do not specify level, use default from config
    for di in dicominfo:
        results.addChar(di[0]+idname, str(di[1])[:min(len(str(di[1])),128)]) # do not specify level, use default from config

