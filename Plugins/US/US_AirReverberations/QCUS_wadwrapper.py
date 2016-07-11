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

__version__ = '01062015'
__author__ = 'aschilham'




import sys
import numpy as np
import os
if not 'MPLCONFIGDIR' in os.environ:
    # using a fixed folder is preferable to a tempdir, because tempdirs are not automatically removed
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import QCUS_lib
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
    """
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
    """
    # 2. Check data format
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(inputfile,headers_only=headers_only,logTag=logTag())

    # 3. Build and populate qcstructure
    qclib = QCUS_lib.US_QC(guimode=False)
    cs = QCUS_lib.USStruct(dcmInfile,pixeldataIn,dicomMode)
    cs.verbose = False # do not produce detailed logging
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

