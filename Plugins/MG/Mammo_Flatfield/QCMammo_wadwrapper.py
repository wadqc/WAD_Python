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

try:
    import pydicom as dicom
except ImportError:
    import dicom
from . import QCMammo_lib
try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib

def logTag():
    return "[QCMammo_wadwrapper] "

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!


##### Series wrappers
def mammoqc_series(data, results, **kwargs):
    """
    QCMammo_UMCU checks:
        Uniformity (5 rois) and SNR (hologic),
        DoseRatio (empirical/calculated from DICOM)
        Artefacts (spots, dead pixels)

    Workflow:
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
        6. Build artefact picture thumbnail
    """
    inputfile = data.series_filelist[0]  # give me a filename

    ## 2. Check data format
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(inputfile,headers_only=False,logTag=logTag())

    ## 3. Build and populate qcstructure
    remark = ""
    qcmammolib = QCMammo_lib.Mammo_QC()
    cs_mam = QCMammo_lib.MammoStruct(dcmInfile,pixeldataIn)
    cs_mam.verbose = False # do not produce detailed logging
    if qcmammolib.NeedsCropping(cs_mam):
        cs_mam.expertmode = True
        qcmammolib.RestrictROI(cs_mam)
        remark = "CROPPED"

    ## 4. Run tests
    # Uniformity check
    error = qcmammolib.Uniformity(cs_mam)
    if error:
        remark += "/ERROR_UNIFORMITY"
    # Contrast L50
    error = qcmammolib.L50Contrast(cs_mam) # doesn't do anything if not L50
    if error:
        remark += "/ERROR_L50CONTRAST"
    # Dose Ratio
    error = qcmammolib.DoseRatio(cs_mam)
    if error:
        remark += "/ERROR_DOSERATIO"
    # Artefacts
    error = qcmammolib.Artefacts(cs_mam)
    if error:
        remark += "/ERROR_ARTEFACTS"

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the
    ## WAD IQ database
    includedlist = [
        'unif_pct',
        'snr_hol',
        'doseratio',
        'art_clusters',
        'expert_inoutoverin',
        'contrast_snr'
    ]
    excludedlist = [
        'verbose',
        'dcmInfile',
        'pixeldataIn',
        'hasmadeplots',
        'means',
        'stdevs',
        'unif',
        'snr_hol',
        'unif_rois',
        'doseratio',
        'art_clusters',
        'art_image',
        'art_borderpx',
        'art_threshold',
        'art_rois',
        'expertmode',
        'expert_roipts',
        'expert_frac',
        'expert_inoutoverin',
        'filtername',
        'scannername',
        'contrast_rois',
        'contrast_mean',
        'contrast_sd'
    ]

    idname = '_'+cs_mam.filtername
    if "SUMPRES" in idname:
        idname = "" # only one image

    results.addChar('NOTE'+idname, remark)

    for elem in cs_mam.__dict__:
        if elem in includedlist:
            newkeys = []
            newvals = []
            try:
                elemval =  cs_mam.__dict__[elem]
                if 'contrast_snr' in elem: # array
                    for ix,snr in enumerate(elemval):
                        newkeys.append('CNR'+str(ix))
                        newvals.append(snr)
                elif 'art_clusters' in elem:
                    newkeys.append(str(elem))
                    newvals.append(len(elemval))
                else:
                    newkeys.append(str(elem))
                    newvals.append(elemval)
            except:
                print(logTag()+"error for",elem)

            tmpdict={}
            for key,val in zip(newkeys,newvals):
                results.addFloat(key+str(idname), val, quantity=str(key))


    ## 6. Build artefact picture thumbnail
    filename = 'test'+idname+'.jpg' # Use jpg if a thumbnail is desired

    #object_naam_pad = outputfile.replace('result.xml','test'+idname+'.jpg') # Use jpg if a thumbnail is desired
    qcmammolib.saveAnnotatedArtefactImage(cs_mam,filename)
    results.addObject('ArtefactImage'+idname,filename)

def mammoheader_series(data,results,params):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Run tests
        2. Build xml output
    """

    try:
        info = params.find("info").text
    except AttributeError:
        info = 'qc' # selected subset of DICOM headers informative for QC testing

    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    ## 1. Run tests
    qcmammolib = QCMammo_lib.Mammo_QC()
    cs = QCMammo_lib.MammoStruct(dcmInfile,None)
    cs.verbose = False # do not produce detailed logging
    dicominfo = qcmammolib.DICOMInfo(cs,info)

    ## find filtername
    idname = "_XX"
    for di in dicominfo:
        if "FilterMaterial" in di[0]:
            if "RHO" in di[1]:
                idname = "_RH"
                break
            elif "MOL" in di[1]:
                idname = "_MO"
                break
            elif "ALU" in di[1]:
                idname = "" # only one image "_SUMPRES"
                break
            elif "SILV" in di[1]:
                idname = "_AG"

    ## 2. Add results to 'result' object
    results.addChar('pluginversion'+idname, str(qcmammolib.qcversion)) # do not specify level, use default from config
    for di in dicominfo:
        results.addChar(di[0]+idname, str(di[1])[:min(len(str(di[1])),128)]) # do not specify level, use default from config

