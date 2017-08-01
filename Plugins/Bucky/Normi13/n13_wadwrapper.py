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
#   20170801: added mirror param
#   20170622: identify more float vars from header
#   20170310: add override params; take average over series
#   20161220: Removed class variables; removed testing stuff
#   20160825: fixes for portable detector
#   20160802: sync with wad2.0
#
#
from __future__ import print_function

__version__ = '20170801'
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
from . import n13_lib
import numpy as np

try: 
    # try local folder
    import wadwrapper_lib
except ImportError:
    # try pyWADlib from plugin.py.zip
    from pyWADLib import wadwrapper_lib

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

def logTag():
    return "[n13_wadwrapper] "

# helper functions
"""
    roomWKZ1 = n13_lib.Room('WKZ1', pid_tw=[70,50],
                               linepairmarkers={'type':'RXT02',xymm0.6':[-83.,-25.],'xymm1.0':[-99.,-8.]},artefactborderpx=[5,5,5,5],
                               detectorname={'SN152495':'Tafel', 'SN152508':'Wand', 'SN1522YG':'Klein1', 'SN151685':'Groot2'})
    <params>
      <roomname>WKZ1</roomname>
      <tablepidmm>70</tablepidmm>
      <wallpidmm>50</wallpidmm>
      <linepair_type>RXT02</linepair_type>
      <linepairs>
        <xymm1.0 x="-99.0" y="-8.0" />
        <xymm0.6 x="-83.0" y="-25.0" />
      </linepairs>
      <artefactborderpx>
        <xmin px="5" />
        <xmax px="5" />
        <ymin px="5" />
        <ymax px="5" />
      </artefactborderpx>
      <detector_name>
        <names SN151685="Groot2" SN1522YG="Klein1" SN152495="Tafel" SN152508="Wand" />
      </detector_name>
    </params>
"""
def override_settings(room, params):
    """
    Look for 'use_' params in to force behaviour of module and disable automatic determination of param.
    """
    try:
        room.pixmm = float(params.find('use_pixmm').text)
    except:
        pass
    try:
        room.mustbeinverted = params.find('use_mustbeinverted').text.lower() in ['1', 'true', 'y', 'yes']
    except:
        pass
    try:
        room.mustbemirrored = params.find('mustbemirrored').text.lower() in ['1', 'true', 'y', 'yes']
    except:
        pass

def _getRoomDefinition(params):
    # Use the params in the config file to construct an Scanner object
    try:
        # a name for identification
        roomname = params.find('roomname').text

        # Need to know the type of linepairs insert
        linepair_type = params.find('linepair_type').text
        if not linepair_type in ['None','RXT02','typ38']:
            raise ValueError('Incorrect linepair type %s'%linepair_type)
        
        # load the locations of markers on the linepair pattern. if these are not given, use the hardcoded values
        linepairmarkers = {}
        try:
            markers = params.find('linepairs')
            if linepair_type == 'RXT02':
                mnames = ['xymm1.0','xymm0.6']
            elif linepair_type == 'typ38':
                mnames = ['xymm1.8','xymm0.6','xymm1.4','xymm4.6']
                
            for mname in mnames:
                marker  = markers.find(mname)
                linepairmarkers[mname] = [ float(marker.attrib['x']), float(marker.attrib['y']) ]
        except:
            print(logTag()+' exact locations of markers on linepair pattern not supplied by config. Using empirical values; please check if these are valid here.')
        linepairmarkers['type'] = linepair_type
        
        # border to exclude
        artefactborderpx = [0,0,0,0]
        try:
            bpxs = params.find('artefactborderpx')
            mnames = ['xmin','xmax','ymin','ymax']
            for i,mname in enumerate(mnames):
                marker  = bpxs.find(mname)
                artefactborderpx[i] = int(marker.attrib['px'])
        except:
            print(logTag()+' no border supplied by config. Using [0,0,0,0].')
            
        # load detector names
        detectorname = {}
        try:
            dets_names = params.find('detector_name')
            names = dets_names.find('names')
            for a,b in names.items():
                detectorname[a] = b
        except:
            print(logTag()+' no explicit detector_name pairs defined in config.')

        # is a fixed setting forced?
        try:
            pidmm = [ float(params.find('pidmm').text) ]
        except:
            try:
                # Source to Detector distance and Patient to Detector distance for wall and table (both in mm)
                pidmm = [ float(params.find('tablepidmm').text), float(params.find('wallpidmm').text) ]
            except:
                try:
                    dummy = float(params.find('use_pixmm').text)
                    pidmm = [-1] # will fill in use_pixmm later
                except:
                    raise ValueError('Must supply "tablepidmm" and "wallpidmm", or "pidmm", or "use_pixmm"')

        try:
            sidmm = [ float(params.find('sidmm').text) ]
        except:
            try:
                sidmm = [ float(params.find('tablesidmm').text), float(params.find('wallsidmm').text) ]
            except:
                sidmm = [-1, -1] # not supplied

        # do we want a suffix added to the results, based on table/wall or detectorname?
        try:
            auto_suffix = params.find('auto_suffix').text.lower() in ['1', 'true', 'y', 'yes']
        except:
            auto_suffix = False
        print(logTag()+' auto_suffix set to ',auto_suffix)

        outvalue    = -1 # not supplied
        
        # no artificial thresholds present or needed
        room =  n13_lib.Room(roomname, outvalue=outvalue,
                               pid_tw=pidmm, sid_tw=sidmm,
                               artefactborderpx=artefactborderpx,
                               linepairmarkers=linepairmarkers,detectorname=detectorname,auto_suffix=auto_suffix)
        override_settings(room, params)
        return room

    except AttributeError as e:
        raise ValueError(logTag()+" missing room definition parameter!"+str(e))


###### Series wrappers
def xrayqc_series(data, results, params):
    """
    n13_UMCU checks:
        XRayEdges
        LowContrast
        DynamicRange
        MTF

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
    
    # select only middle slice for series
    numslices = len(pixeldataIn)
    if dicomMode == wadwrapper_lib.stMode3D:
        nim = int(len(pixeldataIn)/2.)
        dcmInfile   = dcmInfile._datasets[nim]
        pixeldataIn = np.average(pixeldataIn, axis=0)
        dicomMode = wadwrapper_lib.stMode2D

    ## 3. Build and populate qcstructure
    remark = ""
    qclib = n13_lib.XRayQC()
    room = _getRoomDefinition(params)
    cs = n13_lib.XRayStruct(dcmInfile,pixeldataIn,room)
    cs.verbose = False # do not produce detailed logging

    ## 4. Run tests
    error,msg = qclib.QC(cs)

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the WAD IQ database
    if not cs.DetectorSuffix() is None:
        idname = '_'+cs.DetectorSuffix()
    else:
        idname = ''

    ## first Build artefact picture thumbnail
    label = 'normi13'
    filename = '%s%s.jpg'%(label,idname) # Use jpg if a thumbnail is desired

    qclib.saveAnnotatedImage(cs,filename, 'normi13')
    results.addObject('%s%s'%(label,idname),filename)

    labvals = qclib.ReportEntries(cs)
    tmpdict={}
    for elem in labvals:
        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','pos':missing or a number} )
        # if no pos given, the next one will be given
        # if no quantity given, 'name' will be used
        # if no level given, the default will be used
        quan = elem['quantity'] if 'quantity' in elem else str(elem['name'])
        level = elem['level'] if 'level' in elem else None
        rank = elem['rank'] if 'rank' in elem else None
        results.addFloat(elem['name']+str(idname), elem['value'], quantity=quan, level=level,rank=rank)
    results.addFloat('num_slices'+str(idname), numslices)


def xrayqc_uniformity_series(data, results, params):
    """
    n13_uniformity checks:
        Uniformity
        Artefacts
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

    # select only middle slice for series
    numslices = len(pixeldataIn)
    if dicomMode == wadwrapper_lib.stMode3D:
        nim = int(len(pixeldataIn)/2.)
        dcmInfile   = dcmInfile._datasets[nim]
        pixeldataIn = np.average(pixeldataIn, axis=0)
        dicomMode = wadwrapper_lib.stMode2D

    ## 3. Build and populate qcstructure
    remark = ""
    qclib = n13_lib.XRayQC()
    room = _getRoomDefinition(params)
    cs = n13_lib.XRayStruct(dcmInfile,pixeldataIn,room)
    cs.verbose = False # do not produce detailed logging

    ## 4. Run tests
    error,msg = qclib.QCUnif(cs)

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the WAD IQ database
    if not cs.DetectorSuffix() is None:
        idname = '_'+cs.DetectorSuffix()
    else:
        idname = ''

    ## First Build artefact picture thumbnail
    label = 'unif'
    filename = '%s%s.jpg'%(label,idname) # Use jpg if a thumbnail is desired

    qclib.saveAnnotatedImage(cs,filename, 'artefacts')
    results.addObject('%s%s'%(label,idname),filename)

    labvals = qclib.ReportEntries(cs)
    tmpdict={}
    for elem in labvals:
        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','pos':missing or a number} )
        # if no pos given, the next one will be given
        # if no quantity given, 'name' will be used
        # if no level given, the default will be used
        quan = elem['quantity'] if 'quantity' in elem else str(elem['name'])
        level = elem['level'] if 'level' in elem else None
        rank = elem['rank'] if 'rank' in elem else None
        results.addFloat(elem['name']+str(idname), elem['value'], quantity=quan, level=level,rank=rank)
    results.addFloat('num_slices'+str(idname), numslices)


def xrayheader_series(data,results,params):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Read only headers
        2. Run tests
        3. Build xml output
    """
    info = 'qcwad'

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    ## 2. Run tests
    qclib = n13_lib.XRayQC()
    room = _getRoomDefinition(params)

    ## Table or Wall? from distances and sensitivity; for well defined protocols to be defined in DESCRIPTION field
    cs = n13_lib.XRayStruct(dcmInfile,None,room)
    cs.verbose = False # do not produce detailed logging
    dicominfo = qclib.DICOMInfo(cs,info)
    if not cs.DetectorSuffix() is None:
        idname = '_'+cs.DetectorSuffix()
    else:
        idname = ''

    ## 3. Build xml output
    floatlist = [
        'Exposure (mAs)',
        'Exposure (uAs)',
        'DistanceSourceToDetector (mm)',
        'ExposureTime (ms)',
        'ExposureTime (us)',
        'ImageAreaDoseProduct',
        'Sensitivity',
        'kVp',
        'CollimatorLeft',
        'CollimatorRight',
        'CollimatorUp',
        'CollimatorDown',
        'EntranceDose_mGy',
        'RelativeXRayExposure'
    ]
    offset = -26
    results.addChar('pluginversion'+idname, str(qclib.qcversion)) # do not specify level, use default from config
    for elem in dicominfo:
        quan = elem['quantity'] if 'quantity' in elem else str(elem['name'])
        level = elem['level'] if 'level' in elem else None # if not specify level, use default from config
        rank = elem['rank'] if 'rank' in elem else None
        if elem['name'] in floatlist:
            try:
                dummy = float(elem['value'])
            except ValueError:
                elem['value'] = -1
            results.addFloat(elem['name']+str(idname), float(elem['value']), quantity=quan, level=level,rank=rank)
        else:
            results.addChar(elem['name']+str(idname), str(elem['value'])[:min(len(str(elem['value'])),128)], quantity=quan, level=level,rank=rank)

    results.addChar('room'+idname, cs.forceRoom.name) # do not specify level, use default from config
    results.addChar('stand'+idname, cs.DetectorStand(), level=1,rank=offset+1) # do not specify level, use default from config
