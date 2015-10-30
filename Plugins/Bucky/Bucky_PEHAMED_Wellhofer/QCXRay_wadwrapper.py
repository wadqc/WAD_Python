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
# 
#


__version__='20151027'
__author__ = 'aschilha'
import sys
import os
if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

import dicom
import QCXRay_lib
try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

def logTag():
    return "[QCXRay_wadwrapper] "

# helper functions
"""
    roomWKZ1 = Room("WKZ1",outvalue=1023,tablesid=1150,wallsid=2000, tablepid=65, wallpid=50,phantom=lit.stWellhofer)
    <params>
      <roomname>WKZ1</roomname>
      <phantom>wellhofer</phantom>
      <tablesidmm>1150</tablesidmm>
      <tablepidmm>65</tablepidmm>
      <wallsidmm>2000</wallsidmm>
      <wallpidmm>50</wallpidmm>
      <outvalue>1023</outvalue>
      
      <sensitivities>
        <threshold date="20100101" value="35" />
      </sensitivities>
      
      <sdthreshold>40</sdthreshold>
    </params>
"""
def _getRoomDefinition(params):
    # Use the params in the config file to construct an Scanner object
    try:
        # a name for identification
        roomname = params.find('roomname').text

        # phantom name (only pehamed or wellhofer)
        phantoms_supported = ['pehamed','wellhofer']
        phantom = params.find('phantom').text
        if not phantom in phantoms_supported:
            raise ValueError(logTag()+' unsupported phantom %s'%phantom)

        # load the locations of markers on the linepair pattern. if these are not given, use the hardcoded values
        linepairmarkers = {}
        try:
            markers = params.find('linepair_typ38')
            mnames = ['mm1.8','mm0.6','mm1.4','mm4.6']
            for mname in mnames:
                marker  = markers.find(mname)
                linepairmarkers[mname] = [ float(marker.attrib['x']), float(marker.attrib['y']) ]
        except:
            print logTag()+' exact locations of markers on linepair pattern not supplied by config. Using empirical values; please check if these are valid here.'
            
        # Source to Detector distance and Patient to Detector distance for wall and table (both in mm)
        tablepidmm  = float(params.find('tablepidmm').text)
        wallpidmm   = float(params.find('wallpidmm').text)

        outvalue    = -1 # not supplied
        wallsidmm   = -1 # not supplied
        tablesidmm  = -1 # not supplied
        try: # only for FCR
            wallsidmm   = float(params.find('wallsidmm').text)
            tablesidmm  = float(params.find('tablesidmm').text)
            # pixelvalue that defines 'outside phantom' use '-1' to calculate from four cornerpoints
            outvalue    = int(params.find('outvalue').text)
        except:
            pass
        
        
        # for fcr systems there is no dicom tag to indicate wall or table, but a hack on SD or Sensitivity is possible
        try:
            thresholdlist = []
            sensitivities = params.find("sensitivities")
            for threshold in sensitivities.findall("threshold"):
                thresholdlist.append([int(threshold.attrib["date"]),int(threshold.attrib["value"])])
            return QCXRay_lib.Room(roomname, outvalue=outvalue,
                                   tablesid=tablesidmm, wallsid=wallsidmm, 
                                   tablepid=tablepidmm, wallpid=wallpidmm,
                                   phantom=phantom, sens_threshold = thresholdlist,
                                   linepairmarkers=linepairmarkers)
        except:
            pass

        # no sensitivity threshold, so try if threshOnSD exists
        try:
            sdthreshold = float(params.find("sdthreshold").text)
            return QCXRay_lib.Room(roomname, outvalue=outvalue,
                                   tablesid=tablesidmm, wallsid=wallsidmm, 
                                   tablepid=tablepidmm, wallpid=wallpidmm,
                                   phantom=phantom, sdthresh = sdthreshold,
                                   linepairmarkers=linepairmarkers)
        except:
            pass

        # no artificial thresholds present or needed
        return QCXRay_lib.Room(roomname, outvalue=outvalue,
                               tablesid=tablesidmm, wallsid=wallsidmm, 
                               tablepid=tablepidmm, wallpid=wallpidmm,
                               phantom=phantom,linepairmarkers=linepairmarkers)
    except AttributeError,e:
        raise ValueError(logTag()+" missing room definition parameter!"+str(e))


###### Series wrappers
def xrayqc_series(data, results, params):
    """
    QCXRay_UMCU checks:
        Horizontal uniformity
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

    ## 3. Build and populate qcstructure
    remark = ""
    qclib = QCXRay_lib.XRayQC()
    room = _getRoomDefinition(params)
    cs = QCXRay_lib.XRayStruct(dcmInfile,pixeldataIn,room)
    cs.verbose = False # do not produce detailed logging

    ## 4. Run tests
    error,msg = qclib.QC(cs)

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the WAD IQ database
    stand = qclib.TableOrWall(cs)
    idname = '_'+stand

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

    ## 6. Build artefact picture thumbnail
    filename = 'test'+idname+'.jpg' # Use jpg if a thumbnail is desired

    qclib.saveAnnotatedImage(cs,filename)
    results.addObject('AnnotatedImage'+idname,filename)

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
    qclib = QCXRay_lib.XRayQC()
    room = _getRoomDefinition(params)

    ## Table or Wall? from distances and sensitivity; for well defined protocols to be defined in DESCRIPTION field
    cs = QCXRay_lib.XRayStruct(dcmInfile,None,room)
    cs.verbose = False # do not produce detailed logging
    dicominfo = qclib.DICOMInfo(cs,info)
    idname = '_'+qclib.TableOrWall(cs)

    ## 3. Build xml output
    floatlist = [
        'Exposure (mAs)',
        'DistanceSourceToDetector (mm)',
        'ExposureTime (ms)',
        'ImageAreaDoseProduct',
        'Sensitivity',
        'kVp'
    ]
    offset = -25
    results.addChar('pluginversion'+idname, str(qclib.qcversion)) # do not specify level, use default from config
    for elem in dicominfo:
        quan = elem['quantity'] if 'quantity' in elem else str(elem['name'])
        level = elem['level'] if 'level' in elem else None # if not specify level, use default from config
        rank = elem['rank'] if 'rank' in elem else None
        if elem['name'] in floatlist:
            results.addFloat(elem['name']+str(idname), elem['value'], quantity=quan, level=level,rank=rank)
        else:
            results.addChar(elem['name']+str(idname), str(elem['value'])[:min(len(str(elem['value'])),128)], quantity=quan, level=level,rank=rank)

    results.addChar('room'+idname, cs.forceRoom.name) # do not specify level, use default from config
    results.addChar('stand'+idname, qclib.TableOrWall(cs), level=1,rank=offset+1) # do not specify level, use default from config
    results.addBool('boolshit'+idname, True) # do not specify level, use default from config
