import os
import ConfigParser
import dicom
import re
from datetime import datetime

import create_dicom
import storescu


def main(configFile):
    PROCESSED_FILEPATH = "processed_files.txt"
    open(PROCESSED_FILEPATH, 'a').close()
    
    #Read config
    config = ConfigParser.ConfigParser()
    config.read(configFile)

    sourcedir = config.get("FILECONFIG","SRCDIR")
    fileext = config.get("FILECONFIG","EXT")
    
    #Convert comma-separated values into list
    contains = [s.strip() for s in config.get("FILECONFIG","CONTAINS").split(',')]
    not_contains = [s.strip() for s in config.get("FILECONFIG","NOT_CONTAINS").split(',')]

    dcmconfig = {
        'patid':config.get("DCMCONFIG","PATID"),
        'patname':config.get("DCMCONFIG","PATNAME"),
        'studydes':config.get("DCMCONFIG","STUDYDES"),
        'seriesdes':config.get("DCMCONFIG","SERIESDES"),
        'stationname':config.get("DCMCONFIG","STATIONNAME"),}
    tag = dicom.tag.Tag(config.get("DCMCONFIG","TAG").split(','))

    aet =  config.get("SERVERCONFIG","AET")
    destip = config.get("SERVERCONFIG","IP")
    port = int(config.get("SERVERCONFIG","PORT"))
    
    #Processed filepaths
    with open(PROCESSED_FILEPATH, 'r') as f:
        processed = [f.strip() for f in f.readlines()]

    #Find all files in `sourcedir` with extension `fileext`
    allfiles = [fn for fn in os.listdir(sourcedir) if fn.endswith(fileext)]
    
    #Filter allfiles
    filepathlist = []
    for filename in allfiles:
        #Include filenames which contain all substrings in `contains`
        if not all([substr in filename for substr in contains]):
            continue
        #Exclude filenames which contain any substring in `not_contains`
        if any([substr in filename for substr in not_contains]):
            continue
        
        #Include if file hasn't been processed before
        path = os.path.join(sourcedir, filename)
        if path not in processed:
            filepathlist.append(os.path.join(sourcedir, filename))

    for path in filepathlist:
        print os.path.basename(path)
        
        #Get file content
        with open(path) as f:
            payload = ''.join(f.readlines())

        #Get date and time from filepath
        match = re.findall(r'\d+', path)
        dcmconfig['studydate'] = datetime.strptime(match[-2], '%d%m%Y')
        dcmconfig['studytime'] = datetime.strptime(match[-1], '%H%M%S')        

        #Create DCM
        tmpdicom = create_dicom.create_dicom(tag, payload, path, dcmconfig)

        #Send to PACS
        storescu.StoreSCU(aet, destip, port, tmpdicom)

        #Add path to processed_files.txt
        with open(PROCESSED_FILEPATH, 'a') as f:
            f.write("\n"+path)


if __name__ == "__main__":
    main("config.ini")
