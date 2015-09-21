'''
This is a tool to encapsulate a single file in a dicom file and optionally
send it to the WAD server
'''

__author__='DD'
__version__='21092015'

import os,path
import sys
import ConfigParser
import dicom
from datetime import datetime
import create_dicom
import storescu


def main(configFile,filename):
    
    #Read config
    config = ConfigParser.ConfigParser()
    config.read(configFile)

    dcmconfig = {
        'patid':config.get("DCMCONFIG","PATID"),
        'patname':config.get("DCMCONFIG","PATNAME"),
        'studydes':config.get("DCMCONFIG","STUDYDES"),
        'seriesdes':config.get("DCMCONFIG","SERIESDES"),
        'stationname':config.get("DCMCONFIG","STATIONNAME"),}
    tag = dicom.tag.Tag(config.get("DCMCONFIG","TAG").split(','))

    dcmconfig['studydate'] = datetime.strptime(config.get("DCMCONFIG","STUDYDATE"), '%Y%m%d')
    dcmconfig['studytime'] = datetime.strptime(config.get("DCMCONFIG","STUDYTIME"), '%H%M%S')

    aet =  config.get("SERVERCONFIG","AET")
    destip = config.get("SERVERCONFIG","IP")
    port = int(config.get("SERVERCONFIG","PORT"))
    
    #Get file content
    with open(filename) as f:
        payload = ''.join(f.readlines())


    #Create DCM
    tmpdicom = create_dicom.create_dicom(tag, payload, filename, dcmconfig)

    #Send to PACS
    if config.get("SERVERCONFIG","IP") == "": 
        print("No IP specified in config.ini; saving dicom file locally")
        savename = filename+'.dcm'
        tmpdicom.save_as(savename)

    else:
        try:
            storescu.StoreSCU(aet, destip, port, tmpdicom)
        except:
            print('Failed to send file to DCM4CHEE, check server settings?')
        
if __name__ == "__main__":
    
    main("config.ini",sys.argv[1])
