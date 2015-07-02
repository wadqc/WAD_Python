import os,sys
import dicom
from dicom import tag
from dicom.dataset import Dataset, FileDataset
import dicom.UID
import argparse
from datetime import datetime
import StringIO

def extract_date(textfile):
    date = 0
    return date 


def create_dcm(filename,patname='BarcoQC',patid='112233'):


    textfile = open(filename,'rb')

    filename = os.path.splitext(filename)[0]+'.dcm'
    print "Output filename:", filename

    # Populate required values for file meta information
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.104.1' # NM Image Storage misschien beter OT
    file_meta.MediaStorageSOPInstanceUID = "1.2.3" # !! Need valid UID here for real work
    file_meta.ImplementationClassUID = "1.2.3.4" # !!! Need valid UIDs here
   
    print "Setting dataset values..."
   
    # Create the FileDataset instance (initially no data elements, but file_meta supplied)
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0"*128)
   
    # Add the data elements -- not trying to set all required here. Check DICOM standard
    ds.PatientName = patname
    ds.PatientID = patid
    
    #datestring = extract_date(textfile)
    #print datestring
    
    #ds.InstanceCreationDate = datestring
    # Set the transfer syntax
    ds.is_little_endian = True
    ds.is_implicit_VR = True

    #doctag = dicom.tag.Tag(("0042","0011"))
    doctag = dicom.tag.Tag(("0001","9999"))
    ds.add_new(doctag,'OB', textfile.read())

   
    print "Writing test file", filename
    ds.save_as(filename)
    print "File saved."
   


def main():
    parser = argparse.ArgumentParser(description='Commandline arguments to run reconstructor')  
    parser.add_argument('-f',metavar='file',nargs='?',type=str)

    args = vars(parser.parse_args())
    filename = args['f'] #dicom directory

    create_dcm(filename) 

    return


    


if __name__ == "__main__":
    sys.exit(main())
