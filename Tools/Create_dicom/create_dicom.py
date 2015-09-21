import dicom
from dicom.dataset import Dataset, FileDataset


def create_dicom(private_tag, payload, filename, dcmconfig):
    """ Function creates minimal dicom file from scratch with required tags
        and stores payload (string) in the specified private tag.
    """
    
    PatientID = dcmconfig['patid']
    PatientName = dcmconfig['patname']
    StudyDescription = dcmconfig['studydes']
    SeriesDescription = dcmconfig['seriesdes']
    StationName = dcmconfig['stationname']

    StudyDate = dcmconfig['studydate'].strftime('%Y%m%d')
    StudyTime = dcmconfig['studytime'].strftime('%H%M%S')

    # create empty dicomfile
    file_meta = Dataset()

    # Raw Data Storage
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.66'

    # unieke uid's
    file_meta.MediaStorageSOPInstanceUID = dicom.UID.generate_uid()
    file_meta.ImplementationClassUID = dicom.UID.generate_uid()

    ds = FileDataset(filename, {},file_meta = file_meta,preamble="\0"*128)

    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.7' # secondary capture SOP UID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.StudyInstanceUID = dicom.UID.generate_uid()
    ds.SeriesInstanceUID = dicom.UID.generate_uid()

    ds.PatientID = PatientID
    ds.PatientName = PatientName
    ds.StudyDescription = StudyDescription
    ds.SeriesDescription = SeriesDescription
    ds.StationName = StationName
    ds.Modality = 'OT'

    ds.StudyDate =  StudyDate
    ds.SeriesDate =  ds.StudyDate

    ds.ContentDate = ds.StudyDate
    ds.StudyTime = ds.SeriesTime = ds.ContentTime = StudyTime

    ds.add_new(private_tag,'OB', payload)

    return ds
