# functions to send strings as private tags in a dummy dicom file to the server
#
import socket
import uuid
from netdicom import AE, StorageSOPClass, VerificationSOPClass
try:
    import pydicom as dicom
except ImportError:
    import dicom
import os

def sendasdicom(dcmconfig,datastr):
    """
    create a dicom file with the datastr embedded as a private tag, and send it to dcm4chee
      dcmconfig is a dict that should at least contain: dest_aet, src_aet, dest_ip, dest_port
      datastr is a string that should be of even string length (required for private tags content)

    returns output of create_dicom: (Succes,msg) with Succes True or False and msg a string or None
    """
    # fill in defaults if these fields are missing
    if not 'dest_aet' in dcmconfig: dcmconfig['dest_aet'] = 'DCM4CHEE' # WAD-IQC PACS AE TITLE
    if not 'src_aet' in dcmconfig: dcmconfig['src_aet'] = 'PYNETDICOM' # AE TITLE TO STORE IN SRC_AETITLE
    if not 'dest_ip' in dcmconfig: dcmconfig['dest_ip'] = 'localhost' # WAD-IQC PACS IP
    if not 'dest_port' in dcmconfig: dcmconfig['dest_port'] = 11112 # WAD-IQC PACS PORT (no quotes!)
    
    # space padding until string length is even (required for private tags content)
    payload_str  = str(datastr)
    payload_size = len(payload_str)
    if payload_size%2 != 0:
        payload_str = payload_str.ljust(payload_size+1)
    
    ##Create DCM
    # generate unique temporary filename for dicom output
    unique_dcm_filename = '/tmp/'+str(uuid.uuid4())+'.dcm'

    # create dicom file and inject payload
    tmpdicom = create_dicom(dcmconfig,payload_str,unique_dcm_filename)

    #Send to PACS
    return send_dicom(dcmconfig['dest_aet'],dcmconfig['src_aet'],dcmconfig['dest_ip'],dcmconfig['dest_port'],unique_dcm_filename)

def create_dicom(dcmconfig,payload,filename):
    """ Function creates minimal dicom file from scratch with required tags
        and stores payload (string) in the specified private tag.
    """
    ptag = dicom.tag.Tag(dcmconfig['privatetag'].split(','))

    # create empty dicomfile
    file_meta = dicom.dataset.Dataset()

    # Raw Data Storage
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.66'

    # unieke uid's
    file_meta.MediaStorageSOPInstanceUID = dicom.UID.generate_uid()
    file_meta.ImplementationClassUID = dicom.UID.generate_uid()

    ds = dicom.dataset.FileDataset(filename, {},file_meta = file_meta,preamble="\0"*128)

    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.7' # secondary capture SOP UID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.StudyInstanceUID = dicom.UID.generate_uid()
    ds.SeriesInstanceUID = dicom.UID.generate_uid()

    ds.PatientID = dcmconfig['patientid']
    ds.PatientName = dcmconfig['patientname']
    ds.StudyDescription = dcmconfig['studydesc']
    ds.SeriesDescription = dcmconfig['seriesdesc']
    ds.StationName = dcmconfig['stationname']
    ds.Modality = 'OT'
    ds.StudyDate   = dcmconfig['studydate']
    ds.SeriesDate  = dcmconfig['seriesdate']
    ds.ContentDate = dcmconfig['contentdate']
    ds.StudyTime   = dcmconfig['studytime']
    ds.SeriesTime  = dcmconfig['seriestime']
    ds.ContentTime = dcmconfig['contenttime']

    ds.add_new(ptag,'OB', payload)

    ds.save_as(filename)

def send_dicom(calledAET,callingAET,remotehost,remoteport,filename):
    """ Function sends a dicom file to specified host.
        Returns (status, msg):
          Status:
            True  : success
            False : error occurred (added to logfile)
          msg: string or None
    """
    verbose = True

    ts = [
        dicom.UID.ExplicitVRLittleEndian, 
        dicom.UID.ImplicitVRLittleEndian, 
        dicom.UID.ExplicitVRBigEndian
    ]

    # call back
    def OnAssociateResponse(association):
        if verbose: print "[Input2DCM] Association response received"
        pass

    # find random free port on localhost
    tmpsocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    tmpsocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
    tmpsocket.bind(('', 0))
    freelocalport = tmpsocket.getsockname()[1]
    tmpsocket.close()

    # create application entity
    MyAE = AE(callingAET, freelocalport, [StorageSOPClass,  VerificationSOPClass], [], ts)
    MyAE.OnAssociateResponse = OnAssociateResponse

    # remote application entity
    RemoteAE = dict(Address=remotehost, Port=remoteport, AET=calledAET)

    # create association with remote AE
    if verbose: print '[Input2DCM] Request association'
    assoc = MyAE.RequestAssociation(RemoteAE)

    if not assoc:
        msg = 'Could not establish association'
        print '[Input2DCM] %s'%msg
        return (False,msg)

    # perform a DICOM ECHO, just to make sure remote AE is listening
    try:
        status = assoc.VerificationSOPClass.SCU(1)
        if verbose: print '[Input2DCM] DICOM Echo ... done with status "%s"' % status
    except:
        msg = 'DICOM Echo ... problem occurred'
        print '[Input2DCM] %s'%msg
        return (False,msg)

    dcm = dicom.read_file(filename)
    try:
        status = assoc.SCU(dcm, 1)
        if verbose: print '[Input2DCM] DICOM StoreSCU ... done with status "%s"' % status
    except:
        msg = 'DICOM StoreSCU ... problem occurred'
        print '[Input2DCM] %s'%msg
        return (False,msg)

    if verbose: print 'Release association'
    assoc.Release(0)
    os.remove(filename)

    # done
    MyAE.Quit()

    return (True,None)

def write_generic_config(config,filename):
    conf = dict(config) # make a copy for usage
    defaults = {
        'info_name':'Input2DCM',
        'info_description':'Input2DCM Generic Config',
        'info_version': '20150814',
        'actions': [
            {
                'plugin': 'Plugin_development.Generic.generic',
                'function': 'main',
                'default_level': '1',
                'params': {
                    'tag_type':'dict',
                    'private_group':'0x0071',
                    'private_element':'0x9999'
                },
                'limits': {
                    'sensitivity': {
                        'acc_low':'1.23',
                        'acc_high':'2.5',
                        'crit_low':'1.00',
                        'crit_high':'3.00',
                    }
                }
            }
        ]
    }
  
    for key,value in defaults.items():
        if not key in conf: conf[key] = value

    with open(filename,'w') as fout:
        fout.write('<?xml version="1.0" encoding="UTF-8"?>\n')

        fout.write('<!--\n')
        fout.write('  Configuration file written by Input2DCM.py\n')
        fout.write('-->\n')

        fout.write('<Input2DCM_config>\n')

        fout.write('  <name>%s</name>\n'% conf['info_name'])
        fout.write('  <description>%s</description>\n' % conf['info_description'])
        fout.write('  <version>%s</version>\n' % conf['info_version'])
    
        for act in conf['actions']:
            fout.write('  <action>\n')
            fout.write('    <plugin>%s</plugin>\n' % act['plugin'])
            fout.write('    <function>%s</function>\n' % act['function'])
            fout.write('    <default_level>%s</default_level>\n'% act['function'])
        
            fout.write('    <params>\n')
            for pkey,pval in act['params'].items():
                fout.write('        <%s>%s</%s>\n' %(pkey,pval,pkey))
            fout.write('    </params>\n')
        
            fout.write('    <limits>\n')
            for pkey,pval in act['limits'].items():
                fout.write('        <result description="%s">\n' % pkey)
                for ppkey,ppval in pval.items():
                    fout.write('            <%s>%s</%s>\n' %(ppkey,ppval,ppkey))
                fout.write('        </result>\n')
            fout.write('    </limits>\n')
        
            fout.write('  </action>\n')
    
        fout.write('</Input2DCM_config>\n')
    