'''
This plugin parses getinge washer QC reports
'''

import dicom
from dicom import tag

__version__='23062015'


def QC_getinge_main(data,results, **kwargs):



     relevantfile = data.getAllInstances()[0]
     print 'Reading dicomfile', relevantfile
     
             
     dicomobject = relevantfile
     tmpseriesdescription =  dicomobject[tag.Tag("0008","103E")].value

     waarde = dicomobject[tag.Tag("0071","0001")].value
     results.addChar('Result Thermal selfdesinfection',waarde ,level='1')

     waarde = dicomobject[tag.Tag("0071","0003")].value
     results.addChar('Error code',waarde ,level='1')

     waarde = dicomobject[tag.Tag("0071","0005")].value
     results.addChar('Error text',waarde)
     







