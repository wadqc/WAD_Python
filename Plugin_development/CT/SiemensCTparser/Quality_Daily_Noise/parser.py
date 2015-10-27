# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
# Description:
# This plugin parses quality_daily_noise file

__version__ = '20082015'
__author__ = 'TdW'



from dicom import tag
import xml.etree.ElementTree as ET
import lxml.etree as etree

def print_xml(xmlroot):
   for child in xmlroot:
        print '=='*20
        print child.tag, child.attrib, child.text
        
        for subchild in child:
            print '\t', subchild.tag, subchild.attrib, subchild.text

            for value in subchild:
                print '\t\t', value.tag, value.attrib, value.text
                for subvalue in value:
                    print '\t\t\t', subvalue.tag, subvalue.attrib, subvalue.text

                    for subsubvalue in subvalue:
                        print '\t\t\t\t', subsubvalue.tag, subsubvalue.attrib, subsubvalue.text

def parseqcreport(data,results,**kwargs):
    params = kwargs.get('params', None)
    p = {}
    for param in params:
        p[param.tag] = (param.text,param.attrib)

    print p
    
    status_ok=True

    relevantfile = data.getAllInstances()[0]
    xmltext = relevantfile[tag.Tag(p.get('use_private_tag')[0].split(','))]

    root = etree.fromstring(xmltext.value)
    print_xml(root)
    
    datagroup = root.find('DataGroup')
    scandate = datagroup.get('DateTime')
    
    datapackages = root.find(".//DataPackage[@Package_Name='ID_FCT_SUMMARY']")

    summary = datapackages.find(".//DataEntry[@Variable_Name='ID_ErrorCodeSummary']")
    outval = summary.find("StringVariable[@String_Value]").attrib['String_Value']
    results.addFloat('ID_ErrorCodeSummary',outval,level=1)

    measurementsroot = root.find("DataGroup[@Group_Name='ID_GROUP']")
    measurements = measurementsroot.findall(".//DataPackage[@Package_Name='ID_MODE']")

    
    print '--'*20
    for measurement in measurements:
       print 'New measurement'

       currvolt = measurement.find(".//DataEntry[@Variable_Name='ID_VOLT']")
       currcurrent = measurement.find(".//DataEntry[@Variable_Name='ID_CURRENT']")

       print currcurrent.find('LongVariable').attrib['Long_Value']
       print currvolt.find('LongVariable').attrib['Long_Value']

       try:
          print measurement.find(".//DataArray[@Array_Name='ID_RESULT']").find("DataSeries[@Series_Name='ID_SIGMA']").find('.//DoubleVariable').attrib['Double_Value']

       except:
          pass

       print measurement.find(".//DataArray[@Array_Name='ID_RESULT']").find("DataSeries[@Series_Name='ID_VOLT']").find('.//DoubleVariable').attrib['Double_Value']


       
