# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#

__version__ = '20092015'
__author__ = 'DD'


from dicom import tag
import xml.etree.ElementTree as ET
import lxml.etree as etree


def parseqcreport(data,results,**kwargs):

    # 1 - Read parameters from config file:
    params = kwargs.get('params', None)
    p = {}
    for param in params:
        p[param.tag] = (param.text,param.attrib)
    print p
    
    # 2 - Load XML data
    relevantfile = data.getAllInstances()[0]
    xmltext = relevantfile[tag.Tag(p.get('use_private_tag')[0].split(','))]
    root = etree.fromstring(xmltext.value)
    
    # 5 - WRITE CODE TO FIND THE TAGS OF INTEREST HERE:
    # For example mydata = root.find("MyXMLTag[@MyXMLTag_Name='SOME NAME']")

    # 6 - Write results to WAD software
    # For example results.addFloat('MyParameter',value,level=1)
    
