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
# This is a very simple plugin to output all the dicom header tags as results to the WAD server.

__version__='01062015'
__author__='DD'



import dicom
import re

def print_header(data, results, **kwargs):
    '''
    Function receives a data class object and from the header of the first instance it creates a dictionary

    Workflow:
      2. For each tag in the header take the first 500 chars and remove unwanted characters
      3. See if there is a keyword description for the tag, if so store it in variable q
      4. Write out results fixed to level 3.

    '''

    header = data.getAllInstances()[0]

    for key in header.keys():
        stringedkey = str(key).replace('(','').replace(')','').replace(' ','')[:127]
        quantity = dicom.datadict.keyword_for_tag(key)[:127]
        print 'plugin;',key, quantity, re.sub(r'[^\w]','',header[key].repval)[:127]
        s = re.sub(r'[^\w]','',header[key].repval)[:127]
        results.addChar(stringedkey,s,level=2,quantity=quantity)
