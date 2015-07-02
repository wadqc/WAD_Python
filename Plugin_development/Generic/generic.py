# -*- coding: utf-8 -*-
import sys
import dicom
from dicom.dataset import Dataset, FileDataset
from ast import literal_eval

from pyWAD.plugindata import PluginData
from pyWAD.pluginresults import PluginResults



def main(data, results, **kwargs):
    """ Function extracts results from dicom file, using the private dicom tag
        configured in the config-xml [private_group,private_element] and writes
        the results to the results-object for further processing.
        
        For an example, check:
            Tools/Generic/create_dicom.py      (creates a dicom file)
            Testing/Generic/generic_input.xml  (configures where and what to extract from the dicom file)
    """

    parameters = kwargs.get('params', None)

    p = {}
    for param in parameters:
        p[param.tag] = param.text

    tag_type = p.get('tag_type')
    private_group = p.get('private_group')
    private_element = p.get('private_element')

    # currently only one instance supported
    try:
        dcm=data.getAllInstances()[0]
    except:
        print 'Error reading input data!'
        sys.exit(1)

    result_data = dcm[(private_group,private_element)].value

    # probably not needed (but just in case): strip newlines, returns and tabs
    result_data = result_data.replace('\n', '').replace('\r', '').replace('\t','')

    result_dict = literal_eval(result_data)

    data_description = result_dict.keys()[0]
    print 'Data description = %s' %data_description
    data_list = result_dict[data_description]


    def translate_type_to_function(type):
       return {
         'float'  : 'addFloat',
         'char'   : 'addChar',
         'bool'   : 'addBool',
         'object' : 'addObject',
       }.get(type)


    # extracting data from post and adding them to the results object
    for m in data_list:
        #print measurement['description']
        if m.get('type')=='object':
            # convert base64 string to binary file
            outdata = m.get('base64_blob').decode("base64")
            outfile = open(m.get('filename'), "wb")
            outfile.write(outdata)
            outfile.close()
            print 'adding result for', m.get('description')
            results.addObject(m.get('description'),m.get('filename'),m.get('level'),m.get('quantity'),m.get('units'))
        elif len(m.get('value'))>0:
            function_name=translate_type_to_function(m.get('type'))
            if function_name is not None:
                print 'adding result for', m.get('description')
                add_result=getattr(results,function_name)
                add_result(m.get('description'),m.get('value'),m.get('level'),m.get('quantity'),m.get('units'))
            else:
                print "Error: Unrecognized datatype '%s' (allowed: float, bool, object, char)!"%type
            
