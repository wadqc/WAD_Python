# -*- coding: utf-8 -*- (specification required to test special characters)
"""
Example plugin which demonstrates the use of the data, results and (optional)
params objects.

Author: Rob van Rooij
"""

__version__ = "1.0" # Mandatory, is automatically added to the results object

import os
import sys


def testFunction(data, results, **kwargs):
    """Example function.
    """
    print('Number of instances:', len(data.getAllInstances()))
    instance = data.getAllInstances()[0]
    study_date = instance.StudyDate    
    
    #special_characters = 'Caractères spéciaux'  # Python 3
    special_characters = u'Caractères spéciaux'  # Python 2
    
    ### Just testing... See how pydicom manages 
    print('Type     ', type(instance.ExposureTime))
    print('Is int   ', instance.ExposureTime == 650)
    print('Is str   ', instance.ExposureTime == '650')
    
    print('Type     ', type(instance.Rows))
    print('Is int   ', instance.Rows == 4096)
    print('Is str   ', instance.Rows == '4096')
    
    print('Type     ', type(instance.HalfValueLayer))
    print('Is float ', instance.HalfValueLayer == 0.351)
    print('Is str   ', instance.HalfValueLayer == '0.351')

    results.addChar('Plugin location', os.path.realpath(__file__))
    results.addChar('Study date', study_date, level=2)
    results.addChar('Special Characters', special_characters)


def anotherFunction(data, results, params):
    """Example function.
    """
    print('Number of instances:', len(data.getAllInstances()))
    print('Filtered instances:', len(data.getInstanceByTags({
        "SeriesDescription": "Flat Field",
        (0x0008, 0x1030): "Weekly QC",
        "0x0008, 0x1090": "Lorad Selenia",
        }))
    )
    print([len(x) for x in data.getAllSeries(stop_before_pixels=True)])
    print([len(x) for x in data.getSeriesByDescription("Flat Field", force=True)])

    if params:
        print('\nParameters:')
        print(ET_to_XML(params))

        print('\nOr converted to lists:')
        print(ET_to_lists(params))
        print()

    nr_of_instances = len(data.getAllInstances())
    plaatje = 'image.png'

    results.addFloat('Number of instances', nr_of_instances, quantity='Amount', units='instances', level=2)
    results.addObject('Image', plaatje, level=1)
    results.addBool('Some bool', True, level=2)


def ET_to_lists(element):
    """Return a nested list from an ElementTree Element

    Arguments:
        element: xml.etree.ElementTree.Element object
    """
    if not len(element):
        #If the ElementTree object has no children below, return the text contents of the element
        return element.text

    lst = []
    for child in element:
        #Add a sub-list for every child element, add the contents of the child recursively
        lst.append([child.tag, child.attrib, ET_to_lists(child)])

    return lst


def ET_to_XML(element, indent=4, level=0):
    """Return the XML string corresponding to an ElementTree Element

    Arguments:
        element: xml.etree.ElementTree.Element object
        indent (int): how many spaces to indent per level
        level (int): level of indentation
    """

    #Join any attributes into a string of form: ' attr1="val1" attr2="val2'
    attrib_str = "".join([' {}="{}"'.format(k, v) for k, v in element.attrib.items()])

    if not len(element):
        #If this ET element has no children, return string of form <tag attr1="val1" attr2="val2">text</tag>
        return "{0}<{1}{2}>{3}</{1}>".format(
            level * indent * " ",
            element.tag,
            attrib_str,
            element.text
            )

    XML_list = []
    XML_list.append("{}<{}{}>".format(
        level * indent * " ",
        element.tag,
        attrib_str
    ))  # Add parent tag
    for child in element:
        XML_list.append(ET_to_XML(child, indent, level + 1))  # Add children recursively, indented by one
    XML_list.append("{}</{}>".format(level * indent * " ", element.tag))  # Close parent tag

    #Join all the lines with line-breaks before returning
    return "\n".join(XML_list)


if __name__ == "__main__":
    from pyWAD import PluginData, PluginResults

    series_lst = [['./dicom.dcm']]

    plugin_data = PluginData(series_lst)
    plugin_results = PluginResults(default_level=1)

    testFunction(plugin_data, plugin_results)
    anotherFunction(plugin_data, plugin_results, params={}, parameters=None)

    for result in plugin_results:
        print(result)
