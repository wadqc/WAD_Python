#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# PyWAD is open-source software and consists of a set of modules written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes modules for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Changelog:
#   20160901: initial version based on parts of pyWAD code
#
# This program is used to test a plugin (through a config.xml) with a given dataset.
# TODO: make it run for separate series, and separate studies, and separate instances. write to config/series/instance/result.xml, obj.jpg
#
# ./run_module.py -c ~/MyConfigsXML/us_philips_epiq_instance.xml -d ../wad_new_modules/US/AirReverberations/TestSet/StudyEpiqCurve/ -r results_epiq.xml
from __future__ import print_function

__version__ = '20160901'
__author__ = 'aschilham'

import os
import sys
import argparse
from pyWAD.xml_result_generator import getXMLFromResult
from pyWAD.xml_parser import parseConfig
from pyWAD.plugindata import PluginData
from pyWAD.pluginresults import PluginResults
    
def _run(action_list, series_list, out_fn):
    """
    for each action, run the correct python module and function, and add results to file
    """
    result_xml_list = []
    for action in action_list:
        # Import plugin function
        module = __import__(action['plugin'], fromlist=[action['function']])
        function = getattr(module, action['function'])
        
        # Get plugin version
        version = "%s (%s)"%(module.__version__, action['plugin'])

        # Initiate objects
        plugin_data = PluginData(series_list, action['filters'])
        plugin_results = PluginResults(action['default_level'])

        # Evaluate function
        print("----- Evaluating {} in {} -----".format(
            action['function'], action['plugin']))
        function(plugin_data, plugin_results, params=action['params'])

        # Report plugin version
        if plugin_results:
            plugin_results.addChar("Version", version, level=2)

        # Convert results to XML
        for result in plugin_results:
            result_xml_list.append(getXMLFromResult(result, action['limits']))
            print(result)

    if result_xml_list:
        outXML = '\n'.join(['<WAD>'] + result_xml_list + ['</WAD>'])
        with open(out_fn, 'w') as f:
            f.write(outXML)
        
            
def get_series_list(studyfolder):
    """
    parsed_input['serieslist'] is a list of series, where each series is a list of instance filenames
    """
    series_list = []
    for series in os.listdir(studyfolder):
        sfolder = os.path.join(studyfolder, series)
        instance_list = os.listdir(sfolder)
        series_list.append([os.path.abspath(os.path.join(sfolder, fn)) for fn in instance_list])
    return series_list
    
def run(configfile, studyfolder, resultfile):
    """
    a callable function for direct access
    """
    action_list = parseConfig(configfile)
    series_list = get_series_list(studyfolder)
    
    _run(action_list, series_list, resultfile)
    
    
def pyWADinput():
    """
    convert command line flags to needed files and folder.
    extract action_list from config_file and generate series_list from studyfolder.
    return arguments for run.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', help='Dicom study folder')
    parser.add_argument('-r', help='Results file')
    parser.add_argument('-c', help='Config file')
    args = parser.parse_args()

    if args.c is None or args.d is None or args.r is None:
        parser.print_help()
        sys.exit()
        
    return args.c, args.d, args.r

    
if __name__ == "__main__":
    configfile, studyfolder, out_fn = pyWADinput()

    run(configfile, studyfolder, out_fn)
