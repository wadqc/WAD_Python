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
# This is a simple plugin to print the location of the dcm4chee dataset

__version__='01062015'
__author__='DD'

import dicom

def fileLocation(data,results,**kwargs):
    '''
    Function receives a data class object and returns for each series a link to the first instance
    '''

    relevantFileList = data.series_filelist

    
    for series in relevantFileList:
        tmpSeriesDescription = dicom.read_file(series[0]).SeriesDescription
        results.addChar('Series',series[0],level='2',quantity=tmpSeriesDescription)
        

    
