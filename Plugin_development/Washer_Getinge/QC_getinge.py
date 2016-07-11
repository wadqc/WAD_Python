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
     







