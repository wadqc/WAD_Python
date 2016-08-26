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
# 
#

__version__='01062015'
__author__ = 'tdw'


# requirements: 
#   tesseract (ubuntu: apt-get install tesseract-ocr tesseract-ocr-eng
#              macosx: port install tesseract tesseract-eng)
#   pyocr (pip install pyOCR)

import sys
from PIL import Image
from scipy import ndimage
import numpy as np
import dicom
import pyocr
from ast import literal_eval
import pylab as pl
import re
import scipy.misc
# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
    raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")

from pyWAD.plugindata import PluginData
from pyWAD.pluginresults import PluginResults


def extract(data, results, **kwargs):
    """ Function extracts data (text, numbers, images) from dicom secondary capture
        images by using OCR.	    
    """

    parameters = kwargs.get('params', None)
    
    # read in measurement parameters from xml, containing
    # bounding box information required for OCR (upper-left coordinate,
    # width, height as measured with ImageJ (XY swapped!)),
    # level, quantity, units, etc, for extraction and adding to results object.

    measurementlist = []

    xstr = lambda s: '' if s is None else str(s)

    for param in parameters:
       if 'measurement' in param.tag:
          # print param.tag
          measurement = {}
          for subparam in param:
              # print '    ',subparam.tag,'=',subparam.text
              measurement[subparam.tag]=subparam.text

          type=measurement.get('type')
          description=measurement.get('description')
          UL=literal_eval(measurement.get('upper_left_coordinates'))
          width=literal_eval(measurement.get('width'))
          height=literal_eval(measurement.get('height'))
          level=literal_eval(measurement.get('level'))
          quantity=xstr(measurement.get('quantity'))
          units=xstr(measurement.get('units'))
          filename=measurement.get('filename')
          measurementlist.append((param,type,description,UL,width,height,level,quantity,units,filename))


    # check if OCR tools tesseract or cuneiform are available
    tools = pyocr.get_available_tools()
    if len(tools) == 0:
        print("No OCR tool found")
        sys.exit(1)

    tool = tools[0]
    print("Using %s for OCR" % (tool.get_name()))

    # load secondary capture (currently only one instance supported)
    dcm=data.getAllInstances()[0]
    img = dcm.pixel_array

    # restructure dicom RGB data to python RGB data
    img = img.reshape(img.shape[1],img.shape[2],img.shape[0])


    def translate_type_to_function(type):
       return {
         'float'  : 'addFloat',
         'char'   : 'addChar',
         'bool'   : 'addBool',
         'object' : 'addObject',
       }.get(type.lower())


    for measurement in measurementlist:
       subparam=measurement[0]
       type=measurement[1]
       description=measurement[2]
       UL=measurement[3]
       width=measurement[4]
       height=measurement[5]
       level=measurement[6]
       quantity=measurement[7]
       units=measurement[8]
       filename=measurement[9]
  
       # slice-out the relevant part of the image and enlarge to prevent OCR mismatches
       bounding_box=img[UL[1]:UL[1]+height,UL[0]:UL[0]+width,:]

       if type == 'object':
          #pl.imsave(filename,bounding_box)
	  im = scipy.misc.toimage(bounding_box) 
	  im.save(filename)


          results.addObject(description,filename,level,quantity,units)
       else:
          bounding_box=np.round(ndimage.interpolation.zoom(bounding_box, zoom=(6,6,1),order=1))

          # extract numbers/text from bounding box
          value=tool.image_to_string(Image.fromarray(bounding_box))

          # translate OCR-boolean to correct value
          if type.lower() == 'bool':
             if value.lower() in ['yes','ja','1','true','waar','no','nee','0','false','onwaar']:
                 value = value.lower() in ['yes','ja','1','true','waar']

          # strip non-numeric characters from floating number
          if type.lower() == 'float':
             # first strip % and spaces (without warning)
             value=value.replace('%','').replace(' ','')
             # next the other characters
             newvalue=re.sub(r'[^\d.]+', '', value)
             if newvalue!=value:
                 print u"Warning: replaced value {} by {}".format(value,newvalue).encode('utf8')
             value=newvalue

          # if value present and extraction succesful, add to results
          if len(value)>0:
             function_name=translate_type_to_function(type)
             if function_name is not None:
                add_result=getattr(results,function_name)
                try:
                   add_result(description,value,level,quantity,units)
                except Exception as e:
                   print "Exception occurred for parameter '{}': {}".format(description,e)
             else:
                print "Error: Unrecognized datatype '%s' (allowed: float, bool, object, char)!" %type
          else:
             print "Warning: no value for parameter '%s'" %description
