# -*- coding: utf-8 -*-
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
from __future__ import print_function

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!
"""
OCR using tessaract or cuniform through pyOCR.
pytesseract is also fine to use, and does not need zooming (but results are less), while pyocr needs zoom at least 6.0.
Take care to put the bounding box around only txt, and exclude other objects!

On Ubuntu 16.04: apt install python-pyocr tesseract-ocr tesseract-ocr-eng
"""
__version__ = '20160831'
__author__ = 'aschilham'

from PIL import Image
import numpy as np
from scipy import ndimage as scind
import pyocr
import re
import scipy.misc
# sanity check: we need at least scipy 0.10.1 to avoid problems mixing PIL and Pillow
scipy_version = [int(v) for v in scipy.__version__ .split('.')]
if scipy_version[1]<10 or (scipy_version[1] == 10 and scipy_version[1]<1):
    raise RuntimeError("scipy version too old. Upgrade scipy to at least 0.10.1")


def getOCRTool():
    # check if OCR tools tesseract or cuneiform are available
    tools = pyocr.get_available_tools()
    if len(tools) == 0:
        raise RuntimeError("ERROR No OCR tool found")
    tool = tools[0]
    print("[ocr_lib] Using %s for OCR" % (tool.get_name()))
    return tool


def txt2type(txt, type, prefix='',suffix=''):
    """
    If prefix is defined, the length of the string is used to skip the first num characters.
    If suffix is defined, the length of the string is used to skip the last num characters.
    """
    if not type in ['string', 'bool', 'float']:
        raise ValueError('Unknown type %s'%type)

    txt = txt.strip() # removes '\n'
    txt = txt[len(prefix):]
    if len(suffix)>0:
        txt = txt[:-len(suffix)]
    txt = txt.strip()

    if len(txt) == 0:
        raise ValueError('[ocr_lib] ERROR! empty text line!')

    if type.lower() == 'string':
        return txt

    if type.lower() == 'bool':
        return (txt.lower() in ['1', 'true', 'y', 'yes'])
        
    # strip non-numeric characters from floating number
    if type.lower() == 'float':
        # first strip % and spaces and turn comma into period (without warning)
        txt = txt.replace('%','').replace(' ','').replace(',', '.')
        # next the other characters
        newtxt = re.sub(r'[^\d.]+', '', txt)
        if newtxt != txt:
            print(u"[ocr_lib] Warning: replaced value {} by {}".format(txt, newtxt).encode('utf8'))
            txt = newtxt
        return float(txt)

    
def OCR(pixeldata, xywh, zpos=0, ocr_zoom=10, ocr_threshold=0):
    """
    Use pyOCR which for OCR
    ul = upperleft pixel location [x,y]
    ocr_zoom = factor to enlarge image (15)
    ocr_threshold = remove values below this threshold (after inversion)
    """
    tool = getOCRTool()

    # slice-out the relevant part of the image
    x,y,width,height = xywh

    if len(np.shape(pixeldata)) == 3:
        part = pixeldata[zpos][x:x+width,y:y+height]
    elif len(np.shape(pixeldata)) == 2:
        part = pixeldata[x:x+width,y:y+height]
    else:
        raise ValueError('[ocr_lib] Unknown dataformat of %d dimensions'%len(np.shape(pixeldata)))

    # heuristic contrast enhancement: want white is txt, background = 0
    # invert if needed
    edgeval = (np.mean(part[0,:])+np.mean(part[-1,:]))/2
    if edgeval > 128:
        part = edgeval-part
        part[part<0] = 0

    # remove noise/gradient
    part[part<ocr_threshold] = 0

    # enhance contrast
    minval = np.min(part)
    maxval = np.max(part)
    if (maxval-minval)<128:
        part = (part-minval)*(255/(maxval-minval))

    part = np.transpose(part) # input was pyqtgraph-like

    if ocr_zoom is None:
        # enlarge to prevent OCR mismatches; below 20px font height accuracy drops off
        minheight = 200 # this value to prevent pixGenHalftoneMask errors
        minwidth = 600 # this value to prevent pixGenHalftoneMask errors
        if height<minheight or width<minwidth:
            minzoom = int(max([minheight/height,minwidth/width]))+1
            ocr_zoom = minzoom

    if not ocr_zoom is None:
        part = np.round(scind.interpolation.zoom(part, zoom=(ocr_zoom,ocr_zoom),order=1))

    # extract numbers/text from bounding box
    ##import pytesseract
    ##txt = pytesseract.image_to_string(Image.fromarray(part))
    
    txt = tool.image_to_string(Image.fromarray(part))
    return txt, part
