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
from __future__ import print_function

"""
Normi13 analysis, Low Contrast functions:
  o Mean, SNR, STDEV
  o Dynamic Range (max step/min step)
"""
import numpy as np

class LoCoStruct :
    def __init__ (self):
        self.low_rois = []     # circles of low contrast [x0,y0,r] in pix
        self.low_rois_bku = [] # circles of background above low contrast [x0,y0,r] in pix
        self.low_rois_bkd = [] # circles of background above below low contrast [x0,y0,r] in pix

        self.mean_sg = []
        self.sdev_sg = []
        self.mean_bk = [] # corresponding to max cnr
        self.sdev_bk = [] # corresponding to max cnr
        self.low_cnr = [] # max cnr for up and down

def LowContrast(cs):
    """
    for each disc, calc av in box in center and in box 1 cm above it (without Cu grid)
    Calculate CNR: (av(Cu)-av(bk))/sqrt((sd(Cu)^2+sd(bk)^2)/2)
    """
    error = True

    # discs are 1 cm in diameter
    diampx = cs.phantommm2pix(6.)
    radpx  = int(diampx/2.+.5)
    yrefmm = 10. # offset to background

    low_rois = [] # [x0,y0,rad]
    low_rois_bku = [] # [x0,y0,rad] # rois above for bk
    low_rois_bkd = [] # [x0,y0,rad] # rois below for bk
    mean_sg = [] # average in object
    sdev_sg = [] # average in object

    dxmm = 17. # dist in mm between loco objects
    dypx = cs.phantommm2pix(11.) # dist in px between loco objects
    for i in range(6):
        px,py = cs.geom.phantomposmm2pix(-42.5+i*dxmm,-73.)
        low_rois.append([px,py,1.*radpx])
        low_rois_bku.append([px,py-dypx,1.*radpx])
        low_rois_bkd.append([px,py+dypx,1.*radpx])

    imwid,imhei = np.shape(cs.pixeldataIn)
    mean_bku = []
    mean_bkd = []
    sdev_bku = []
    sdev_bkd = []
    calcsets = [
        (low_rois, mean_sg, sdev_sg),
        (low_rois_bku, mean_bku, sdev_bku),
        (low_rois_bkd, mean_bkd, sdev_bkd)
    ]
    for rois, means, sdevs in calcsets:
        for r in rois:
            [px,py,radpx] = [int(rr) for rr in r]
            xmin = max(0,px-radpx)
            ymin = max(0,py-radpx)
            xmax = min(px+radpx, imwid-1)
            ymax = min(py+radpx, imhei-1)
            if ymax<ymin: # field too small, no low contrast elements
                error = True
                return
            # select data inside circle only 
            smallimage = cs.pixeldataIn[xmin:xmax+1, ymin:ymax+1]

            x,y = np.indices((xmax-xmin+1, ymax-ymin+1))
            mask = ((x+xmin-px)**2 + (y+ymin-py)**2 ) < radpx**2
            means.append(np.mean(smallimage[mask]))
            sdevs.append(np.std(smallimage[mask]))

    low_cnr = []
    mean_bk = []
    sdev_bk = []
    # calculate for each object the max contrast to noise wrt upper or lower background. Keep the values for the max cnr
    for i in range(len(low_rois_bkd)):
        cnr_up = (mean_bku[i]-mean_sg[i])/np.sqrt(0.5*(sdev_sg[i]**2+sdev_bku[i]**2))
        cnr_do = (mean_bkd[i]-mean_sg[i])/np.sqrt(0.5*(sdev_sg[i]**2+sdev_bkd[i]**2))
        if cnr_up>cnr_do:
            low_cnr.append(cnr_up)
            mean_bk.append(mean_bku[i])
            sdev_bk.append(sdev_bku[i])
        else:
            low_cnr.append(cnr_do)
            mean_bk.append(mean_bkd[i])
            sdev_bk.append(sdev_bkd[i])

        if cs.verbose:
            print("mean fg/bk=",mean_sg[i],"/",mean_bk[i])
            print("sdev fg/bk=",sdev_sg[i],"/",sdev_bk[i])
            print('cnr = ',low_cnr[i])

    # copy to loco struct
    cs.loco.low_rois = low_rois 
    cs.loco.low_rois_bku = low_rois_bku
    cs.loco.low_rois_bkd = low_rois_bkd 

    cs.loco.low_cnr = low_cnr
    cs.loco.mean_sg = mean_sg
    cs.loco.mean_bk = mean_bk
    cs.loco.sdev_sg = sdev_sg
    cs.loco.sdev_bk = sdev_bk

    error = False
    return error
