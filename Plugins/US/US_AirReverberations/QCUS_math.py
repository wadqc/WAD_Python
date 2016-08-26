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
import scipy.ndimage as scind
import numpy as np
from scipy import optimize

class CircleFit:
    # Coordinates of the 2D points
    x = []
    y = []
    residu = -1
    
    def __init__(self,data_xy):
        self.x = [xy[0] for xy in data_xy]
        self.y = [xy[1] for xy in data_xy]
        self.residu = -1.
        
    def calc_R(self, xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((self.x-xc)**2 + (self.y-yc)**2)

    def f_2(self,c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = self.calc_R(*c)
        return Ri - Ri.mean()
   
    def fit(self):
        method_2 = "leastsq"
        
        center_estimate = np.mean(self.x),np.mean(self.y)
        center_2, ier = optimize.leastsq(self.f_2, center_estimate)
        
        xc_2, yc_2 = center_2
        Ri_2       = self.calc_R(*center_2)
        R_2        = Ri_2.mean()
        residu_2   = sum((Ri_2 - R_2)**2)

        self.residu = residu_2
        return center_2,R_2
    

def localnormalization(pSrc, sigma,bksigma = None):
    """
    Local normalization: [ i(x,y)-mean(x,y) ]/stdev(x,y)
    Can be approximated as [ I- Gauss{I,sigma}] / sqrt[ Gauss{I-Gauss{I,sigma},sigma}^2]
    """
    blurIm = scind.gaussian_filter1d(pSrc,sigma,order=0)
    devIm = pSrc-blurIm

    sdIm  = np.sqrt(scind.gaussian_filter1d(devIm**2,sigma,order=0))
    sdIm[sdIm<1.e-6]=1. # prevent div by zero

    locnormIm = devIm/sdIm
    return locnormIm

def movingaverage(data1d, window_size):
    if window_size <2:
        return data1d
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(data1d, window, 'same')    

