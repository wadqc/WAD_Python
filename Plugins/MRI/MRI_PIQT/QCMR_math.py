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

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:29:57 2013

AffineTransformation Originally written by Jarno Elonen <elonen@iki.fi>, November 2007, released into the Public Domain
http://elonen.iki.fi/code/misc-notes/affine-fit/

Slightly rewritten by aschilha
"""
import numpy as np
from scipy.optimize import curve_fit
import scipy.ndimage as scind
import scipy.stats


class AffineTransformation:
    """Result object that represents the transformation
       from affine fitter."""
    M = []
    dim = 0
    def __init__ (self, matrix,pdim):
        self.M = matrix
        self.dim = pdim

    def To_Str(self):
        res = ""
        for j in range(self.dim):
            str = "x%d' = " % j
            for i in range(self.dim):
                str +="x%d * %f + " % (i, self.M[i][j+self.dim+1])
            str += "%f" % self.M[self.dim][j+self.dim+1]
            res += str + "\n"
        print self.SimpleMatrix()
        return res

    def SimpleMatrix(self):
        mat = []
        for j in range(self.dim):
            row = []
            for i in range(self.dim):
                row.append(self.M[i][j+self.dim+1])
            row.append(self.M[self.dim][j+self.dim+1])
            mat.append(row)
        return mat
        
    def Transform(self, pt):
        res = [0.0 for a in range(self.dim)]
        for j in range(self.dim):
            for i in range(self.dim):
                res[j] += pt[i] * self.M[i][j+self.dim+1]
            res[j] += self.M[self.dim][j+self.dim+1]
        return res

    def getShift(self):
        shift = [ self.M[self.dim][j+self.dim+1] for j in xrange(self.dim)]
        return shift
        
    def getScaling(self):
        mat = self.SimpleMatrix()
        scaling = [np.sqrt(mat[0][0]**2+mat[0][1]**2),np.sqrt(mat[1][0]**2+mat[1][1]**2)]
        return scaling

    def getRotation(self):
        mat = self.SimpleMatrix()
        rot = np.arctan2(mat[0][0],mat[0][1])
        return rot

    def RigidTransform(self, pt):
        mat = self.SimpleMatrix()
        scaling = self.getScaling()
        row = mat[0]
        nx = row[0]/scaling[0]*pt[0]+row[1]/scaling[0]*pt[1]+row[2]
        row = mat[1]
        ny = row[0]/scaling[1]*pt[0]+row[1]/scaling[1]*pt[1]+row[2]
        
        return [nx,ny]

def Affine_Fit( from_pts, to_pts ):
    """Fit an affine transformation to given point sets.
      More precisely: solve (least squares fit) matrix 'A'and 't' from
      'p ~= A*q+t', given vectors 'p' and 'q'.
      Works with arbitrary dimensional vectors (2d, 3d, 4d...).

      Written by Jarno Elonen <elonen@iki.fi> in 2007.
      Placed in Public Domain.

      Based on paper "Fitting affine and orthogonal transformations
      between two sets of points, by Helmuth Spath (2003)."""

    """
    Example:
    from_pt = ((1,1),(1,2),(2,2),(2,1)) # a 1x1 rectangle
    to_pt = ((4,4),(6,6),(8,4),(6,2))   # scaled x 2, rotated 45 degrees and translated
    
    trn = Affine_Fit(from_pt, to_pt)
    
    print "Transformation is:"
    print trn.To_Str()
    
    err = 0.0
    for i in range(len(from_pt)):
        fp = from_pt[i]
        tp = to_pt[i]
        t = trn.Transform(fp)
        print ("%s => %s ~= %s" % (fp, tuple(t), tp))
        err += ((tp[0] - t[0])**2 + (tp[1] - t[1])**2)**0.5
    
    print "Fitting error = %f" % err
    """
    q = from_pts
    p = to_pts
    if len(q) != len(p) or len(q)<1:
        print "from_pts and to_pts must be of same size."
        return False

    dim = len(q[0]) # num of dimensions
    if len(q) < dim:
        print "Too few points => under-determined system."
        return False

    # Make an empty (dim) x (dim+1) matrix and fill it
    c = [[0.0 for a in range(dim)] for i in range(dim+1)]
    for j in range(dim):
        for k in range(dim+1):
            for i in range(len(q)):
                qt = list(q[i]) + [1]
                c[k][j] += qt[k] * p[i][j]

    # Make an empty (dim+1) x (dim+1) matrix and fill it
    Q = [[0.0 for a in range(dim)] + [0] for i in range(dim+1)]
    for qi in q:
        qt = list(qi) + [1]
        for i in range(dim+1):
            for j in range(dim+1):
                Q[i][j] += qt[i] * qt[j]

    # Ultra simple linear system solver. Replace this if you need speed.
    def gauss_jordan(m, eps = 1.0/(10**10)):
      """Puts given matrix (2D array) into the Reduced Row Echelon Form.
         Returns True if successful, False if 'm' is singular.
         NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
         Written by Jarno Elonen in April 2005, released into Public Domain"""
      (h, w) = (len(m), len(m[0]))
      for y in range(0,h):
        maxrow = y
        for y2 in range(y+1, h):    # Find max pivot
          if abs(m[y2][y]) > abs(m[maxrow][y]):
            maxrow = y2
        (m[y], m[maxrow]) = (m[maxrow], m[y])
        if abs(m[y][y]) <= eps:     # Singular?
          return False
        for y2 in range(y+1, h):    # Eliminate column y
          c = m[y2][y] / m[y][y]
          for x in range(y, w):
            m[y2][x] -= m[y][x] * c
      for y in range(h-1, 0-1, -1): # Backsubstitute
        c  = m[y][y]
        for y2 in range(0,y):
          for x in range(w-1, y-1, -1):
            m[y2][x] -=  m[y][x] * m[y2][y] / c
        m[y][y] /= c
        for x in range(h, w):       # Normalize row y
          m[y][x] /= c
      return True

    # Augement Q with c and solve Q * a' = c by Gauss-Jordan
    M = [ Q[i] + c[i] for i in range(dim+1)]
    if not gauss_jordan(M):
        print "Error: singular matrix. Points are probably coplanar."
        return False

    # Make a result object
    transformation = AffineTransformation(M,dim)
    return transformation
    

def FindCenters2D(pts,datain,distpx,discpx,minimod=False):
    error = True
    searchrad = int(distpx+.5)
    searchrad = max(1,searchrad)

    # scale by 4x4 by copies for 1/4 pix accuracy
    multip = 4
    data = np.kron(datain,np.ones((multip,multip)))

    # smoothing to get rid of noise and give max respons over avg disc
    sigma = multip*discpx/2.
    data = scind.gaussian_filter(data.astype(float), sigma,mode='constant')
    searchrad *= multip
    """"
    # scale by 2x2 by copies
    a = np.array([[1, 1],
                  [0, 1]])
    n = 2
    np.kron(a, np.ones((n,n)))
    """
    
    widthpx = np.shape(data)[0] ## width/height in pixels
    heightpx = np.shape(data)[1]

    for y in range(len(pts)):
        for x in range(len(pts[y])):
            rp = pts[y][x]
            if(len(rp)==0):
                continue
            x0 = multip*rp[0]
            minx = max(0,x0-searchrad)
            maxx = min(widthpx-2,x0+searchrad)
            y0 = multip*rp[1]
            miny = max(0,y0-searchrad)
            maxy = min(heightpx-2,y0+searchrad)
            cropped = data[minx:maxx+1,miny:maxy+1]
            if(minimod == True):
                (x1,y1) = np.unravel_index(cropped.argmin(),cropped.shape)
            else:
                (x1,y1) = np.unravel_index(cropped.argmax(),cropped.shape)
            x1 += minx
            y1 += miny
            rp[0] = x1/multip
            rp[1] = y1/multip  

    error = False
    return error,pts


# Define model function to be used to fit to the data 
def gauss(x, *p):
    A, mu, sigma,c = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))+c

def GaussianFit(data):
    error = True
    
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [np.max(data)-np.min(data), len(data)/2., 1.,data[0]]

    pos = range(len(data))
    coeff, var_matrix = curve_fit(gauss, pos, data, p0=p0)
    error = False
    return error, coeff

def LinearFit(ydata):
    error = True
    pos = range(len(ydata))
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(pos,ydata)

    return error,(slope, intercept, r_value, p_value, std_err)

def AreaUnderCurve(freqs,merits,maxfreq=-1,normalized=True):
    area = 0.
    if(len(freqs) == 0 or len(merits) == 0):
        return area
    for k in range(0,len(freqs)-1):
        if maxfreq >0. and freqs[k]<= maxfreq:
            y1 = merits[k]
            y2 = merits[k+1]
            dx = freqs[k+1]-freqs[k]
            if(y2<y1):
                y2 = merits[k]
                y1 = merits[k+1]
            area += dx*(y1+(y2-y1)/2.)

        if maxfreq >0. and freqs[k]< maxfreq and freqs[k+1]>maxfreq:
            dx = maxfreq-freqs[k]
            y1 = merits[k]
            y2 = y1+(maxfreq-freqs[k])/(freqs[k+1]-freqs[k])*(merits[k+1]-merits[k])

            if(y2<y1):
                yswap = y2
                y2 = merits[k]
                y1 = yswap
            area += dx*(y1+(y2-y1)/2.)

    if normalized:
        if maxfreq < 0.:
            area /= freqs[-1]
        else:
            area /= maxfreq

    return area
