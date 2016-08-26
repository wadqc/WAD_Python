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
import numpy as np
import os

def _importProfile(fname='profile.tsv'):
    return np.loadtxt(os.path.join('temp',fname), delimiter='\t')

def _exportProfile(data,fname='profile.tsv'):
    np.savetxt(os.path.join('temp',fname), data, delimiter='\t')

def _importNDArray(fname='ndarray.npy'):
    return np.load(os.path.join('temp',fname))

def _exportNDArray(data,fname='ndarray.npy'):
    np.save(os.path.join('temp',fname), data)

