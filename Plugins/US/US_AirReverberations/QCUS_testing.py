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

