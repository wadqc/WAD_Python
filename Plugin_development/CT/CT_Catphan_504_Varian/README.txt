This module is completely based on pylinac CBCT: http://pylinac.readthedocs.org/
As pylinac is developed for Python3 some additional actions have to be taken to make it work:
- Install the minimum requirements for pylinac3to2 (with pip). Warning! These versions are higher than the current (22-12-2015) versions used in PyWAD, use at own risk (Numpy 1.8.2, Scipy 0.14.1):
    numpy == 1.9
    scipy == 0.15
    matplotlib == 1.3.1
    pydicom == 0.9.9
    Pillow == 2.5
- Install additional requirement: future==0.15.2	
- Download pylinac3to2 from: https://github.com/Wenze/pylinac/tree/3to2
- Install pylinac using setup.py
- Disable figure plotting for Matlab by editing 'matplotlibrc', change backend to: backend      : Agg 
	file location: start python, import matplotlib, matplotlib.matplotlib_fname(), output: '/home/foo/.config/matplotlib/matplotlibrc'
	Note: This should not be necessary when using show=False in pylinac. This works in IDE but unfortunately not in WAD

