from distutils.core import setup
import py2exe
import sys, os

sys.path.append("source")
sys.argv.append('py2exe')

opts = {
    'build': {'build_base': "build/"},
    'py2exe': {"dist_dir": "bin/",
               "compressed": True,
               "optimize": 2,
               "bundle_files" : 1,
               "includes" :[],
               "excludes": ['_gtkagg', '_tkagg', '_agg2', 
                            '_cairo', '_cocoaagg',
                            '_fltkagg', '_gtk', '_gtkcairo', 'tcl',
                            'Tkconstants', 'Tkinter', 'numpy'],
               "dll_excludes": ['w9xpopen.exe'],
              }
       }

setup(console=["source/File2DicomSend.py"],
      options=opts,
      zipfile = None,
)
