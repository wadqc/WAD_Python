#!/usr/bin/env python
import sys
import os.path as path
import os
import zipfile

# name of the zipfile to be created; this name should be <root-plugin.py>.zip
outfile = 'pywadplugin.py.zip'

# First rule: exclude these directories completely
excludedirs = [
    '.git',
    'Tools',
    'Testing',
    'Documentation',
    '1', 'wpr', # some personal rubbish
]

# Second rule: exclude all filenames that match these exactly (in any folder)
excludefiles = [
    'error.log'
]

# Third rule: only include filenames that match these extensions (in any folder)
includeextensions = [
    '.py','.exe'
]

class SimpleArgs():
    """
    Class to handle unsophisticated commandline arguments, and to isolate it from the main programme
    """
    # define the allowed modes
    modes = ['list','zip']

    def __init__(self):
        pass

    def usage(self,arglist=[]):
        """
        display the usage of this file and exit
        """
        print "Usage: %s <command>" %arglist[0]
        print "  <command> is one of:",
        for mo in self.modes:
            print mo,
        print ""
        if len(arglist) > 1:
            print "%s is not a valid command."%arglist[1]
        sys.exit()
    
    def getmode(self,arglist):
        # Get the arguments list 
        mode = self.modes[0]
        
        if len(arglist) != 2 :
            self.usage(arglist)
        else:
            if arglist[1] in self.modes:
                mode = arglist[1]
            else:
                usage(arglist)
        return mode

#------------
## main
mode = SimpleArgs().getmode(sys.argv)

# make sure the working dir is the Tools folder in pywad2
os.chdir(os.path.dirname(__file__)) 
# change to parent dir
os.chdir('..') 
rootdir = os.getcwd()

# create zipfile
if mode == 'zip':
    zf = zipfile.ZipFile(outfile, 'w',zipfile.ZIP_DEFLATED)
    
# find all files that match the given criteria
count = 0
for subdir, dirs, files in os.walk(rootdir,topdown=True):
    dirs[:] = [d for d in dirs if d not in excludedirs]
    files[:] = [f for f in files if f not in excludefiles]
    for fname in files:
        for ext in includeextensions:
            if fname.endswith(ext):
                filename = os.path.join(subdir, fname)
                daf = os.path.relpath(filename, start=rootdir)
                if mode == 'zip':
                    zf.write(filename,daf)
                print daf
                count += 1
                    
# close zipfile
if mode == 'zip':
    print "Done. Packed %d files in %s"%(count,os.path.join(rootdir,outfile))
    zf.close()
else:
    print "%d files meet the criteria. Just listing, no zip file created."%count
