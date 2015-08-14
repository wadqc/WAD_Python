import os
import sys
import re
from itertools import chain

importdict = {}
package_dirs  = ('../../pyWAD','../../pyWADLib','../../Plugins','../../Tools')
localfiles = []

for dirname,  dirnames,  filenames in chain.from_iterable(os.walk(path) for path in package_dirs):
    for filename in filenames:
        if os.path.splitext(filename)[-1] == '.py':
            if filename == '__init__.py':
                pass
            elif filename == 'CheckDependencies.py':
                pass
            else:
                tmpfile = os.path.join(dirname, filename)
                localfiles.append(os.path.splitext(filename)[0])
                text = open(tmpfile).readlines()

                for line in text:
                    if len(line.split()) > 0:

                        match1 = '^(\s)*import(\s)(\w)'
                        match2 = '^(\s)*from(\s)*(\w)(\s)*import'

                        def directimport(line):
                            for tmpmod in [elem.split('as')[0].strip() for elem in line.split('import')[1].split(',')]:

                                
                                if tmpmod in importdict.keys():
                                    importdict[tmpmod].append((filename,dirname))
                                else:
                                    importdict[tmpmod] = [(filename,dirname)]


                        def fromimport(line):
                            for tmpmod in [elem.split('as')[0].strip() for elem in line.split('import')[1].split(',')]:

                                
                                if tmpmod in importdict.keys():
                                    importdict[tmpmod].append((filename,dirname))
                                else:
                                    importdict[tmpmod] = [(filename,dirname)]

                        if re.match(match1,line):                        
                            directimport(line)
                        elif re.match(match2,line):
                            fromimport(line)

import importlib
for key in importdict.keys():
    if key not in localfiles:
        
        try:
            importlib.import_module(key)
        except:
            print 'Warning could not load',key, 'calling routine:', importdict[key]





