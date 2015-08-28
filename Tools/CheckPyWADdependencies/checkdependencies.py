import os, re

importdict = {}

for dirname,  dirnames,  filenames in os.walk('../Plugins'):
#    for subdirname in dirnames:
#        print os.path.join(dirname, subdirname)
    for filename in filenames:
        if os.path.splitext(filename)[-1] == '.py':
            if filename == '__init__.py':
                pass
            else:
                print filename

                tmpfile = os.path.join(dirname, filename)
                text = open(tmpfile).readlines()

                for line in text:
                    if len(line.split()) > 0:
                        
                        # There are three import possibilities:
                        # 1) import A
                        # 2) import A,B
                        # 3) from A import B

                        match1 = '^(\s)+import+(\s)(\w)'
                        match2 = '^(\s)+from+(\s)(\w)(\s)+import+(\s)(\w)'
                        match3 = '^(\s)+import+(\s)(\w)(\s)(,(\s)(\w))*'
                        match4 = '^(\s)+from+(\s)(\w)(\s)+import+(\s)(\w)(,(\s)(\w))*'


                        if line.split()[0] ==  'import':
                            #print 'import: ',line

                            for elem in line.split(','):
                                tmpmod = elem.split()[0]
                                
                                if tmpmod in importdict.keys():
                                    importdict[tmpmod].append(filename)
                                else:
                                    importdict[tmpmod] = [filename]
                            

for key in importdict.keys():
    print key



