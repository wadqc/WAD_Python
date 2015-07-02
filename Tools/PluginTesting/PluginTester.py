#!/usr/bin/env python
import os,sys
import ConfigParser
import dicom
from datetime import datetime
from loaddicom import LoadDicom
import subprocess
import datetime
import os
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import argparse
import math


# CHANGELOG:
#  20150626: added sanity check; fixed instance level

def instancexml(number,filename):
    # returns a xml-like string containing an series object
    return " <instance><number>%s</number><filename>%s</filename></instance>"%(number,filename)

def seriesxml(number,description,instance):
    # returns a xml-like string containing an series object
    desc = description.replace('<','').replace('>','')
    return "<series><number>%s</number><description>%s</description>%s</series>"%(number,desc,instance)

def inputxml(configpath,outputpath,idpar,name,uid,series,level='series'):
    # returns a input.xml skeleton as a xml-like string
    xml = '''<?xml version="1.0" encoding="UTF-8"?><WAD>
           <version>
           </version>
           <analysemodule_cfg>%s</analysemodule_cfg>
           <analysemodule_output>%s</analysemodule_output>
           <analyselevel>%s</analyselevel>
           <patient>
              <id>%s</id>
              <name>%s</name>
              <study>
                 <uid>%s</uid>
                 <description>
                 </description>
                   %s
              </study>
           </patient>
        </WAD>
        '''%(configpath,outputpath,level,idpar,name,uid,series)
    return xml


def createInputXML(configxml,datapath,outputdir,level):
    '''
    From a study hierarchy an inputXML file is created. This file can be used directly to test the pywad framework using
    python pywadplugin.py <inputXML>.
    '''


    tmpoutputdir = outputdir

    serieslst = LoadDicom(datapath)
    inputxmllst = []

    if level=='study':
        seriesxmllst = []
        for key in serieslst.keys():
            tmpseries = serieslst[key]
            tmpinstance = ''.join([instancexml(*instance) for instance in tmpseries['filelst']])
            seriesxmllst.append(seriesxml('1',key,tmpinstance)) #to do add correct seriesnumber
        _seriesxml = ''.join(seriesxmllst)
            
        inputxmllst.append(inputxml(configxml,os.path.join(tmpoutputdir,'result.xml'),'idpar','name','uid',_seriesxml,level=level))

    elif level=='series':
        for key in serieslst.keys():
            tmpseries = serieslst[key]
            tmpinstance = ''.join([instancexml(*instance) for instance in tmpseries['filelst']])
            tmpseriesxml =  seriesxml('1',key,tmpinstance) #to do add correct seriesnumber
            inputxmllst.append(inputxml(configxml,os.path.join(tmpoutputdir,'result.xml'),'idpar','name','uid',tmpseriesxml,level=level))
        
    elif level=='instance':
        for key in serieslst.keys():
            tmpseries = serieslst[key]
            for instance in tmpseries['filelst']:
                tmpseriesxml =  seriesxml('1',key,instancexml(*instance)) #to do add correct seriesnumber
                inputxmllst.append(inputxml(configxml,os.path.join(tmpoutputdir,'result.xml'),'idpar','name','uid',tmpseriesxml,level=level))
        
    return inputxmllst


def _determinelevel(configfile):
    '''
    determine analysis_level from location of configfile:
    if it is in a subfolder 'series' or 'study' or 'instance', then that is the 
    analysis_level. Else, assume it is 'series'
    '''
    # /some/folder/series/config_x.xml
    levelfolder = os.path.split(os.path.split(configfile)[0])[1]
    level = 'series' if not levelfolder in ['study','instance'] else levelfolder
    return level

def _testlist_personal(pywadpath,dataroot,configdict):
    # returns a list of user tests, based on defined plugins

    errorlist = {'udata':[],'uconfig':[],'uinput':[]}
    testlist = {} # { plugin: {id: {input:[], config:[], data:[]} } }
    print '\nBuilding list of user tests...'
    for plugin,settings in configdict.iteritems():
        if plugin == 'INIT':
            continue

        # /some/folder/series/config_x.xml
        level = settings['level'] if 'level' in settings else _determinelevel(settings['config'])
            
        if not level in ['study','instance','series']:
            errorlist['ulevel'].append('%s cannot determine analysis_level'%levelfolder)
            print 'ERROR',errorlist['udata'][-1]
            continue

        if not os.path.exists(os.path.abspath(settings['data'])):
            errorlist['udata'].append('%s does not exist'%os.path.abspath(settings['data']))
            print 'ERROR',errorlist['udata'][-1]
            continue

        if not os.path.exists(os.path.abspath(settings['config'])):
            errorlist['uconfig'].append('%s does not exist'%os.path.abspath(settings['config']))
            print 'ERROR',errorlist['uconfig'][-1]
            continue

        testid = os.path.splitext(os.path.split(settings['config'])[1])[0]
        if not plugin in testlist:
            testlist[plugin] = {}
        if not testid in testlist[plugin]:
            testlist[plugin][testid] = {'config': [], 'data': [], 'level':[]}

        testlist[plugin][testid]['config'].append(settings['config'])
        testlist[plugin][testid]['level'].append(level)
        testlist[plugin][testid]['data'].append(settings['data'])
        
    return testlist,errorlist
    
def _testlist_default(pywadpath,dataroot,pluginroot='Plugins'):
    # returns a list of tests to perform for default plugins

    print '\nBuilding list of %s tests...'%pluginroot
    if not pluginroot.endswith('/'): pluginroot += '/'
    plugindir = os.path.join(pywadpath,pluginroot) # note the trailing '/'
    # traverse root directory, and list directories as dirs and files as files
    plugins = []
    for root, dirs, files in os.walk(plugindir):
        subdir = root.split(plugindir)[1]
        #print (len(path) - 1) *'---' , os.path.basename(root)       
        for f in files:
            if f.endswith('.py') and f != '__init__.py':
                #print subdir, f ,os.path.join(plugindir,subdir)
                plugins.append(subdir)
    plugins = list(set(plugins))

    # check if Testing folder structure is in accordance with plugins:
    print '\nChecking folder consistency of Testing and Plugins for %s...'%pluginroot
    testingdir = os.path.join(pywadpath,'Testing') 
    errorlist = {'testing':[]}
    testlist = {} # { plugin: {id: {config:[], data:[], 'level':[]} } }
    for plugin in plugins:
        testpath = os.path.join(testingdir,plugin)
        if not os.path.exists(testpath):
            errorlist['testing'].append('%s does not exist'%testpath)
            continue
        checkpaths = [testpath,
                      os.path.join(testpath,'series'),
                      os.path.join(testpath,'study'),
                      os.path.join(testpath,'instance'),
                      ]
        for cp in checkpaths:
            if os.path.exists(cp):
                for f in os.listdir(os.path.abspath(cp)):
                    if f.endswith('.xml') and f.startswith('config_'):
                        testid = f.replace('.xml','').replace('config_','')
                        if not plugin in testlist:
                            testlist[plugin] = {}
                        if not testid in testlist[plugin]:
                            testlist[plugin][testid] = {'config': [], 'data': [], 'level':[]}
        
                        testlist[plugin][testid]['config'].append(os.path.join(cp,f))
                        testlist[plugin][testid]['level'].append(_determinelevel(testlist[plugin][testid]['config'][-1]))

        if not plugin in testlist:
            errorlist['testing'].append('%s has no corresponding config_.xml files'%(testpath))
                
    for e in errorlist['testing']:
        print 'ERROR',e
        
    print '\nChecking folder consistency of Testing and Data for %s...'%pluginroot
    errorlist['data'] = []
    # check data folder structure
    droplist = []
    for plugin,testdict in testlist.iteritems():
        for test,bla in testdict.iteritems():
            testpath = os.path.join(os.path.join(dataroot,plugin),'dicom_%s'%test)
            if not os.path.exists(testpath):
                errorlist['data'].append('%s does not exist'%testpath)
            elif len(os.listdir(testpath)) == 0:
                errorlist['data'].append('%s is empty'%testpath)
            else:
                testlist[plugin][test]['data'].append(testpath)
        if len(testlist[plugin][test]['data']) == 0:
            del testlist[plugin][test]
        if len(testlist[plugin]) == 0:
            droplist.append(plugin)

    testlist = { k : v for k,v in sorted(testlist.iteritems()) if not k in droplist }

    for e in errorlist['data']:
        print 'ERROR',e

            
    return testlist,errorlist

def sanitycheck(args):
    '''
    reguires a folder structure 
       <root>/Plugins/<some path>/<plugin name>
       <root>/Testing/<some path>/<plugin name>
    The Testing plugin folder should at least one of:
       series/config_<id>.xml (plugin configuration file)
       instance/config_<id>.xml (plugin configuration file)
       study/config_<id>.xml (plugin configuration file)
    In addition there needs to be a <testdata> folder with a folder dicom_<id>, in which test data is present for this plugin.

    Example: 
      plugin: <root>/Plugins/CT/CT_Philips_QuickIQ
      plugin files: ....
      testing: <root>/Testing/CT/CT_Philips_QuickIQ
      testing files: series/config_body.xml,series/config_head.xml
      testdata: <dataroot>/CT/CT_Philips_QuickIQ/dicom_head/ and <dataroot>/CT/CT_Philips_QuickIQ/dicom_body/
    '''
    #Read config
    config = ConfigParser.ConfigParser()
    config.read(args.inifile)

    sections = config.sections()
    configdict = {}
    for section in sections:
        options = config.options(section)
        configdict[section] = dict([(option,config.get(section,option)) for option in options])
    
    pywadpath = args.packagedir if not args.packagedir is None else configdict['INIT']['packagedir']
    exename = configdict['INIT']['exe']
    dataroot = args.dataroot if not args.dataroot is None else configdict['INIT']['dataroot']
    skipdefault = config.getboolean('INIT','skipdefault') if config.has_option('INIT','skipdefault') else False
    add_development = config.getboolean('INIT','adddevelopment') if config.has_option('INIT','adddevelopment') else False

    testlist,errorlist = {},{}
    if not skipdefault:
        testlist_default,errorlist_default = _testlist_default(pywadpath, dataroot,'Plugins')
        testlist = dict(testlist,**testlist_default)
        errorlist = dict(errorlist,**errorlist_default)
        
    if add_development:
        testlist_default,errorlist_default = _testlist_default(pywadpath, dataroot,'Plugin_development')
        testlist = dict(testlist,**testlist_default)
        errorlist = dict(errorlist,**errorlist_default)

    if len(sections) > 1:
        testlist_user,errorlist_user = _testlist_personal(pywadpath,dataroot,configdict)
        testlist = dict(testlist,**testlist_user)
        errorlist = dict(errorlist,**errorlist_user)
        
    # create an output folder with as a name the current date
    outputdir = os.path.abspath(os.path.join(configdict['INIT']['outputdir'],datetime.date.isoformat(datetime.date.today())))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    numplug,numtests = 0,0
    for plugin,testdict in testlist.iteritems():
        numplug += 1
        numtests += len(testdict.keys())
            
    print ('\nAttempting to run %d plugins, with %d tests...'%(numplug,numtests))

    #generate a batch of inputfiles
    results,errorlist['xml'] = batchtest(exename,pywadpath, testlist, outputdir)
    
    for e in errorlist['xml']:
        print 'ERROR',e
    
    
    with open(os.path.join(outputdir,'codetest.csv'),'w') as f:
        f.write('input,returncode,message\n')
        for elem in results:
            f.write('%s,%d,"%s"\n'%(elem['input'],elem['returncode'],str(elem['err']).replace('"','').replace('\n', ' ')))

        f.write('\nERROR,type\n')
        for key,errors in sorted(errorlist.iteritems()):
            for error in sorted(errors):
                f.write('"%s","%s"\n'%(error,key))
            
    print '\noutput in ',outputdir
    
    
def batchtest(exe, exepath, testlist, outputdir):
    '''
    Tries to make proper input files for all tests specified in testlist,
    tries to run the exe on each created input file
    returns results,errorlist
     exe = 'pywadplugin' # name of exe
     exepath = '/folder/pywad2' # path to folder with exe
     testlist = {} # { plugin: {id: {input:[], config:[], data:[]} } }
     results = [] # list of results for each test
     errorlist = [] # list of errors creating XML
    '''
    batchlist = []
    errorlist = []
    for plugin,testdict in sorted(testlist.iteritems()):
        for test,bla in sorted(testdict.iteritems()):
            for conf,dat,level in zip(bla['config'],bla['data'],bla['level']):
                tmpinpxml= createInputXML(os.path.abspath(conf),os.path.abspath(dat),'replace_me_dummy',level)
                if len(tmpinpxml) == 0:
                    errorlist.append('no valid data for creation of input file for %s/%s'%(plugin,test))
                    continue
                
                for ix,elem in enumerate(tmpinpxml):
                    xmldir = os.path.join(outputdir,plugin,level,'%s.%s'%(test,('%d'%ix).zfill(3)))
                    if not os.path.exists(xmldir):
                        os.makedirs(xmldir)
                    elem = elem.replace('replace_me_dummy',xmldir)
                    xmlout = os.path.join(xmldir,'input.xml')
                    with open(xmlout,'w') as f:
                        f.write(elem)
                        batchlist.append({'plugin':plugin, 
                                          'test':test, 
                                          'config':conf, 
                                          'level':level, 
                                          'input':os.path.abspath(xmlout)})
    

    # now run the plugin for each created input file
    numprocs = len(batchlist)
    digits = 1 if numprocs<10 else int(math.ceil(math.log(numprocs,10)+.5))
    
    print ('Using %d datasets for testing...'%numprocs)
    results = []
    for ix,proc in enumerate(batchlist):
        print '%s/%d %s: %s (%s)'%(str(ix+1).zfill(digits),numprocs,proc['plugin'],proc['test'],proc['level']),
        out,err,returncode = None,None,1
        p = subprocess.Popen(['python',os.path.abspath(os.path.join(exepath,exe)),proc['input']],
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=os.path.abspath(exepath))

        out, err =  p.communicate()
        returncode = p.returncode
        if returncode==0 and err: # force runtimewarnings and other stderr output to result in errors!
            returncode = 999

        if returncode == 0:
            print 'OK'
        elif returncode == 999:
            print 'WARNING: ',err
        else:
            print 'ERROR: ',err

        if args.verbose:
            print out
            
        output = proc.copy()
        #output['out'] = out # to be added later for more ingenious checks

        output['err'] = err
        output['returncode'] = returncode
        results.append(output)
    return results,errorlist
   
if __name__ == "__main__":
    current_uid = os.getuid()
    wanted_user = 'www-data'
        
    _inifile = 'testconfig_default.ini'
    _dataroot = None
    _packagedir = None
    parser = argparse.ArgumentParser(description='Tester for pyWAD2 modules')
    parser.add_argument('-i','--inifile',
                        default=_inifile,type=str,
                        help='the configuration file to be used [%s]'%_inifile,dest='inifile')
    parser.add_argument('-d','--data',
                        default=None,type=str,
                        help='replace the dataroot in the inifile value with this folder [%s]'%str(_dataroot),dest='dataroot')
    parser.add_argument('-p','--package',
                        default=None,type=str,
                        help='replace the pywad2 packagedir in the inifile value with this folder [%s]'%str(_packagedir),dest='packagedir')
    parser.add_argument('-v','--verbose',
                        default=None,type=bool,
                        help='Verbose')   

    args = parser.parse_args()

    msg = ''
    try:
        msg += "For proper testing, this python script should be run as user %s\n"%wanted_user
        msg += "and the output folder defined in the inifile %s should be writable by user %s\n"%(args.inifile,wanted_user)
        msg += "  sudo chown -R www-data TestingOutput\n"
        msg += "  sudo -u %s python PluginTester.py\n"%wanted_user
        import pwd
        if pwd.getpwnam(wanted_user).pw_uid == current_uid:
            msg = 'Good! Running as %s'%wanted_user
    except:
        msg += 'However, user %s does not exist on this system...\n'%wanted_user

    print msg
 
    sanitycheck(args)

