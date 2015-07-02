#!/usr/bin/python
""" Changelog:
20140703: initial version
"""
import sys
import MySQLdb


import MySQLdb
import sys,os,shutil
from operator import itemgetter

import NEMAIQ
import NEMASUV

#from pyWAD import PluginData
#    series_lst = [['./dicom.dcm']]
#    plugin_data = PluginData(series_lst)

class get_earl_qc():

        #WAD structuur
        # 
        # 1) selector_tabel_fk -> earl_iq en earl_suv keys
        # 2) gewenste_processen tabel bevat selector_fk en geeft alle processsen. Status = 10 uit deze tabel betekent dat het afgerond is.
        # 3) De result tabellen bevatten het gewenste processen key nr.

        # 4) Door te zoeken op grootheid procstat en waarde 0 of "unprocessed" dan vinden we alle onverwerkte gewenste processen.   # 5) De benodigde dicom file vinden we uit gewenste_processen.analyse_input_fk -> filename+path -> open xml en lees source uit

        ## Select a study to process. There are 2 options:
        ## 1) If it is a NEMAIQ we run the NEMAIQ tool. 
        ## 2) If it is a SUV we run the NEMASUV tool. 


        ## Run tool and give manual input


    ## Read database and find all initialized EARL studies
    def __init__(self,iqselector='EARL_IQ',suvselector='EARL_SUV'):

        self.db = MySQLdb.Connection(host="127.0.0.1",user="iqc",passwd="TY8BqYRdn3Uhzq8T",db="iqc" )
        self.earliqlist = []
        self.earlsuvlist = []

        self.earlsuvselector = suvselector
        self.earliqselector = iqselector    

        with self.db: 
            self.cur = self.db.cursor()

            self.cur.execute("SELECT * FROM analysemodule_output")
            self.analyseoutput = self.cur.fetchall()

        self.selection={}
        self.qresults = self.findEarlQC()
        
        self.curr_results={}
        

        ## Update the database results with the 


    def earllookup(self,idx):

        print ('0000\t')
        self.cur.execute("SELECT * from gewenste_processen where analysemodule_output_fk = %s"%str(idx))

        tmpchar = int(self.cur.fetchall()[0][1])

        self.cur.execute("SELECT * from selector where pk = %s"%str(tmpchar))


        return self.cur.fetchall()[0][1]

    def print_selectors(self):
        print self.selectorrows

    

    def print_lastearl(self):
        print '-'*20
        print 'List of last QC dates:'


        for row in self.earlsuvlist + self.earliqlist:

            tmpresults = []
            for elem in self.analyseoutput:
                
                if row in elem[1]:
                    tmpresults.append(elem)

            #print 'Aanwezige resultaten:',tmpresults
            datelist = []
            for elem in tmpresults:
                datelist.append([elem[1],elem])#[int(elem[1].split('/')[1]),int(elem[0])])

            print "datelist: ", sorted(datelist) #.sort(key=itemgetter(0))
        return sorted(datelist)


    def keyfromresult(self, dictionary, value):
        print '----', dictionary
        tmpdict = {}
        for key in dictionary:
            for tmp in dictionary[key]:
                tmpdict[str(tmp[3])] = key

        
        print 'tmpdictkye', tmpdict.keys(),value, type(value) ,type(tmpdict.keys()[0])

        if value in tmpdict.keys():
            return tmpdict[value]
        return None


    def print_earliq(self):
        print 'EARL IQ results:', self.earliqlist
        return

    def print_earlsuv(self):
        print 'EARL SUV results:', self.earlsuvlist
        return
        
    def get_input_pars(self,commandline=True):

        inputdict = {}
        

        if commandline == True:
            
            print "\n Input parameters: "
            #inputdict['mode'] = raw_input("mode: ")
            #inputdict['axorder'] = raw_input("axorder:") 
            #inputdict['sporder']= raw_input("sporder:")
            inputdict['sphereac']= raw_input("sphereac: (default 0)") or 0
            inputdict['ctimespac']= raw_input("ctimespac: ") or '000000'
            inputdict['bgac']=raw_input("bgac: (default 0)") or 0
            inputdict['bgactime']=raw_input("bgactime:") or '000000'
            inputdict['spremac']=raw_input("spremac: (default 0)") or 0
            inputdict['spremactime']=raw_input("spremactime:") or '000000'
            inputdict['bgremac']=raw_input("bgremac: (default 0)") or 0
            inputdict['bgremactime']=raw_input("bgremactime:") or '000000'
            inputdict['scantime']=raw_input("scantime: ") or '000000'

        return inputdict

    def get_earliq(self,filename,mode,axorder,sporder):
        inputdict = self.get_input_pars()

        self.curr_results =  NEMAIQ.update_wad(filename, mode, axorder, sporder, **inputdict)


        #sphereac, ctimespac, bgac, bgactime, spremac, spremactime, bgremac, bgremactime ,scantime)

    def get_earlsuv(self,filename,mode,axorder,sporder):
        inputdict = self.get_input_pars()
        inputdict['volume']=9.293
        inputdict['dimensions']=[10,10]

        self.curr_results = NEMASUV.update_wad(filename, mode, axorder, sporder, **inputdict)

#bgac, dimensions = [15,5], volume = '9200', caltime=120000 , remainder=0, remtime= 120000)




        
    def write_results(self):


        print "updating float results"
        
        

        self.cur.execute("SELECT * from resultaten_floating where gewenste_processen_fk = %s"%str(self.selection['gewenste_pk']))
        tmpfloats = self.cur.fetchall()


        for elem in tmpfloats:
            print 'parameter',elem[2],'value', elem[7], 'updating to' ,self.curr_results[elem[2]]
            mysqlstr = "update resultaten_floating set waarde='%s' where (gewenste_processen_fk=%s) and (omschrijving='%s')"%(self.curr_results[elem[2]],str(self.selection['gewenste_pk']),str(elem[2]))

            try:
                self.cur.execute(mysqlstr)
            except:
                print "Failed to execute: %s"%mysqlstr

        print "updating char self.curr_results"
        self.cur.execute("SELECT * from resultaten_char where gewenste_processen_fk = %s"%str(self.selection['gewenste_pk']))
        tmpchar = self.cur.fetchall()
        for elem in tmpchar:
            print 'parameter',elem[2],'value', elem[5], 'updating to'  , self.curr_results[elem[2]]
            
            mysqlstr = "update resultaten_char set waarde='%s' where (gewenste_processen_fk=%s) and (omschrijving='%s')"%(self.curr_results[elem[2]],str(self.selection['gewenste_pk']),str(elem[2]))
            try:
                self.cur.execute(mysqlstr)
            except:
                print "Failed to execute %s"%mysqlstr

        print "updating image results"
        self.cur.execute("SELECT * from resultaten_object where gewenste_processen_fk = %s"%str(self.selection['gewenste_pk']))
        tmpfig = self.cur.fetchall()
        


        for elem in tmpfig:
            if elem[2] == 'processor.log':
                pass
            else:

                print 'parameter',elem[2],'value', elem[5], 'updating to' , self.curr_results[elem[2]]
            
                self.curr_results[elem[2]].savefig(elem[5])


        print "setting status to processed"
    
        self.cur.execute("update resultaten_char set waarde='processed' where (gewenste_processen_fk=%s) and omschrijving='status'"%str(self.selection['gewenste_pk']))
        self.db.commit()




    def findEarlQC(self):

        print '-'*20
        print 'The following earl results have been found:'

        results = {}
        

        for row in [self.earlsuvselector,self.earliqselector]:
            
            results[row] = []

            tmpresults = []
            for elem in self.analyseoutput:
                if row in elem[1]:
                    tmpresults.append(elem)

            #print 'Aanwezige resultaten:',tmpresults

            datelist = []
            for elem in tmpresults:
                datelist.append([int(elem[1].split('/')[1]),int(elem[0])])

            #print "datelist: ", sorted(datelist) #.sort(key=itemgetter(0))
                
            for elem in sorted(datelist):
                tmppk =  elem[1]
                self.cur.execute("SELECT * from gewenste_processen where analysemodule_output_fk = %s"%tmppk)
                output = self.cur.fetchall()[0]
                


                self.cur.execute("SELECT * from resultaten_char where gewenste_processen_fk = %s"%output[0])
                limits = self.cur.fetchall()

                
                for subelem in limits:
                    if subelem[2] == 'status':
                        results[row].append([output[5],subelem[2],subelem[5],subelem[1]])

        print '--RESULTATEN--:'
        for key in results.keys():
            print key
            for nr, subresult in enumerate(results[key]):
                 print '\t',nr, subresult


        return results
        '''

        # we select a result
        print 'Voer het gewenste procesnr in:'
        gewenste_pk = raw_input()
        print 'Je voerde in %s, indien correct voer 0 in'%gewenste_pk

        tmpstatus = raw_input()
        if tmpstatus == '0':
            pass
        else:
            sys.exit()
        
        self.selection['gewenste_pk'] = gewenste_pk

        print "--@@--"
        print 'gewenste_pk',gewenste_pk
        print 'resultaten',results
        print self.keyfromresult(results, gewenste_pk)

        self.selection['type'] =  self.keyfromresult(results, gewenste_pk) #'EARL IQ' #replace with inverese dict based on selection




        # from the gewenste processen key we obtain the analysemodule_input key
        self.cur.execute("SELECT * from gewenste_processen where pk = %s"%str(gewenste_pk))
        tmpoutkey = self.cur.fetchall()[0][8]
        print 'tmpselectorkey',tmpoutkey
        #from the analyse_output key we obtain the file

        self.cur.execute("SELECT * from analysemodule_input where pk = %s"%str(tmpoutkey))
        tmpoutput = self.cur.fetchall()[0]

        #set the base dir

        tmpxml = os.path.join('/opt',tmpoutput[2],tmpoutput[1])
        import xml.etree.ElementTree as ET

        tree = ET.parse(tmpxml)
        root = tree.getroot()

        pat = tree.findall('patient')[0]
        stud = pat.findall('study')[0]
        ser = stud.findall('series')[0]
        ins = ser.findall('instance')[0]
        fil = ins.getchildren()[1].text
            


#,self.cur.fetchall()[0][2])
        print 'resultaat voor keuze', fil

        seriesdir = os.path.split(fil)[0]

        #now that we have the series we can launch the tool
        #first we ask the necessary pars

        m = 'max'
        order = 'zxy'
        ort =  'acw'


        print 'Selection type:', self.selection['type']
            
        if self.selection['type'] == self.earliqselector:
            self.get_earliq(filename=seriesdir,mode=m,axorder=order,sporder=ort)


        elif self.selection['type'] == self.earlsuvselector:
            self.get_earlsuv(filename=seriesdir,mode=m,axorder=order,sporder=ort)


        self.write_results()
        '''

def main():
    test = get_earl_qc()

   # test.print_earliq()
   # test.print_earlsuv()

if __name__ == "__main__":
    sys.exit(main())

