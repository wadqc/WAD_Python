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
    def __init__(self):

        self.db = MySQLdb.Connection(host="127.0.0.1",user="iqc",passwd="TY8BqYRdn3Uhzq8T",db="iqc" )
        self.earliqlist = []
        self.earlsuvlist = []
        

        with self.db: 
            self.cur = self.db.cursor()

            self.cur.execute("SELECT * FROM analysemodule_output")
            self.analyseoutput = self.cur.fetchall()

        self.get_earlsuv()
        self.get_earliq()


## Update the database results with the 


    def print_selectors(self):
        print self.selectorrows

    def print_lastearl(self):
        print '-'*20
        print 'List of last QC dates:'

        for row in ['EARL SUV','EARL IQ']:

            tmpresults = []
            for elem in self.analyseoutput:
                
                if row in elem[1]:
                    tmpresults.append(elem)

            #print 'Aanwezige resultaten:',tmpresults
            datelist = []
            for elem in tmpresults:
                datelist.append([int(elem[1].split('/')[1]),int(elem[0])])

            print "datelist: ", sorted(datelist) #.sort(key=itemgetter(0))


    def print_earliq(self):
        print 'EARL IQ results:', self.earliqlist
        return

    def print_earlsuv(self):
        print 'EARL SUV results:', self.earlsuvlist
        return

    def get_earlsuv(self):
        return

    def get_earliq(self):

        print '-'*20
        print 'The following earl results have been found:'

        results = {}
        

        for row in ['NM_EARL_SUV','NM_EARL_IQ']:
            
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

        selection = {}

        


        # we select a result
        #print 'selecting 2014-7-5 15:48:47 gewenste_pk = 30'
        print 'Voer het gewenste procesnr in:'
        gewenste_pk = raw_input()
        print 'Je voerde in %s, indien correct voer 0 in'%gewenste_pk

        tmpstatus = raw_input()
        if tmpstatus == '0':
            pass
        else:
            sys.exit()

        
        selection['gewenste_pk'] = gewenste_pk
        selection['type'] = 'NM_EARL_SUV'

        # from the gewenste processen key we obtain the analysemodule_input key
        self.cur.execute("SELECT * from gewenste_processen where pk = %s"%str(gewenste_pk))
        tmpoutkey = self.cur.fetchall()[0][8]
        print 'tmpselectorkey',tmpoutkey
        #from the analyse_output key we obtain the file

        self.cur.execute("SELECT * from analysemodule_input where pk = %s"%str(tmpoutkey))
        tmpoutput = self.cur.fetchall()[0]

        #set the base dir

        tmpxml = os.path.join('D:\WAD-software',tmpoutput[2],tmpoutput[1])
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
        bga = 10.0
        spa = 1.0

            
        if selection['type'] == 'NM_EARL_IQ':
                results['EARL IQ'] = NEMAIQ.update_wad(filename = seriesdir, mode = m, axorder = order, sporder = ort,bgac = bga, sphereac = spa)
        elif selection['type'] == 'NM_EARL_SUV':
                results['EARL SUV'] = NEMASUV.update_wad(filename = seriesdir, mode = m, axorder = order, sporder = ort,bgac = bga, dimensions = [15,5], volume = '9200', caltime=120000 , remainder=0, remtime= 120000)

        print "write to database"                
        for key in results.keys():
            print key, results[key]

        print "updating float results"
        self.cur.execute("SELECT * from resultaten_floating where gewenste_processen_fk = %s"%str(selection['gewenste_pk']))
        tmpfloats = self.cur.fetchall()
        for elem in tmpfloats:
            print 'parameter',elem[2],'value', elem[7], 'updating to' , results[elem[2]]


        print "updating char results"
        self.cur.execute("SELECT * from resultaten_char where gewenste_processen_fk = %s"%str(selection['gewenste_pk']))
        tmpchar = self.cur.fetchall()
        for elem in tmpchar:
            print 'parameter',elem[2],'value', elem[5], 'updating to' , results[elem[2]]


        print "updating image results"
        self.cur.execute("SELECT * from resultaten_object where gewenste_processen_fk = %s"%str(selection['gewenste_pk']))
        tmpfig = self.cur.fetchall()
        
        for elem in tmpfig:
            print 'parameter',elem[2],'value', elem[5], 'updating to', results[elem[2]]


        print "setting status to processed"
        self.cur.execute("update resultaten_char set waarde='processed' where gewenste_processen_fk=%s"%str(gewenste_pk))
        self.db.commit()
    


def main():
    test = get_earl_qc()
    test.print_earliq()
    test.print_earlsuv()

if __name__ == "__main__":
    sys.exit(main())

