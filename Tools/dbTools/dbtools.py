#dbtools.py
"""
Module for query/maintenance of the dcm4chee and icq databases.

Changelog:
20150617: Fixed crash on empty value
20150408: Status string clean up
20150401: Ignore pps_start, always prefer series_date; option to ignore pps_start, which makes a difference for MR
20150320: Added num_expected to .period
20150318: Added status query (now same functionality as in dbTools.py)

TODO:
"""
import datetime
import time
from dateutil.relativedelta import relativedelta
import MySQLdb as mdb
from operator import attrgetter
from ConfigParser import SafeConfigParser

qcversion = 20150617

class Lit:
    stNone = 'NONE'
    stResultOK = 'OK'
    stResultERROR = 'ERROR!'
    stResultWARNING = 'WARNING!'
    stResultLATE = 'LATE'
    
    db_iqc = 'iqc'

    table_gewenste_processen = 'gewenste_processen'
    table_studies = 'study'
    table_series = 'series'
    table_instances = 'instance'
    table_patient = 'patient'
    table_files = 'files'
    table_selector = 'selector'
    table_status_omschrijving = 'status_omschrijving'
    
    table_char = 'resultaten_char'
    table_float = 'resultaten_floating'
    table_object = 'resultaten_object'
    table_bool = 'resultaten_boolean'
    table_status = 'resultaten_status'
    
    table_analyseinput = 'analysemodule_input'
    table_analyseoutput = 'analysemodule_output'
    table_collectorstatusseries  = 'collector_series_status'
    table_collectorstatusstudies = 'collector_study_status'

class QCStatus:
    """
    Class for qcstatus of a machine
    """
    def __init__(self,machine,ininame=None):
        self.machine = machine
        self.selectors = []
        self.qcIgnores = []
        self.modality = []
        self.frequency = []
        self.status_short = []
        self.status_full = []
        self.date = None
        self.period = None # ('frequency',oldestdate,num_done,num_expected)

        self.overdue = False
        self.production = 0.
        
        if not ininame is None:
            self.configRead(ininame)
        
    def getHeaders(self):
        """
        Returns a list of names of entries and corresponding attribute
        """
        headers = [
            ('Machine','machine'),
            ('Modality','modality'),
            ('DateTime','date'),
            ('ShortStatus','status_short'),
            ('Overdue','overdue'),
            ('Frequency','frequency'),
            ('Production%','production')
        ]
        return headers

    def configError(self,ininame):
        """
        Displays an error reading the given config file.
        """
        error = "ERROR: configfile %s corrupt or entry non-existent. Should look like:\n[%s]\nModality = CT\n" \
              "Selectors = Selector1|Selector2\nFrequency=bi-weekly\nIgnoreDescriptions=name1|name2\n\n"%(ininame,self.machine)
        raise ValueError("[configRead] Error: %s" % (error))
    
    def configRead(self,ininame):
        """Reads given config """
        parser = SafeConfigParser()
        parser.read(ininame)
        
        # check if machine is defined in config
        section = self.machine
        if not parser.has_section(section):
            self.configError(ininame)

        now = datetime.datetime.now()
        frequencies = [
            ('weekly', now-datetime.timedelta(days=7)),
            ('bi-weekly',now-datetime.timedelta(days=14)),
            ('daily',now-datetime.timedelta(days=1)),
            ('monthly',now-relativedelta(months=1))
        ]
        self.period = ('semi-anual',now-relativedelta(months=6))
        
        # check if all info of this machine is available
        candidates = ['Modality','Selectors','Frequency','IgnoreDescriptions']
        for can in candidates:
            if not parser.has_option(self.machine,can):
                configError(ininame)

        # build the list of selectors
        self.selectors = parser.get(self.machine,'Selectors').split('|')

        # read the modality
        self.modality = parser.get(self.machine,'Modality')

        # build the list of ignores
        qcIgnores = parser.get(self.machine,'IgnoreDescriptions').split('|')
        if len(qcIgnores) == 1 and qcIgnores[0] == '':
            qcIgnores = []

        # read the frequency of testing
        freq = parser.get(self.machine,'Frequency')
        try:
            lf = freq.lower()
            found = False
            for ff in frequencies:
                if lf == ff[0]:
                    self.frequency = ff
                    found = True
            if not found:
                raise ValueError("[configRead] Error: frequency %s not recognized" % (freq))
        except Exception,e:
            print e

class Entry:
    """
    Class for entries in the iqc and dcm4chee database
    """
    def __init__(self,seriesid=None,studyid=None,stationname=None,srcaet=None,modality=None,date=None,description=None,gpid=None,status=None,note=None):
        self.gpid = gpid
        self.stationname = stationname
        self.modality = modality
        self.date = date
        self.seriesid = seriesid
        self.studyid = studyid
        self.srcaet = srcaet
        self.description = description
        self.label = None # identifier for _head or _table or _...
        self.status = status
        self.qcstatus = None
        self.note = note

    def getHeaders(self):
        """
        Returns a list of names of entries and corresponding attribute
        """
        headers = [
            ('SeriesId','seriesid'),
            ('StudyId','studyid'),
            ('DateTime','date'),
            ('Modality','modality'),
            ('StationName','stationname'),
            ('SrcAET','srcaet'),
            ('Description','description'),
            ('Label','label'),
            ('GewensteProcesId','gpid'),
            ('Status','status'),
            ('QCStatus','qcstatus'),
            ('Note','note'),
        ]
        return headers

class DBIO:
    """
    Class with all lower level sql actions
    """
    def __init__(self,host,user,pswd):
        """
        Needs hostname, username and pswd of the mysql database to be used later
        """
        self.host = host
        self.user = user
        self.pswd = pswd
        self.ignore_pps_start = False

    def ignorePPS(self,val):
        self.ignore_pps_start = val

    def connectdb(self,db,host=None,user=None,pswd=None):
        """
        open a connection to given database db
        """
        if host is None:
            host = self.host
        if user is None:
            user = self.user
        if pswd is None:
            pswd = self.pswd

        try:
            con = mdb.connect(host,user,pswd,db) # host, username, pswd, db
            cur = con.cursor(mdb.cursors.DictCursor) # dictionary format, i.e. access by column name
            dbversion='0'
            if db == Lit.db_iqc:
                cur.execute("select * from config")
                rows = cur.fetchall()
                for row in rows:
                    if row['property'] == "Version_Database":
                        dbversion = str(row['value'])
        except mdb.Error, e:
            raise ValueError("[connectdb] Error %d: %s" % (e.args[0],e.args[1]))

        return con,cur,dbversion

    def closedb(self,con):
        """
        Porperly close db connection
        """
        if con:
            con.close()

    def _labelsFromGP(self,db):
        gp_label ={}

        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[_labelsFromGP] Error: %s" % (e.args[0]))

        labeldescriptions = [ 'Anatomy', 'stand_', 'FilterMaterialLT']
        cur.execute("SELECT * FROM %s" % (Lit.table_char))
        rows = cur.fetchall()
        for row in rows:
            gpid = row['gewenste_processen_fk']
            if gp_label.has_key(gpid):
                continue
            oms = row['omschrijving']
            ix = [ i for i,elem in enumerate(labeldescriptions) if elem in oms ]
            if len(ix)>0:
                gp_label[gpid] = row['waarde']

        ### close db
        self.closedb(con)
        
        return gp_label

        
    def addLabelToSeries(self,db,entries,gp_label):
        # If an extra identifier (like Filter or Antomy) is given use that as a label

        for entry in entries:
            if gp_label.has_key(entry.gpid):
                entry.label = gp_label[entry.gpid]
        return entries
    
    def splitOnLabels(self,entries):
        """
        If entries have labels, split in dictionary
        """
        result = {}
        for entry in entries:
            lab = entry.label if not entry.label is None else Lit.stNone
            if result.has_key(lab):
                result[lab].append(entry)
            else:
                result[lab] = [entry]

        return result
        
    def addDateTimeToEntries(self,db,entries,fromgp=False):
        """
        Add Time to entries from series of gp
        """
        # connect to db
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[addDateTimeToEntries] Error: %s" % (e.args[0]))

        for entry in entries:
            # add datetime
            if fromgp: # the original list is from gewenste_processen, so first fill datetime with series
                cur.execute("SELECT * FROM %s WHERE pk=%d" % (Lit.table_series,entry.seriesid))
                row = cur.fetchone()
                if row is None:
                    continue
                entry.date = row['pps_start']
                if self.ignore_pps_start:
                    entry.date = None
                entry.studyid=row['study_fk']
                entry.stationname=row['station_name']
                entry.srcaet=row['src_aet']
                entry.modality=row['modality']
                entry.description=row['series_desc']

            if entry.date is None:
                cur.execute("SELECT content_datetime FROM %s WHERE series_fk=%d" % (Lit.table_instances,entry.seriesid))
                instance_rows = cur.fetchall()
                seriesdate = min ( [ row['content_datetime'] for row in instance_rows] )
                entry.date = seriesdate
            if entry.date is None:
                cur.execute("SELECT study_datetime FROM %s WHERE pk=%d" % (Lit.table_studies,entry.studyid))
                row = cur.fetchone()
                entry.date = row['study_datetime']
            if entry.description is None:
                cur.execute("SELECT study_desc FROM %s WHERE pk=%d" % (Lit.table_studies,entry.studyid))
                row = cur.fetchone()
                entry.description = row['study_desc']

        ### close db
        self.closedb(con)
        
        return entries
        
    def readSeries(self,db,note=None):
        """
        Make a list of entries of series in given db
        """
        # connect to db
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[readSeries] Error: %s" % (e.args[0]))

        results = []
        cur.execute("SELECT * FROM %s" %Lit.table_series)
        rows = cur.fetchall()
        results = [ Entry(seriesid=row['pk'], studyid=row['study_fk'], stationname=row['station_name'], 
                         srcaet=row['src_aet'], modality=row['modality'],description=row['series_desc'],
                         date=row['pps_start'], gpid=None,note=note) for row in rows ]
        ### close db
        self.closedb(con)

        # add missing info
        results = self.addDateTimeToEntries(db,results,fromgp=False)        
        
        return results

    def readGewenst(self,db,note=None,selector=None):
        """
        Make a list of entries of series in given db
        """
        # connect to db
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[readGewenst] Error: %s" % (e.args[0]))

        results = []
        if selector is None:
            cur.execute("SELECT * FROM %s" %Lit.table_gewenste_processen)
        else:
            cur.execute("SELECT * FROM %s WHERE selector_fk=%d" %(Lit.table_gewenste_processen,selector))

        rows = cur.fetchall()
        results = [ Entry(seriesid=row['series_fk'], studyid=row['study_fk'],date=row['creation_time'], gpid=row['pk'],status=row['status'],note=note) for row in rows ] 

        ### close db
        self.closedb(con)
        
        return results

    def readSelectors(self,db):
        """
        Make a dictionary of selectornames to ids
       """
        # connect to db
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[readGewenst] Error: %s" % (e.args[0]))

        results = {}
        cur.execute("SELECT * FROM %s" % 'selector')
        rows = cur.fetchall()
        for row in rows:
            results[row['name']] = row['pk']

        ### close db
        self.closedb(con)
        
        return results

    def readModalitySpecs(self,db):
        # read list of unique modalities in iqc.selector_series
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[readModalitySpecs] Error: %s" % (e.args[0]))

        modalities = []
        descriptions = {}
        aets = {}
        stationnames = {}
        cur.execute("SELECT * FROM %s" % 'selector_series')
        rows = cur.fetchall()
        for row in rows:
            mo = row['modality']
            if not mo in modalities:
                modalities.append(mo)
                descriptions[mo] = []
                aets[mo] = []
                stationnames[mo] = []
    
            desc = row['series_desc']
            if desc is None:
                desc = Lit.stNone
            if not desc in descriptions[mo]:
                descriptions[mo].append(desc)
    
            desc = row['station_name']
            if desc is None:
                desc = Lit.stNone
            if not desc in stationnames[mo]:
                stationnames[mo].append(desc)
        
            desc = row['src_aet']
            if desc is None:
                desc = Lit.stNone
            if not desc in aets[mo]:
                aets[mo].append(desc)
                
        ### close db
        self.closedb(con)
        
        return modalities, descriptions,aets,stationnames
    
    def qcStatusEntry(self,db,entries,qcIgnores):
        """
        For entries as a dictionary
        """
        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[readModalitySpecs] Error: %s" % (e.args[0]))

        for key,entry in entries.iteritems():
            # check status of analysis; only look for results if there are any
            if entry.status <5 :
                # report "Waiting for analysis...."
                entry.qcstatus = 'waiting for analysis...'
                continue

            if entry.status == 10 :
                # report "ERROR! Analysis failed!"
                entry.qcstatus = Lit.stResultERROR+' analysis failed'
                continue
    
            if entry.status == 20 :
                # report "ERROR! Analysis removed!"
                entry.qcstatus = Lit.stResultERROR+' analysis removed'
                continue
    
            # look if all results adhere to limits
            check_error = ''
            check_warning = ''

            # first check char limits
            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (Lit.table_char,entry.gpid))
            rows = cur.fetchall()
            for row in rows:
                # does it have a criterion?
                criterium = row['criterium']
                if criterium is None or len(criterium) == 0:
                    continue
                # should this value be ignored?
                omschrijving = row['omschrijving']
                ignore = False
                for ign in qcIgnores:
                    if ign in omschrijving:
                        ignore = True
                        break
                if ignore:
                    continue
    
                if not row['waarde'] == criterium:
                    if check_error != '':
                        check_error += ':'
                    check_error += omschrijving
    
            # check float limits
            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (Lit.table_float,entry.gpid))
            rows = cur.fetchall()
            for row in rows:
                kr_onder = row['grens_kritisch_onder']
                kr_boven = row['grens_kritisch_boven']
                if kr_boven is None or kr_onder is None:
                    continue
                ac_onder = row['grens_acceptabel_onder']
                ac_boven = row['grens_acceptabel_boven']
                if ac_boven is None or ac_onder is None:
                    continue

                ignore = False
                omschrijving = row['omschrijving']
                for ign in qcIgnores:
                    if ign in omschrijving:
                        ignore = True
                        break
                if ignore:
                    continue
    
                val = row['waarde']
                if val < kr_onder or val> kr_boven:
                    if check_error != '':
                        check_error += ':'
                    check_error += omschrijving
                elif val < ac_onder or val> ac_boven:
                    if check_warning != '':
                        check_warning += ':'
                    check_warning += omschrijving

            if check_error != '':
                entry.qcstatus = Lit.stResultERROR+' '+check_error
            if check_warning != '':
                if not entry.qcstatus is None:
                    entry.qcstatus += ': '
                else:
                    entry.qcstatus = ''
                entry.qcstatus += Lit.stResultWARNING+' '+check_warning
            if entry.qcstatus is None:
                entry.qcstatus = Lit.stResultOK
        
        ### close db
        self.closedb(con)
        
        return entries
        
    def _latestTestSelectors(self,qcstatus,gp_lab=None,db=Lit.db_iqc):
        """
        Returns a dictionary of the entries with for each selector the series with the latest acquisition date
        """

        ###1 open db iqc
        try:
            con,cur,dbversion = self.connectdb(host=self.host,db=db)
        except ValueError,e:
            raise ValueError("[latestTest] Error: %s" % (e.args[0]))

        sel2id = self.readSelectors(db)
        
        results = {}
        inperiod = 0
        for sel in qcstatus.selectors:
            # find the id of the selector
            try:
                selid = sel2id[sel]
            except ValueError,e:
                raise ValueError("[latestTest] Error: %s" % (e.args[0]))
            
            # read all sels from gewenste processen. 
            gps = self.readGewenst(db,selector=selid)

            # For each gewenste processen add acc time
            gps = self.addDateTimeToEntries(db, gps,fromgp=True)

            # select only relevant (latest 6 months)
            gps = [ g for g in gps if g.date >= qcstatus.period[1] ]
            if not gp_lab is None:
                gps = self.addLabelToSeries(db, gps,gp_lab)

            gpss = self.splitOnLabels(gps)

            # select entry with latest datetime
            for key,gps in gpss.iteritems():
                if len(gps)>0:
                    results[sel+key] = max(gps,key=attrgetter('date'))
                    inperiod += len(gps) # sum(1 for x in gps if x.date >= qcstatus.period[1])
                    
        if len(results.keys()) >1:
            inperiod /= len(results.keys())
        qcstatus.period = (qcstatus.period[0],qcstatus.period[1],inperiod)
        results = self.qcStatusEntry(db, results, qcstatus.qcIgnores)

        ### close db
        self.closedb(con)
        
        return results
        
    def last_day_of_month(self,any_day):
        # calculates the last day of month
        next_month = any_day.replace(day=28) + datetime.timedelta(days=4)  # this will never fail
        return (next_month - datetime.timedelta(days=next_month.day)).day

    def _production(self,qcstatus):
        """
        Calculates percentage of tests done over expected
        """
    
        # 1. Find same date 6 months ago
        dt_now = datetime.date.today()
        y1 = dt_now.year
        m1 = dt_now.month-6 # semi-anual
        if m1 < 1:
            y1 -= 1
            m1 += 12
        max_day_of_month = self.last_day_of_month( datetime.date(year=y1,month=m1,day=1) )
        dt_start = datetime.date(year=y1,month=m1,day=min(dt_now.day,max_day_of_month))

        # 2. Find number of days/weeks/months is reporting period
        delta_days = (dt_now - dt_start).days
        delta_weeks = delta_days/7
        delta_months = dt_now.month - dt_start.month+ 12*(dt_now.year - dt_start.year)

        # 3. Determine expected number of processes
        freq = qcstatus.frequency[0]
        if 'monthly' in freq:
            num_expected = delta_months
        elif 'bi-weekly' in freq:
            num_expected = delta_weeks/2
        elif 'weekly' in freq:
            num_expected = delta_weeks
        elif 'daily' in freq:
            num_expected = delta_days
        
        qcstatus.period = (qcstatus.period[0],qcstatus.period[1],qcstatus.period[2],num_expected) # append num_expected
        comment = '%s %d %d/%d'%(qcstatus.machine,len(qcstatus.selectors),qcstatus.period[2],num_expected)
        qcstatus.production = 100.*qcstatus.period[2]/num_expected
        return qcstatus,comment

    def latestTestMachine(self,machine,ininame,gp_lab=None,db=Lit.db_iqc):
        #print '[latestTestMachine]',machine
        if gp_lab is None:
            gp_lab = self._labelsFromGP(db)
        
        status = QCStatus(machine,ininame)
        results = self._latestTestSelectors(status,gp_lab,db)
        if len(results)>0:
            status.date = max(results.itervalues(),key=attrgetter('date')).date
            status.overdue = status.date<status.frequency[1]

        status.status_full = ''
        status.status_short = ''
        errors = []
        warnings = []
        for key,entry in results.iteritems():
            if entry.qcstatus == Lit.stResultOK:
                continue
            if Lit.stResultERROR in entry.qcstatus:
                errors.extend(entry.qcstatus.split(Lit.stResultERROR)[1].split(Lit.stResultWARNING)[0].strip(' :').split(':'))
            if Lit.stResultWARNING in entry.qcstatus:
                warnings.extend(entry.qcstatus.split(Lit.stResultWARNING)[1].strip(' :').split(':'))
        
        if len(errors)>0:
            status.status_short = "%d ERRORS"%len(errors)
            status.status_full = Lit.stResultERROR+' '+':'.join(errors)
        if len(warnings)>0:
            if len(status.status_short)>0:
                status.status_short += ', '
                status.status_full += ', '
            status.status_short += "%d WARNINGS"%len(warnings)
            status.status_full += Lit.stResultWARNING+' '+':'.join(warnings)

        #print status.machine, status.date, status.status_full
        #print status.machine, status.date, status.status_short,status.overdue,status.period[2]

        status,comment = self._production(status)
        return status,comment
