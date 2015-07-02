#!/usr/bin/python
"""
Changelog:
20150617: Fix crash on empty value; SeekAndDestroy (SrcAET) added; restructured SeekAndDestroy for SrcAET,Modality,StationName
20150413: SeekAndDestroy (modality)
20150407: SeekAndDestroy (selector); fix empty stationname
20150317: Fix crash on empty value; added SrcAET info
20150305: Fix for empty selectors crash (EPIStab)
20150210: Fix for crash on trying to add a None to qctable
20141010: Fix remove all traces of one StationName; added remove all traces of one GewensteID; addAnatomy voor CT2
20141008: Add option to remove all traces of one StationName
20140627: added MIR excerpt
20140618: scaling between -1 and 1
20140617: Added Current Status Report
20140616: Allow "None" for station_name or modality; stacked graphs; critical lines
20140613: distinguishable colors in report; sorting on date in report
20140611: reportgenerator added
20140610: ability to select items to show/hide in report
20140606: added FilterName to description of MG studies
20140603: added reporting
20140528: added 'purge' option to Database to delete XML files
20140521: initial version

TODO:
    o Herschrijf onSeekandDestroyStation net als onSeekandDestroyGewenst
    o Bouw optie om data te corrigeren (wrong or no stationname, etc)
"""
import argparse
from collections import OrderedDict
from pyqtgraph.Qt import QtGui,QtCore
import pyqtgraph as pg
import MySQLdb as mdb
import sys
from ConfigParser import SafeConfigParser
import os
import datetime
import time
import functools
from reporter import Reporter

def datetime2seconds(dt):
    return time.mktime(dt.timetuple())

class reportentry:
    gpid = None
    seriesid = None
    stationname = None
    status = None
    selector = None
    modality = None
    date = None
    description = None # concat of series_description and FilterName if exists in resultaten_char
    analysisresult = None
    analysisdate = None
    srcaet = None
    def __init__(self,gpid=None,stationname=None,modality=None,date=None,selector=None,status=None, description=None,seriesid=None,analysisresult=None,analysisdate=None,srcaet=None):
        self.gpid = gpid
        self.stationname = stationname
        self.modality = modality
        self.date = date
        self.selector = selector
        self.status = status
        self.description = description
        self.seriesid = seriesid
        self.analysisresult = analysisresult
        self.analysisdate = analysisdate
        self.srcaet = srcaet

class dbTool(QtGui.QMainWindow):
    qcversion = 20150617
    verbose = False

    stModeProcessor = 'processor'
    stModeReport = 'report'
    stModeDestroy = 'destroy'

    runmode = None
    host = None
    user = None
    pswd = None
    rootxml = None

    iqcdb = 'iqc'
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

    status_reset = 0

    qctable = None
    qcgraph = None
    qcgroupbox = None
    qcplots = None # dictionary of measurement_names (plot)
    qctext = None
    statusLabel = None
    selectedModality = None
    selectedStation = None
    selectedSelector = None
    selectedStatus = None
    selectedGewenst = None

    selectedReportPeriod = None
    selectedReportQCPeriodicity = None
    selectedReportDateTimestart = None
    selectedReportDateTimeuntil = None
    reportentries = {}
    reporter = None

    status_id_text = {}
    selector_id_text = {}

    truncate_tables = [
        'analysemodule_input',
        'analysemodule_output',
        'collector_series_status',
        'collector_study_status',
        'files',
        'gewenste_processen',
        'instance',
        'patient',
        'resultaten_boolean',
        'resultaten_char',
        'resultaten_floating',
        'resultaten_object',
        'resultaten_status',
        'series',
        'study',
        ]

    stNone = 'NONE'
    stResultOK = 'Ok'
    stResultOutOfRange = 'Out of Range'
    stResultNoRanges = 'No limits'
    stResultError = 'Error'
    destroylist = []

    def pswdFromDialog(self):
        diag = QtGui.QInputDialog()
        title = "Enter password"
        text = "Password for mysql user %s at host %s: " % (self.user,self.host)
        s = None
        while s is None:
            qstring, ok = diag.getText(self, QtCore.QString(title), QtCore.QString(text), mode=QtGui.QLineEdit.Password)
            s = str(qstring)
            if ok is False: # user pressed Cancel
                return None
            if s == '':     # user entered nothing
                s = None
        return s

    def configError(self,ininame):
        print "ERROR: configfile %s corrupt or non-existent. Should look like:\n[iqc]\nhost = localhost\n" \
              "user = root\npswd = ***\n\nwith valid values for host and user.\nOptionally xmlroot=/opt can be supplied\n" \
              "if pswd is blank then a dialog will ask for a valid password\n" \
              "alternatively pswd2 = xxx can be given, with xxx the base64encoded password"% ininame
        QtGui.qApp.quit()
        exit()

    def configRead(self,ininame):
        parser = SafeConfigParser()
        parser.read(ininame)
        section = 'iqc'
        if not parser.has_section(section):
            self.configError(ininame)
        else:
            candidates = ['host','user','pswd']
            for can in candidates:
                if not parser.has_option(section,can):
                    self.configError(ininame)
        self.host = parser.get(section,'host')
        self.user = parser.get(section,'user')
        self.pswd = parser.get(section,'pswd')
        if parser.has_option(section,'xmlroot'):
            self.xmlroot = parser.get(section,'xmlroot')
        else:
            self.xmlroot = None

        if self.pswd is '':
            import base64
            pswd2 = parser.get(section,'pswd2')
            if pswd2 is '':
                self.pswd = self.pswdFromDialog()
            else:
                try:
                    self.pswd = base64.b64decode(pswd2)
                except:
                    self.pswd = self.pswdFromDialog()
                    print 'add line pswd2 = %s to dbTool.ini' % base64.b64encode(self.pswd)

    def __init__(self, parent=None):
        super(dbTool, self).__init__(parent)
        """Config"""
        self.configRead('dbTool.ini')

        """Constructor"""
        self.qw = QtGui.QWidget()
        self.setCentralWidget(self.qw)
        self.layout = QtGui.QGridLayout()

        self.statusBar()
        self.statusBar().showMessage('Ready')

        self.setGeometry(1, 1, 700, 512)
        self.setWindowTitle("dbTool v"+str(self.qcversion));
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        exitAction = QtGui.QAction("&Quit", self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(QtGui.qApp.quit)
        fileMenu.addAction(exitAction)

        processorMenu = menubar.addMenu("&Processor")
        procresetAction = QtGui.QAction("Reset: Select which processor status to reset", self)
        procresetAction.triggered.connect(self.processorReset)
        processorMenu.addAction(procresetAction)

        databaseactions = [
            ("Truncate: remove all results and datasets",self.databaseTruncatePopUp),
            ("Purge: Truncate database and remove XML files",self.popupPurge),
            ("Seek and Destroy: remove all iqc database traces of a given Src-AET (will start reimport from dcm4chee)",functools.partial(self.seekandDestroy,mode='Src-AET')),
            ("Seek and Destroy: remove all iqc database traces of a given StationName (will start reimport from dcm4chee)",functools.partial(self.seekandDestroy,mode='StationName')),
            ("Seek and Destroy: remove all iqc database traces of a given GewenstID (will start reimport from dcm4chee)",functools.partial(self.seekandDestroy,mode='GewenstID')),
            ("Seek and Destroy: remove all iqc database traces of a given Selector (will start reimport from dcm4chee)",functools.partial(self.seekandDestroy,mode='Selector')),
            ("Seek and Destroy: remove all iqc database traces of all selectors of a given Modality (will start reimport from dcm4chee)",functools.partial(self.seekandDestroy,mode='Modality')),
        ]
        databaseMenu = menubar.addMenu("&Database")
        for ra in databaseactions:
            dbAction = QtGui.QAction(ra[0], self)
            dbAction.triggered.connect(ra[1])
            databaseMenu.addAction(dbAction)

        reportactions = [
            ("Periodic status report",self.periodicStatus),
            ("Current status report per selector+description",functools.partial(self.currentStatus,mode='selector+description')),
            ("Current status report per selector",functools.partial(self.currentStatus,mode='selector')),
            ("Demo status report",self.reportDemo),
            ("MIR Bucky report",self.reportMIRBucky),
            ]
        reportMenu = menubar.addMenu("&Reports")
        for ra in reportactions:
            reportAction = QtGui.QAction(ra[0], self)
            reportAction.triggered.connect(ra[1])
            reportMenu.addAction(reportAction)

        # create progress bar
        self.statusLabel = QtGui.QLabel()
        self.statusBar().addPermanentWidget(self.statusLabel)

    def clearLayout(self, layout,layoutOnly=False):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clearLayout(item.layout())

        if not layoutOnly:
            self.qcgraph = None
            self.qctable = None

    def connectdb(self,host=None,user=None,pswd=None,db=None):
        if host is None:
            host = self.host
        if user is None:
            user = self.user
        if pswd is None:
            pswd = self.pswd
        if db is None:
            db = self.iqcdb

        try:
            con = mdb.connect(host,user,pswd,db) # host, username, pswd, db
            cur = con.cursor(mdb.cursors.DictCursor) # dictionary format, i.e. access by column name
            cur.execute("select * from config")
            rows = cur.fetchall()
            dbversion='0'
            for row in rows:
                if row['property'] == "Version_Database":
                    dbversion = str(row['value'])
        except mdb.Error, e:
            print "[connectdb] Error %d: %s" % (e.args[0],e.args[1])
            sys.exit(1)
        return con,cur,dbversion

    def closedb(self,con):
        if con:
            con.close()

    def processorShowSelected(self,selected):
        columnsToShow=[
            'status',
            'selector_fk',
            'creation_time'
            ]
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        selrows = []
        headers = []
        for pkid in selected:
            cur.execute("select * from %s where pk=%d" % (self.table_gewenste_processen,pkid))
            rows = cur.fetchall()
            selrows.append(rows[0])  # should be only 1
            if len(headers) == 0:
                for col in rows[0].keys():
                    if col in columnsToShow:
                        headers.append(col)

        headers.append('Modality')
        headers.append('StationName')
        headers.append('SrcAET')
        headers.append('Description')
        headers.append('GewenstID')
        if 1> 0 or self.runmode == self.stModeReport:
            headers.append('DateTime')

        self.qctable.setRowCount(len(selrows))
        if len(selected) > 0:
            self.qctable.setColumnCount(len(headers))
            rowid = 0
            self.qctable.setHorizontalHeaderLabels(headers)
            for row in selrows:
                colid = 0
                gpid = row['pk']
                for hd in headers:
                    if self.status_id_text and hd == 'status':
                        valstr = self.status_id_text[row[hd]]
                        self.reportentries[gpid].status = valstr
                    elif hd == 'selector_fk':
                        valstr = self.selector_id_text[row[hd]]
                        if gpid in self.reportentries:
                            self.reportentries[gpid].selector = valstr
                    elif hd == 'Modality':
                        valstr = self.reportentries[gpid].modality
                    elif hd == 'StationName':
                        valstr = self.reportentries[gpid].stationname
                    elif hd == 'SrcAET':
                        valstr = self.reportentries[gpid].srcaet
                    elif hd == 'Description':
                        valstr = self.reportentries[gpid].description
                    elif hd == 'DateTime':
                        valstr = self.reportentries[gpid].date.strftime('%Y/%m/%d %H:%M:%S')
                    elif hd == 'GewenstID':
                        valstr = str(self.reportentries[gpid].gpid)
                    else:
                        valstr = str(row[hd])
                    if not valstr is None:
                        self.qctable.setItem(rowid, colid, QtGui.QTableWidgetItem(valstr))
                    colid += 1
                rowid += 1
        ### close db
        self.closedb(con)

        self.statusLabel.setText("Selected %d processes." %len(selrows))

    def onActivatedProcessorSelector(self,text):
        if self.verbose:
            print "[onActivatedProcessorSelector]",text,self.selector_id_text[str(text)]
        self.selectedSelector = self.selector_id_text[str(text)]

        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeDestroy:
            self.onSeekAndDestroySelector()
            return
        if self.runmode == self.stModeReport:
            self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        if self.runmode == self.stModeReport:
            self.onActivatedReportChanged()

    def onActivatedProcessorStatus(self,text):
        if self.verbose:
            print "[onActivatedProcessorStatus]",str(text),self.status_id_text[str(text)]
        self.selectedStatus = self.status_id_text[str(text)]

        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeReport:
            self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        if self.runmode == self.stModeReport:
            self.onActivatedReportChanged()

    def onActivatedGewenstProces(self,text):
        if self.verbose:
            print "[onActivatedGewenstProces]",text
        station = str(text)
        self.selectedGewenst = station

        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeDestroy:
            self.onSeekAndDestroyGewenst()
            return

    def onActivatedProcessorStation(self,text):
        if self.verbose:
            print "[onActivatedProcessorStation]",text
        station = str(text)
        self.selectedStation = station

        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeDestroy:
            self.onSeekAndDestroyStation()
            return

        if self.runmode == self.stModeReport:
            self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        if self.runmode == self.stModeReport:
            self.onActivatedReportChanged()

    def onActivatedProcessorSrcAET(self,text):
        if self.verbose:
            print "[onActivatedProcessorSrcAET]",text
        station = str(text)
        self.selectedAET = station

        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeDestroy:
            self.onSeekAndDestroyStation(mode='SrcAET')
            return

        if self.runmode == self.stModeReport:
            self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        if self.runmode == self.stModeReport:
            self.onActivatedReportChanged()

    def onActivatedProcessorModality(self,text):
        if self.verbose:
            print "[onActivatedProcessorModality]",text
        mod = str(text)
        self.selectedModality = mod
        if self.qctable is None or self.statusLabel is None:
            return
        if self.runmode == self.stModeDestroy:
            self.onSeekAndDestroyStation(mode='Modality')
            return

        if self.runmode == self.stModeReport:
            self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        if self.runmode == self.stModeReport:
            self.onActivatedReportChanged()

    def uniqifyList(self,multilist):
        return sorted(OrderedDict.fromkeys(multilist).keys())

    def onSeekAndDestroySelector(self):
        self.onSeekAndDestroyGewenst(self,selector=self.selectedSelector)

    def onSeekAndDestroyGewenst(self,act=None,selector=None):
        self.destroylist = [(None,None)] # fill this one up at the end with XMLfiles
        if self.qctext:
            self.qctext.clear()

        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        """ Plan: Remove all traces from iqc, afterward remove xml files; data in dcm4che should stay intact, so rebuild should happen
        affected tables:
            01. analysemodule_input: /opt/XML/analysemodule_input/*.xml files for removal
            01. analysemodule_output: /opt/XML/analysemodule_input/<Selector Name> for removal
            2. collector_series_status: gekoppeld met series via series_fk
            3. collector_study_status: gekoppeld met study via study_fk
            4. files: gekoppeld met instance via instance_fk
            01. gewenste_processen: gekoppeld met selector/study/series/instance/analysemodule_input/analysemodule_output via xxx_fk
            01. instance: gekoppeld met series via series_fk
            3. patient
            01. series: gekoppeld met study via study_fk, heeft station_name
            01. study: gekoppeld met patient via patient_fk
            6. resultaten_boolean: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_char: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_floating: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_object: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_status: gekoppeld aan gewenste_processen via gewenste_processen_fk
        needs for xml:
            Selector Name

        FlowChart:
            01. Maak lijsten Processen, Series, Studies, Instances, AnalyseInput, AnalyseOutput uit 'gewenste_processen'
            02. Als instance, controlleer dan of alle instances uit serie weg, zo ja dan ook serie weg
            02. Als series, controlleer dan of alle series uit study weg, zo ja dan ook study weg
            02. Als study, controlleer dan of alle studies uit patient weg, zo ja dan ook patient weg
            03. Van boven naar beneden lijsten completeren: Pat->stud->series->instances->files
            04. Dan toevoegen van alles dat indirect uit Studies,Series,Processen volgt
            05. list of xml files /dirs to remove
        """
        # 01. Maak lijsten Processen, Series, Studies, Instances, AnalyseInput, AnalyseOutput uit 'gewenste_processen'
        daProcessen = []
        daInput     = []
        daOutput    = []
        daInstances = []
        daSeries    = []
        daStudies   = []
        if selector is None:
            cur.execute("SELECT * FROM %s where pk=%s" % (self.table_gewenste_processen,self.selectedGewenst))
        else:
            cur.execute("SELECT * FROM %s where selector_fk=%s" % (self.table_gewenste_processen,selector))
            
        rows_processen = cur.fetchall()
        checklist = [
            ('pk',                     daProcessen,self.table_gewenste_processen),
            ('study_fk',               daStudies,  self.table_studies),
            ('series_fk',              daSeries,   self.table_series),
            ('instance_fk',            daInstances,self.table_instances),
            ('analysemodule_input_fk', daInput,    self.table_analyseinput),
            ('analysemodule_output_fk',daOutput,   self.table_analyseoutput),
        ]
        for row in rows_processen:
            for key,lis,tab in checklist:
                if row[key]:
                    lis.append(row[key])
        daProcessen = self.uniqifyList(daProcessen)
        daStudies   = self.uniqifyList(daStudies)
        daSeries    = self.uniqifyList(daSeries)
        daInstances = self.uniqifyList(daInstances)
        daInput     = self.uniqifyList(daInput)
        daOutput    = self.uniqifyList(daOutput)
        for key,lis,tab in checklist:
            if selector is None:
                feedback  = "1. Found %d entries in %s for %s" %(len(lis),tab,self.selectedGewenst)
            else:
                feedback  = "1. Found %d entries in %s for %s" %(len(lis),tab,selector)
            if self.qctext:
                self.qctext.appendPlainText(feedback)
            else:
                print feedback

        # 2. Als instance, controlleer dan of alle instances uit serie weg, zo ja dan ook serie weg
        # 2. Als series, controlleer dan of alle series uit study weg, zo ja dan ook study weg
        # 2. Als study, controlleer dan of alle studies uit patient weg, zo ja dan ook patient weg
        daPatients = []
        checklist = [
            #(tabel_laagnivo, lijst_laagnivo, key_hoognivo, lijst_hoognivo
            (self.table_instances, daInstances, 'series_fk',  daSeries),
            (self.table_series,    daSeries,    'study_fk',   daStudies),
            (self.table_studies,   daStudies,   'patient_fk', daPatients),
        ]
        for tab_lo,list_lo, key_hi,list_hi in checklist:
            high_lows = {}

            cur.execute("SELECT * FROM %s" %tab_lo)
            rows_processen = cur.fetchall()
            for row in rows_processen:
                if row['pk'] in list_lo:
                    if row[key_hi]:
                        high_lows[row[key_hi]]=[ ]
            for row in rows_processen:
                if row[key_hi]:
                    if row[key_hi] in high_lows:
                        high_lows[row[key_hi]].append(row['pk'])

            for hi,los in high_lows.iteritems():
                numdel = 0
                for lo in los:
                    if lo in list_lo:
                        numdel += 1
                    else:
                        print "%s sub %d of %d in %s is not scheduled for delete!" % (key_hi,lo,hi,tab_lo)
                if numdel == len(los):
                    list_hi.append(hi)
        daPatients = self.uniqifyList(daPatients)
        daStudies  = self.uniqifyList(daStudies)
        daSeries   = self.uniqifyList(daSeries)
        feedback  = "2. Selected %d Patients" % len(daPatients)
        feedback  += "\n2. Selected %d Studies" % len(daStudies)
        feedback  += "\n2. Selected %d Series" % len(daSeries)
        feedback  += "\n2. Selected %d Instances" % len(daInstances)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 3. Van boven naar beneden lijsten completeren: Pat->stud->series->instances->files
        # 4. Dan toevoegen van alles dat indirect uit Studies,Series,Processen volgt
        daFiles = []
        daCollectorStudyStatus  = []
        daCollectorSeriesStatus = []
        daBools  = []
        daChars  = []
        daFloats  = []
        daObjects = []
        daStats   = []
        checklist = [
            #(tabel_laagnivo, lijst_laagnivo, key_hoognivo, lijst_hoognivo
            (self.table_studies,   daStudies,   'patient_fk',  daPatients),
            (self.table_series,    daSeries,    'study_fk',    daStudies),
            (self.table_instances, daInstances, 'series_fk',   daSeries),
            (self.table_files,     daFiles,     'instance_fk', daInstances),
            (self.table_collectorstatusstudies, daCollectorStudyStatus,    'study_fk',  daStudies),
            (self.table_collectorstatusseries,  daCollectorSeriesStatus,   'series_fk', daSeries),
            (self.table_bool,   daBools, 'gewenste_processen_fk', daProcessen),
            (self.table_char,   daChars, 'gewenste_processen_fk', daProcessen),
            (self.table_float,  daFloats, 'gewenste_processen_fk', daProcessen),
            (self.table_object, daObjects, 'gewenste_processen_fk', daProcessen),
            (self.table_status, daStats, 'gewenste_processen_fk', daProcessen),
        ]
        for tab_lo, list_lo, key_hi,list_hi in checklist:
            cur.execute("SELECT * FROM %s" %tab_lo)
            rows_processen = cur.fetchall()
            for row in rows_processen:
                if row[key_hi]:
                    if row[key_hi] in list_hi:
                        list_lo.append(row['pk'])
        daPatients  = self.uniqifyList(daPatients)
        daStudies   = self.uniqifyList(daStudies)
        daSeries    = self.uniqifyList(daSeries)
        daInstances = self.uniqifyList(daInstances)
        daFiles     = self.uniqifyList(daFiles)
        daCollectorStudyStatus  = self.uniqifyList(daCollectorStudyStatus)
        daCollectorSeriesStatus = self.uniqifyList(daCollectorSeriesStatus)
        daBools   = self.uniqifyList(daBools)
        daChars   = self.uniqifyList(daChars)
        daFloats  = self.uniqifyList(daFloats)
        daObjects = self.uniqifyList(daObjects)
        daStats = self.uniqifyList(daStats)
        feedback  = "3. Corrected to %d Patients" % len(daPatients)
        feedback  += "\n3. Corrected to %d Studies" % len(daStudies)
        feedback  += "\n3. Corrected to %d Series" % len(daSeries)
        feedback  += "\n3. Corrected to %d Instances" % len(daInstances)
        feedback  += "\n3. Selected %d Files" % len(daFiles)
        feedback  += "\n3. Selected %d CollectorStudyStatus" % len(daCollectorStudyStatus)
        feedback  += "\n3. Selected %d CollectorSeriesStatus" % len(daCollectorSeriesStatus)
        feedback  += "\n3. Selected %d Bools" % len(daBools)
        feedback  += "\n3. Selected %d Chars" % len(daChars)
        feedback  += "\n3. Selected %d Floats" % len(daFloats)
        feedback  += "\n3. Selected %d Objects" % len(daObjects)
        feedback  += "\n3. Selected %d Status" % len(daStats)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 5. list of xml files /dirs to remove
        daXML = []
        cur.execute("SELECT * FROM %s" % self.table_analyseinput)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['pk'] in daInput:
                daXML.append("%s%s"%(row['filepath'],row['filename']))
        feedback = "4. Found %d/%d %s for %s" %(len(daXML),len(rows_processen),"inputXMLs",self.selectedStation)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        len0 = len(daXML)
        cur.execute("SELECT * FROM %s" % self.table_analyseoutput)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['pk'] in daOutput:
                daXML.append("%s%s"%(row['filepath'],os.path.dirname(row['filename'])))
        feedback = "4. Found %d/%d %s for %s" %(len(daXML)-len0,len(rows_processen),"outputXMLs",self.selectedStation)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # Make destroylist
        self.destroylist = []
        self.destroylist.append((None,daXML))
        self.destroylist.append((self.table_patient,daPatients))
        self.destroylist.append((self.table_studies,daStudies))
        self.destroylist.append((self.table_series,daSeries))
        self.destroylist.append((self.table_instances,daInstances))
        self.destroylist.append((self.table_files,daFiles))
        self.destroylist.append((self.table_collectorstatusstudies,daCollectorStudyStatus))
        self.destroylist.append((self.table_collectorstatusseries,daCollectorSeriesStatus))
        self.destroylist.append((self.table_bool,daBools))
        self.destroylist.append((self.table_char,daChars))
        self.destroylist.append((self.table_float,daFloats))
        self.destroylist.append((self.table_object,daObjects))
        self.destroylist.append((self.table_status,daStats))

        self.destroylist.append((self.table_gewenste_processen,daProcessen))
        self.destroylist.append((self.table_analyseinput,daInput))
        self.destroylist.append((self.table_analyseoutput,daOutput))

        ### close db
        self.closedb(con)
        QtGui.QApplication.restoreOverrideCursor()

        return


    def onSeekAndDestroyStation(self,act=None,mode='StationName'):
        self.destroylist = [(None,None)] # fill this one up at the end with XMLfiles
        if self.qctext:
            self.qctext.clear()

        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        """ Plan: Remove all traces from iqc, afterward remove xml files; data in dcm4che should stay intact, so rebuild should happen
        affected tables:
            5. analysemodule_input: /opt/XML/analysemodule_input/*.xml files for removal
            5. analysemodule_output: /opt/XML/analysemodule_input/<Selector Name> for removal
            2. collector_series_status: gekoppeld met series via series_fk
            3. collector_study_status: gekoppeld met study via study_fk
            4. files: gekoppeld met instance via instance_fk
            5. gewenste_processen: gekoppeld met selector/study/series/instance/analysemodule_input/analysemodule_output via xxx_fk
            2. instance: gekoppeld met series via series_fk
            3. patient
            1. series: gekoppeld met study via study_fk, heeft station_name
            1. study: gekoppeld met patient via patient_fk
            6. resultaten_boolean: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_char: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_floating: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_object: gekoppeld aan gewenste_processen via gewenste_processen_fk
            6. resultaten_status: gekoppeld aan gewenste_processen via gewenste_processen_fk
        needs for xml:
            Selector Name

        FlowChart:
            1. Maak lijst 'Series' en 'Studies' uit series met juiste station_name
            2. Gebruik 'Series': voor lijst 'Instances' uit instance en 'Collector_series_status' uit collector_series_status
            3. Gebruik 'Studies': voor lijst 'Patients' uit patient en 'Collector_study_status' uit collector_study_status
            4. Gebruik 'Instances': voor lijst 'Files' uit files
            5. Gebruik 'Series' en 'Studies' en 'Instances': voor lijsten 'Gewenste_processen' en 'Analysemodule_input','Analysemodule_output' uitgewenste_processen
            6. Gebruik 'Gewenste_processen': voor lijsten 'Resultaten_boolean','Resultaten_char','Resultaten_floating','Resultaten_object' uit
            7.
        """
        # 1: Maak lijst 'Series' en 'Studies' uit series met juiste station_name
        if mode == 'StationName': # station_name Mode
            if self.selectedStation == self.stNone:
                cur.execute("SELECT * FROM %s WHERE station_name IS NULL" %self.table_series)
            else:
                cur.execute("SELECT * FROM %s WHERE station_name='%s'" %(self.table_series,self.selectedStation))
            rows_processen = cur.fetchall()
            thisthing = self.selectedStation
        elif mode == 'Modality': #modality mode
            if self.selectedModality == self.stNone:
                cur.execute("SELECT * FROM %s WHERE modality IS NULL" %self.table_series)
            else:
                cur.execute("SELECT * FROM %s WHERE modality='%s'" %(self.table_series,self.selectedModality))
            rows_processen = cur.fetchall()
            thisthing = self.selectedModality
        elif mode == 'SrcAET': #modality mode
            if self.selectedAET == self.stNone:
                cur.execute("SELECT * FROM %s WHERE src_aet IS NULL" %self.table_series)
            else:
                cur.execute("SELECT * FROM %s WHERE src_aet='%s'" %(self.table_series,self.selectedAET))
            rows_processen = cur.fetchall()
            thisthing = self.selectedAET
        else:
            feedback = "Unknown mode %s"%mode
                
            if self.qctext:
                self.qctext.appendPlainText(feedback)
            else:
                print feedback
            return
        
        thisthing = '%s=%s'%(mode,thisthing)    
        feedback = "Found %d rows in %s for %s" %(len(rows_processen),self.table_series,thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        daSeries = []
        daStudies = []
        for row in rows_processen:
            daSeries.append(row['pk'])
            if row['study_fk']:
                daStudies.append(row['study_fk'])
        daSeries = self.uniqifyList(daSeries)
        self.destroylist.append( (self.table_series,daSeries) )
        feedback = "1. Found %d series for %s" %(len(daSeries),self.selectedStation)
        feedback += "\n1. Found %d/%d studies for %s" %(len(daStudies),len(rows_processen),thisthing)

        #selecteer voor delete alleen studies, die geen series meer hebben na delete series
        studseries = {}
        cur.execute("SELECT * FROM %s" %self.table_series)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['study_fk']:
                if row['study_fk'] in daStudies:
                    if row['study_fk'] in studseries:
                        studseries[row['study_fk']].append(row['pk'])
                    else:
                        studseries[row['study_fk']]=[ row['pk'] ]

        daStudies =[]
        for stud,sers in studseries.iteritems():
            numdel = 0
            for ser in sers:
                if ser in daSeries:
                    numdel += 1
                else:
                    print "Series %d of study %d is not scheduled for delete!" % (ser,stud)
            if numdel == len(sers):
                daStudies.append(stud)

        daStudies = self.uniqifyList(daStudies)
        self.destroylist.append( (self.table_studies,daStudies) )
        feedback += "\n1. Corrected to %d studies only for %s" %(len(daStudies),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 2. Gebruik 'Series': voor lijst 'Instances' uit instance en 'Collector_series_status' uit collector_series_status
        daInstances = []
        cur.execute("SELECT * FROM %s" % self.table_instances)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['series_fk']:
                if row['series_fk'] in daSeries:
                    daInstances.append(row['pk'])
        daInstances = self.uniqifyList(daInstances)
        self.destroylist.append( (self.table_instances,daInstances) )
        feedback= "2. Found %d/%d instances for %s" %(len(daInstances),len(rows_processen),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        daCollectorSeriesStatus = []
        cur.execute("SELECT * FROM %s" %self.table_collectorstatusseries)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['series_fk']:
                if row['series_fk'] in daSeries:
                    daCollectorSeriesStatus.append(row['pk'])
        daCollectorSeriesStatus = self.uniqifyList(daCollectorSeriesStatus)
        self.destroylist.append( (self.table_collectorstatusseries,daCollectorSeriesStatus) )
        feedback = "2. Found %d/%d collectorstatusseries for %s" %(len(daCollectorSeriesStatus),len(rows_processen),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 3. Gebruik 'Studies': voor lijst 'Patients' uit patient en 'Collector_study_status' uit collector_study_status
        daPatients = []
        cur.execute("SELECT * FROM %s" %self.table_studies)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['patient_fk']:
                if row['pk'] in daStudies:
                    daPatients.append(row['patient_fk'])
        daPatients = self.uniqifyList(daPatients)
        feedback = "3. Found %d/%d patients for %s" %(len(daPatients),len(rows_processen),thisthing)

        #selecteer voor delete alleen patients, die geen studies meer hebben na delete studies
        patstuds = {}
        for row in rows_processen:
            if row['patient_fk']:
                if row['patient_fk'] in daPatients:
                    if row['patient_fk'] in patstuds:
                        patstuds[row['patient_fk']].append(row['pk'])
                    else:
                        patstuds[row['patient_fk']]=[ row['pk'] ]

        daPatients =[]
        for pat,studs in patstuds.iteritems():
            numdel = 0
            for stud in studs:
                if stud in daStudies:
                    numdel += 1
            if numdel == len(studs):
                daPatients.append(pat)

        feedback += "\n3. Corrected to %d patients for only %s" %(len(daPatients),thisthing)
        self.destroylist.append( (self.table_patient,daPatients) )
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        daCollectorStudyStatus = []
        cur.execute("SELECT * FROM %s" % self.table_collectorstatusstudies)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['study_fk']:
                if row['study_fk'] in daStudies:
                    daCollectorStudyStatus.append(row['pk'])
        daCollectorStudyStatus = self.uniqifyList(daCollectorStudyStatus)
        self.destroylist.append( (self.table_collectorstatusstudies,daCollectorStudyStatus) )
        feedback = "3. Found %d/%d collectorstatusstudies for %s" %(len(daCollectorStudyStatus),len(rows_processen),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 4. Gebruik 'Instances': voor lijst 'Files' uit files
        daFiles = []
        cur.execute("SELECT * FROM %s" %self.table_files)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['instance_fk']:
                if row['instance_fk'] in daInstances:
                    daFiles.append(row['pk'])
        daFiles = self.uniqifyList(daFiles)
        self.destroylist.append( (self.table_files,daFiles) )
        feedback = "4. Found %d/%d files for %s" %(len(daFiles),len(rows_processen),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 5. Gebruik 'Series' en 'Studies' en 'Instances': voor lijsten 'Gewenste_processen' en 'Analysemodule_input','Analysemodule_output' uitgewenste_processen
        daProcessen = []
        daInput = []
        daOutput = []
        cur.execute("SELECT * FROM %s" % self.table_gewenste_processen)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['instance_fk']:
                if row['instance_fk'] in daInstances:
                    daProcessen.append(row['pk'])
                    if row['analysemodule_input_fk']:
                        daInput.append(row['analysemodule_input_fk'])
                    if row['analysemodule_output_fk']:
                        daOutput.append(row['analysemodule_output_fk'])
            elif row['study_fk']:
                if row['study_fk'] in daStudies:
                    daProcessen.append(row['pk'])
                    if row['analysemodule_input_fk']:
                        daInput.append(row['analysemodule_input_fk'])
                    if row['analysemodule_output_fk']:
                        daOutput.append(row['analysemodule_output_fk'])
            elif row['series_fk']:
                if row['series_fk'] in daSeries:
                    daProcessen.append(row['pk'])
                    if row['analysemodule_input_fk']:
                        daInput.append(row['analysemodule_input_fk'])
                    if row['analysemodule_output_fk']:
                        daOutput.append(row['analysemodule_output_fk'])
        daProcessen = self.uniqifyList(daProcessen)
        daInput = self.uniqifyList(daInput)
        daOutput = self.uniqifyList(daOutput)
        self.destroylist.append( (self.table_gewenste_processen,daProcessen) )
        self.destroylist.append( (self.table_analyseinput,daInput) )
        self.destroylist.append( (self.table_analyseoutput,daOutput) )
        feedback  = "5. Found %d/%d processen for %s" %(len(daProcessen),len(rows_processen),thisthing)
        feedback += "\n5. Found %d/%d inputs for %s" %(len(daInput),len(rows_processen),thisthing)
        feedback += "\n5. Found %d/%d outputs for %s" %(len(daOutput),len(rows_processen),thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        # 6. Gebruik 'Gewenste_processen': voor lijsten 'Resultaten_boolean','Resultaten_char','Resultaten_floating','Resultaten_object' uit
        daBools = []
        daChars = []
        daFloats = []
        daObjects = []
        daStatus = []
        resultlist = [
            (self.table_bool, daBools, 'booleans'),
            (self.table_char, daChars, 'chars'),
            (self.table_float, daFloats, 'floats'),
            (self.table_object, daObjects, 'objects'),
            (self.table_status, daStatus, 'status'),
        ]
        for (tab,lis,desc) in resultlist:
            cur.execute("SELECT * FROM %s" % tab)
            rows_processen = cur.fetchall()
            for row in rows_processen:
                if row['gewenste_processen_fk']:
                    if row['gewenste_processen_fk'] in daProcessen:
                        lis.append(row['pk'])
            lis = self.uniqifyList(lis)
            self.destroylist.append( (tab,lis) )
            feedback = "6. Found %d/%d %s for %s" %(len(lis),len(rows_processen),desc,thisthing)
            if self.qctext:
                self.qctext.appendPlainText(feedback)
            else:
                print feedback

        # 7. list of xml files /dirs to remove
        daXML = []
        cur.execute("SELECT * FROM %s" % self.table_analyseinput)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['pk'] in daInput:
                daXML.append("%s%s"%(row['filepath'],row['filename']))
        feedback = "7. Found %d/%d %s for %s" %(len(daXML),len(rows_processen),"inputXMLs",thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        len0 = len(daXML)
        cur.execute("SELECT * FROM %s" % self.table_analyseoutput)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            if row['pk'] in daOutput:
                daXML.append("%s%s"%(row['filepath'],os.path.dirname(row['filename'])))
        feedback = "7. Found %d/%d %s for %s" %(len(daXML)-len0,len(rows_processen),"outputXMLs",thisthing)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback
        self.destroylist[0]=(None,daXML)

        ### close db
        self.closedb(con)
        QtGui.QApplication.restoreOverrideCursor()

        return

    def onActivatedProcessorChanged(self):
        self.reportentries = {}
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        # 1-2: build list with (gp_is.series ids)
        cur.execute("SELECT * FROM %s" % self.table_gewenste_processen)
        rows_processen = cur.fetchall()

        gp_series_fks = []
        for row in rows_processen:
            if self.selectedSelector==self.selector_id_text['[Selector]'] or self.selectedSelector == row['selector_fk']:
                if self.selectedStatus is None or self.selectedStatus==self.status_id_text['[Status]'] or self.selectedStatus == row['status']:
                    gpid = row['pk']
                    sfk = []
                    if row['series_fk']:
                       sfk.append(row['series_fk'])
                    elif row['study_fk']:
                        cur.execute("SELECT pk FROM %s WHERE study_fk=%d" % (self.table_series,row['study_fk']))
                        rows_series = cur.fetchall()
                        for roww in rows_series:
                            sfk.append(roww['pk'])
                    elif row['instance_fk']:
                        cur.execute("SELECT series_fk FROM %s WHERE pk=%d" % (self.table_instances,row['instance_fk']))
                        rows_instances = cur.fetchall()
                        sfk.append(rows_instances[0]['series_fk']) # can be only one result

                    if len(sfk) >0:
                        gp_series_fks.append([gpid,sfk])

        # 3: check modality and stationname in all series per gp_id; if match, add to reset_list
        reset_list = []
        for gpid,skfs in gp_series_fks:
            for skf in skfs:
                cur.execute("SELECT * FROM %s WHERE pk=%d" % (self.table_series,skf))
                series_rows = cur.fetchall()
                if len(series_rows) == 0:
                    continue
                if self.selectedModality=="[Modality]" or (self.selectedModality == self.stNone and series_rows[0]['modality'] is None) or self.selectedModality == series_rows[0]['modality']: # can be only one result
                    if self.selectedStation=="[StationName]" or (self.selectedStation == self.stNone and series_rows[0]['station_name'] is None) or self.selectedStation == series_rows[0]['station_name']: # can be only one result
                        if series_rows[0]['station_name'] is None:
                            series_rows[0]['station_name'] = self.stNone
                        if series_rows[0]['modality'] is None:
                            series_rows[0]['modality'] = self.stNone
                        if series_rows[0]['src_aet'] is None:
                            series_rows[0]['src_aet'] = self.stNone
                        entry = reportentry(stationname=series_rows[0]['station_name'],modality=series_rows[0]['modality'],gpid=gpid,description=series_rows[0]['series_desc'],seriesid=skf,srcaet=series_rows[0]['src_aet'])

                        if not entry.description:
                            entry.description = ''
                        # add datetime
                        cur.execute("SELECT content_datetime FROM %s WHERE series_fk=%d" % (self.table_instances,skf))
                        instance_rows = cur.fetchall()
                        seriesdate = min ( [ row['content_datetime'] for row in instance_rows] )
                        entry.date = seriesdate

                        if self.runmode == self.stModeReport:
#                            cur.execute("SELECT content_datetime FROM %s WHERE series_fk=%d" % (self.table_instances,skf))
#                            instance_rows = cur.fetchall()
#                            seriesdate = min ( [ row['content_datetime'] for row in instance_rows] )
                            checkseriesdate = datetime.date(year=seriesdate.year,month=seriesdate.month,day=seriesdate.day)
                            if self.selectedReportDateTimestart <= checkseriesdate <= self.selectedReportDateTimeuntil:
                                reset_list.append(gpid)
                                entry.date = seriesdate
                                self.reportentries[gpid] = entry
                        else:
                            reset_list.append(gpid)
                            self.reportentries[gpid] = entry
                        break

        ### close db
        self.closedb(con)

        self.addFilterNames(self.reportentries)

        self.processorShowSelected(reset_list)
        return reset_list

    def onActivatedReportPeriod(self,text):
        if self.verbose:
            print "[onActivatedReportPeriod]",text
        self.selectedReportPeriod = str(text)

        if self.qctable is None or self.statusLabel is None:
            return
        self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        self.onActivatedReportChanged()

    def onActivatedReportQCPeriodicity(self,text):
        if self.verbose:
            print "[onActivatedReportQCPeriodicity]",text
        self.selectedReportQCPeriodicity = str(text)

        if self.qctable is None or self.statusLabel is None:
            return
        self.determineReportPeriod()
        self.onActivatedProcessorChanged()
        self.onActivatedReportChanged()

    def determineReportPeriod(self):
        """
        Workflow:
                1. Find number of days/weeks/months is reporting period
        """
        ## 1. Find number of days/weeks/months is reporting period
        if '[Period to Report]' in self.selectedReportPeriod:
            self.selectedReportDateTimestart = datetime.date.min
            self.selectedReportDateTimeuntil = datetime.date.max
        else:
            i = datetime.date.today()
            y1 = i.year
            m1 = i.month
            d1 = i.day
            if 'Last 3 months' in self.selectedReportPeriod:
                dm = 3
                dy = 0
                d0 = d1
            elif 'Last 6 months' in self.selectedReportPeriod:
                dm = 6
                dy = 0
                d0 = d1
            elif 'Last 12 months' in self.selectedReportPeriod:
                dm = 0
                dy = 1
                d0 = d1
            elif 'Previous year' in self.selectedReportPeriod:
                y1 = i.year -1
                m1 = 12
                dy = 0
                dm = 11
                d0 = 1
                d1 = 31
            elif 'This year' in self.selectedReportPeriod:
                dy = 0
                dm = m1-1
                d0 = 1
            #how many?
            y0 = y1 - dy
            m0 = m1 - dm
            if m0 < 1:
                y0 -= 1
                m0 += 12

            self.selectedReportDateTimestart = datetime.date(year=y0,month=m0,day=d0)
            self.selectedReportDateTimeuntil = datetime.date(year=y1,month=m1,day=d1)

    def onActivatedReportChanged(self):
        """
        Workflow:
                2. Find processes of correct modality, stationname, selector and in given reporting period
                3. Determine expected number of processes
        """
        ## 1. Find number of days/weeks/months is reporting period
        delta_days = (self.selectedReportDateTimeuntil - self.selectedReportDateTimestart).days
        delta_weeks = delta_days/7
        delta_months = self.selectedReportDateTimeuntil.month - self.selectedReportDateTimestart.month+ 12*(self.selectedReportDateTimeuntil.year - self.selectedReportDateTimestart.year)
#        print "days/weeks/months",delta_days, delta_weeks,delta_months

        ## 2. Find processes of correct modality, stationname, selector and in given reporting period

        ## 3. Determine expected number of processes
        if '[' in self.selectedReportPeriod or '[' in self.selectedReportQCPeriodicity:
            return
        num_expected = 0
        if 'Monthly' in self.selectedReportQCPeriodicity:
            num_expected = delta_months
        elif 'Bi-Weekly' in self.selectedReportQCPeriodicity:
            num_expected = delta_weeks/2
        elif 'Weekly' in self.selectedReportQCPeriodicity:
            num_expected = delta_weeks
        elif 'Daily' in self.selectedReportQCPeriodicity:
            num_expected = delta_days
#        print "num expected",num_expected

        if '[' in self.selectedModality and '[' in self.selectedStation and self.selector_id_text['[Selector]'] == self.selectedSelector:
            return

        ## Find number of machines:
        machines = [v.stationname for v in self.reportentries.values()]
        machines =  OrderedDict.fromkeys(machines).keys()

        print 'Expected %d QC entries in %s for each selector of %d machines' %(num_expected,self.selectedReportPeriod,len(machines))
        overal_expected = 0
        overal_total = 0
        overal_notprocessed = 0
        overal_retakes = 0
        overal_errors = 0
        for machine in machines:
            cleanentries = {}
            sel_count = {}
            sel_notprocessed = {}
            sel_retakes = {}
            total_count = 0 # total processables for this machine
            total_notprocessed = 0    # total non-processables for this machine
            total_retakes = 0
            total_errors = 0
            selectors = []
            descriptions = {}
            for sel in self.reportentries.values():
                if sel.stationname == machine:
                    selectors.append(sel.selector)
                    if sel.selector in sel_count:
                        sel_count[sel.selector] += 1
                    else:
                        sel_count[sel.selector] = 1
                    total_count += 1

                    if sel.status == 'Error':
                        if sel.selector in sel_notprocessed:
                            sel_notprocessed[sel.selector] += 1
                        else:
                            sel_notprocessed[sel.selector] = 1
                        total_notprocessed += 1
                        continue

                    # count number of distinct descriptions (filters)
                    if not sel.description in descriptions:
                        descriptions[sel.description] = 1
                    else:
                        descriptions[sel.description] += 1

                    dat = sel.date.strftime('%Y/%m/%d')
                    uid = '%s_%s_%s_%s' %(sel.stationname,sel.selector,dat,sel.description)
                    if not uid in cleanentries:
                        cleanentries[uid] = sel
                    if sel.date > cleanentries[uid].date:
                        print "Retake ",sel#%s: %s %s %s" %(uid,dat,sel.date,cleanentries[uid].date)
                        cleanentries[uid] = sel
                        if sel.selector in sel_retakes:
                            sel_retakes[sel.selector] += 1
                        else:
                            sel_retakes[sel.selector] = 1
                        total_retakes += 1
            print descriptions, len(descriptions)
            selectors =  OrderedDict.fromkeys(selectors).keys()
            for sel in selectors:
                if not sel in sel_count:
                    sel_count[sel] = 0
                if not sel in sel_notprocessed:
                    sel_notprocessed[sel] = 0
                if not sel in sel_retakes:
                    sel_retakes[sel] = 0
            errors,scaleddata = self.checkResults(cleanentries,selectors) # per machine, per selector
            self.graphResults(scaleddata)
            ## Add this selector as section to report
            self.addToReport(machine=machine,selectors=selectors,scaleddata=scaleddata,errors=errors,
                             verwacht=num_expected*len(descriptions),counts=sel_count,bads=sel_notprocessed,retakes=sel_retakes)
            # summary per machine
            total_errors = 0
            if len(selectors) == 0:
                print '%s: %d selectors found in %s' %(machine,len(selectors),self.selectedReportPeriod)
            else:
                print '%s: %d selectors; found %d (%d bad, %d retakes) of %d expected tests: net %.1f%%' \
                      %(machine,len(selectors),total_count,total_notprocessed,total_retakes,
                        num_expected*len(selectors)*len(descriptions),100.*(total_count-total_notprocessed-total_retakes)/(num_expected*len(selectors)*len(descriptions)))
                for sel in selectors:
                    for key,val in errors[sel].iteritems():
                        total_errors += val
                print '     : %d problems:'%(total_errors),
                for sel in selectors:
                    for key,val in errors[sel].iteritems():
                        print key,':',val,';',
                print ''
            reportname=''
            reportname += self.selectedReportDateTimestart.strftime('%Y%m%d')
            reportname += '_'+self.selectedReportDateTimeuntil.strftime('%Y%m%d')
            reportname += '_'+machine+'.pdf'
            self.finishReport(reportname,machine=machine,verwacht=num_expected*len(descriptions),count=total_count,
                              bads=total_notprocessed,retakes=total_retakes,errors=total_errors)

            # summary over all machines
            overal_expected += num_expected*len(selectors)*len(descriptions)
            overal_total += total_count
            overal_notprocessed += total_notprocessed
            overal_retakes += total_retakes
            overal_errors += total_errors

        if overal_expected == 0:
            print 'Overal: expected %d tests in %s' %( overal_expected,self.selectedReportPeriod)
        else:
            print 'Overal: found %d (%d bad, %d retakes) of %d expected tests (net %.1f%% done) showing %d errors' \
                  %( overal_total,overal_notprocessed,overal_retakes,overal_expected,100.*(overal_total-overal_notprocessed-overal_retakes)/overal_expected,overal_errors )

    def graphResults(self,scaleddata):
        """
        Expected input: dictionary of all correctly processed entries for one stationname
        For each selector, make a graph of all scaled values if limits given

        """
        if self.qcgraph:
            self.qcgraph.clear()
            if self.qcgraph.getPlotItem().legend:
                self.qcgraph.getPlotItem().legend.items = []
            else:
                self.qcgraph.addLegend()
            self.qcgraph.setXRange(
                datetime2seconds(self.selectedReportDateTimestart),
                datetime2seconds(self.selectedReportDateTimeuntil),
                padding=.0 )
            self.qcplots = {}

        if self.qcgroupbox:
            self.clearLayout(self.qcgroupbox.layout(),layoutOnly=True)

        ix = 0
        linewidth = 2.
        boxcols = 4
        for sel in scaleddata.keys():
            if self.qcgraph:
                self.qcgraph.setTitle(sel)
            selector = scaleddata[sel]
            mintime = None
            maxtime = None
            for name in sorted(selector.keys()):
                data = sorted(selector[name])
                if mintime:
                    mintime = min(mintime,data[0][0])
                    maxtime = max(maxtime,data[-1][0])
                else:
                    mintime = data[0][0]
                    maxtime = data[-1][0]

                penstyle = 1
                if ix < 9*5:
                    penstyle = ix/9+1


                if self.qcgraph: # and 'doseratio' in name:
                    if len(data) == 1:
                        self.qcplots[name] = self.qcgraph.plot(zip(*data)[0],zip(*data)[1], symbol='o',name=name,symbolPen=pg.mkPen(width=linewidth, color=pg.intColor(ix)))
                    else:
                        self.qcplots[name] = self.qcgraph.plot(zip(*data)[0],zip(*data)[1],name=name,pen=pg.mkPen(width=linewidth, color=pg.intColor(ix),style=penstyle))
                    if self.qcgroupbox:
                        nwButton = QtGui.QCheckBox(name)
                        nwButton.setChecked(True)
                        nwButton.stateChanged.connect(functools.partial(self.togglePlot,name))
                        self.qcgroupbox.layout().addWidget(nwButton,ix/boxcols,ix%boxcols)
                    ix += 1

            if self.qcgraph: # and 'doseratio' in name:
                limitpen = pg.mkPen(width=linewidth*2, color=(255,255,255),style=QtCore.Qt.DashLine)
                name = 'CritMAX'
                self.qcplots[name] = self.qcgraph.plot([mintime,maxtime],[1.,1.],pen=limitpen,name='CritMAX')
                name = 'CritMIN'
                self.qcplots[name] = self.qcgraph.plot([mintime,maxtime],[-1.,-1.],pen=limitpen,name='CritMIN')

    def togglePlot(self,name,state):
        if state == QtCore.Qt.Checked:
            self.qcplots[name].show()
            for sample,label in self.qcgraph.getPlotItem().legend.items:
                if label.text == name:
                    label.show()
                    sample.show()
                    break
        else:
            self.qcplots[name].hide()
            for sample,label in self.qcgraph.getPlotItem().legend.items:
                if label.text == name:
                    label.hide()
                    sample.hide()
                    break

    def addFilterNames(self,entries,doubleid=False):
        # If an extra identifier (like Filter or Antomy) is given use that instead of description
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        gps = {}
        for key,entry in entries.iteritems():
            if not entry.gpid in gps:
                gps[entry.gpid] = []
            gps[entry.gpid].append(key)

        cur.execute("SELECT * FROM %s WHERE omschrijving LIKE '%%%s%%'" % (self.table_char,'FilterMaterial'))
        rowsF = cur.fetchall()
        cur.execute("SELECT * FROM %s WHERE omschrijving LIKE '%%%s%%'" % (self.table_char,'Anatomy'))
        rowsP = cur.fetchall()

        fnames = {}
        for row in rowsF:
            uid = row['gewenste_processen_fk']
            if not uid in fnames:
                fnames[uid] = []
            fnames[uid].append(row['waarde'])
        for row in rowsP:
            uid = row['gewenste_processen_fk']
            if not uid in fnames:
                fnames[uid] = []
            fnames[uid].append(row['waarde'])


        for uid,vals in fnames.iteritems():
            if uid in gps:
                for en in gps[uid]:
                    entries[en].description = ' '+'#'.join(sorted(vals))

        ### close db
        self.closedb(con)


    def checkResults(self,entries,selectors=None):
        """
        Check all results of given entries, identified by gewensteprocessen_id.
        If any of the tables resultaten_char, resultaten_floating is outside specs, put in report
        """
        if selectors:
            ## 1. for each selector make dictionary of number of errors
            errors = {}
            ## 2. for each selector make dictionary of scaled measurements vs date
            scaleddata = {}
            for sel in selectors:
                errors[sel] = {}
                scaleddata[sel] = {}

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        for en in entries.values():
            if en.status == 'Error':
                en.analysisresult = self.stResultError
                continue
            hascrits = False
            en.analysisresult = self.stResultNoRanges
            sel = en.selector
            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (self.table_char,en.gpid))
            rows = cur.fetchall()
            for row in rows:
                if row['criterium']:
                    hascrits = True
                    if selectors and (row['waarde'] != row['criterium']):
                        if row['omschrijving'] in errors[sel]:
                            errors[sel][row['omschrijving']] += 1
                        else:
                            errors[sel][row['omschrijving']] = 1

            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (self.table_float,en.gpid))
            rows = cur.fetchall()
            for row in rows:
                krmax = row['grens_kritisch_boven']
                krmin = row['grens_kritisch_onder']
                if not krmax is None and not krmin is None:
                    hascrits = True
                    if selectors:
                        measurement_value = row['waarde']
                        measurement_name = row['omschrijving']
                        if measurement_value is None:
                            continue
                        scaled_value = (measurement_value-(krmax+krmin)/2.)/((krmax-krmin)/2.)
                        if not measurement_name in scaleddata[sel]:
                            scaleddata[sel][measurement_name]=[]
                        scaleddata[sel][measurement_name].append( (datetime2seconds(en.date),scaled_value,measurement_value) )

                        if scaled_value<-1 or scaled_value>1:
                            en.analysisresult = self.stResultOutOfRange
                            if measurement_name in errors[sel]:
                                errors[sel][measurement_name] += 1
                            else:
                                errors[sel][measurement_name] = 1
            if en.analysisresult == self.stResultNoRanges and hascrits:
                en.analysisresult = self.stResultOK
        ### close db
        self.closedb(con)
        if selectors:
            return errors,scaleddata
        return



    def resetListedProcessor(self):
        if not self.reportentries:
            return

        if self.qctable is None:
            return

        if self.qctable.rowCount() == 0:
            return

        if self.qctable.columnCount() == 0:
            return

        reset_list = []
        for en in self.reportentries.values():
            reset_list.append( en.gpid )

        if len(reset_list) == 0:
            return

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        # 4: for all gp_id in reset_list, set status to "reset"
        nrows = 0
        for gpid in reset_list:
            cur.execute("UPDATE %s set status=%d where %s.pk=%d" % (self.table_gewenste_processen,self.status_reset,self.table_gewenste_processen,gpid))
            nrows += cur.rowcount

        self.statusLabel.setText("Did reset of %d processes" %nrows)
        con.commit()

        ### close db
        self.closedb(con)

    def mysqlDeleteRows(self,table,value_list,key='pk'):
        #DELETE FROM your_table WHERE id_users=1 AND id_product=2
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        # for all values of given key in given table do delete row
        nrows = 0
        for val in value_list:
            cur.execute("DELETE FROM %s WHERE %s=%s" % (table,key,str(val)))
            nrows += cur.rowcount

        self.statusLabel.setText("Did delete %d rows of table %s" %(nrows,table))
#        print "Did delete %d rows of table %s" %(nrows,table)
        con.commit()

        ### close db
        self.closedb(con)

    def xmlDeleteFandF(self,dalist):
        error = 1
        path = self.pathXMLInputOutput()
        if path is None:
            print "ERROR! Cannot determine root path for XML files"
            return error

        outname = os.path.join(os.getcwd(),"dbtool_remove_list.txt")
        with open(outname, 'w') as myfile:
            for fname in dalist:
                myfile.write("rm -rf '%s'\n" %(os.path.join(path,fname)))
        try:
            os.popen("pkexec /bin/bash %s" %outname)
        except:
            print "You need to install the policykit to use this command: apt-get install policykit-1"
            return error

        error = 0
        feedback = "Deleted %d XML input files and XML output folders." %len(dalist)
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback

        return error

    def destroySelected(self):
        if len(self.destroylist)<1:
            return

        QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        if not self.destroylist[0][0]: # None is reserved for XML filelist
            error = self.xmlDeleteFandF(self.destroylist[0][1])
        else:
            feedback = "ERROR! Nothing deleted, destroylist is malformed."
            if self.qctext:
                self.qctext.appendPlainText(feedback)
            else:
                print feedback
            QtGui.QApplication.restoreOverrideCursor()
            return

        if error>0:
            feedback = "ERROR! Nothing deleted, XML files could not be removed."
            if self.qctext:
                self.qctext.appendPlainText(feedback)
            else:
                print feedback
            QtGui.QApplication.restoreOverrideCursor()
            return

        num =0
        for tab,lis in self.destroylist:
            if tab:
                num += len(lis)
                self.mysqlDeleteRows(tab,lis)

        feedback = "Deleted %d rows in iqc database tables." %num
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback
        self.statusLabel.setText("%s" %(feedback))

        self.destroylist = []
        feedback = "Done! all traces of %s gone from the iqc database. DCM4CHEE database is unaltered."%self.selectedStation
        if self.qctext:
            self.qctext.appendPlainText(feedback)
        else:
            print feedback
        QtGui.QApplication.restoreOverrideCursor()

    def seekandDestroy(self,mode):
        self.runmode = self.stModeDestroy

        # clear GUI
        self.clearLayout(self.layout)

        # table showing selection
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)
        self.qctable.horizontalHeader().sectionClicked.connect(self.qctable.sortByColumn)
        self.qctext = QtGui.QPlainTextEdit()
        self.qctext.setReadOnly(True) # disable editing
        self.qctext.setTextInteractionFlags(self.qctext.textInteractionFlags() | QtCore.Qt.TextSelectableByKeyboard) # keep selectable

        # Reset Button
        btnDestroy = QtGui.QPushButton("Destroy all traces of selected StationName/ID/Selector/Modality")
        btnDestroy.clicked.connect(self.destroySelected)

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)
        validbdversions = ['1.1.0']
        if not dbversion in validbdversions:
            self.qctext.clear()
            self.qctext.appendPlainText("Database Version %s is not supported for Search & Destroy!" % dbversion)
            self.layout.addWidget(self.qctext,0,0)
            self.qw.setLayout(self.layout)
            ### close db
            self.closedb(con)
            return

        if mode == 'GewenstID':
            # dropdown with gewenste_processen ids
            cmGewenst = QtGui.QComboBox(self)
            cmGewenst.addItem("[GewenstID]")
            self.selectedGewenst = cmGewenst.itemText(0)
            selected = []
            cur.execute("select pk from %s" % (self.table_gewenste_processen))
            rows = cur.fetchall()
            for row in rows:
                if row['pk'] is None:
                    selected.append(self.stNone)
                else:
                    selected.append(row['pk'])
            names = sorted(OrderedDict.fromkeys(selected).keys())
            for ho in names:
                cmGewenst.addItem(str(ho))
            cmGewenst.activated[str].connect(self.onActivatedGewenstProces)
        elif mode == 'Selector':
            # dropdown with gewenste_processen ids
            cmSelectors = QtGui.QComboBox(self)
            cmSelectors.addItem("[Selector]")
            self.selectedSelector = cmSelectors.itemText(0)
            self.selector_id_text = {}
            selected = []
            cur.execute("select * from %s" % (self.table_selector))
            rows = cur.fetchall()
            minid = None
            for row in rows:
                selected.append(row['name'])
                if minid is None:
                    minid = row['pk']
                else:
                    minid = min (minid,row['pk'])
                self.selector_id_text[row['name']] = row['pk']
                self.selector_id_text[row['pk']] = row['name']
            names = sorted(OrderedDict.fromkeys(selected).keys())
            for ho in names:
                cmSelectors.addItem(ho)
            cmSelectors.activated[str].connect(self.onActivatedProcessorSelector)
        elif mode == 'Modality':
            # dropdown with modalities
            cmModality = QtGui.QComboBox(self)
            cmModality.addItem("[Modality]")
            self.selectedModality = cmModality.itemText(0)
            ### fill ComboBoxes
            selected = []
            cur.execute("select modality from %s" % (self.table_series))
            rows = cur.fetchall()
            for row in rows:
                if row['modality'] is None:
                    selected.append(self.stNone)
                else:
                    selected.append(row['modality'])
            names = sorted(OrderedDict.fromkeys(selected).keys())
            for ho in names:
                cmModality.addItem(ho)
            cmModality.activated[str].connect(self.onActivatedProcessorModality)
        elif mode == 'StationName':
            # dropdown with stationnames
            cmStations = QtGui.QComboBox(self)
            cmStations.addItem("[StationName]")
            self.selectedStation = cmStations.itemText(0)
            ### fill ComboBoxes
            selected = []
            cur.execute("select station_name from %s" % (self.table_series))
            rows = cur.fetchall()
            for row in rows:
                if row['station_name'] is None:
                    selected.append(self.stNone)
                else:
                    selected.append(row['station_name'])
            names = sorted(OrderedDict.fromkeys(selected).keys())
            for ho in names:
                cmStations.addItem(ho)
            cmStations.activated[str].connect(self.onActivatedProcessorStation)
        else:
            # dropdown with SrcAET
            cmSrcAET = QtGui.QComboBox(self)
            cmSrcAET.addItem("[SrcAET]")
            self.selectedAET = cmSrcAET.itemText(0)
            ### fill ComboBoxes
            selected = []
            cur.execute("select src_aet from %s" % (self.table_series))
            rows = cur.fetchall()
            for row in rows:
                if row['src_aet'] is None:
                    selected.append(self.stNone)
                else:
                    selected.append(row['src_aet'])
            names = sorted(OrderedDict.fromkeys(selected).keys())
            for ho in names:
                cmSrcAET.addItem(ho)
            cmSrcAET.activated[str].connect(self.onActivatedProcessorSrcAET)

        ### close db
        self.closedb(con)

        ### build GUI
        if mode == 'GewenstID':
            self.layout.addWidget(cmGewenst,0,0)
        elif mode == 'Selector':
            self.layout.addWidget(cmSelectors,0,0)
        elif mode == 'Modality':
            self.layout.addWidget(cmModality,0,0)
        elif mode == 'StationName':
            self.layout.addWidget(cmStations,0,1)
        else:
            self.layout.addWidget(cmSrcAET,0,1)

        self.layout.addWidget(self.qctext,1,0,1,2)
        self.layout.addWidget(btnDestroy,2,1)
        self.qw.setLayout(self.layout)

    def processorReset(self):
        self.runmode = self.stModeProcessor
        # clear GUI
        self.clearLayout(self.layout)

        # table showing selection
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)
        self.qctable.horizontalHeader().sectionClicked.connect(self.qctable.sortByColumn)

        # dropdown with modalities, stationnames, selector
        cmModality = QtGui.QComboBox(self)
        cmModality.addItem("[Modality]")
        self.selectedModality = cmModality.itemText(0)

        cmStations = QtGui.QComboBox(self)
        cmStations.addItem("[StationName]")
        self.selectedStation = cmStations.itemText(0)

        zerotext = "[Selector]"
        cmSelectors = QtGui.QComboBox(self)
        cmSelectors.addItem(zerotext)

        zerotext = "[Status]"
        cmStatus = QtGui.QComboBox(self)
        cmStatus.addItem(zerotext)

        # Reset Button
        btReset = QtGui.QPushButton("Reset processor status for\nprocesses shown in table")
        btReset.clicked.connect(self.resetListedProcessor)

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        ### fill ComboBoxes
        selected = []
        cur.execute("select modality from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            if row['modality'] is None:
                selected.append(self.stNone)
            else:
                selected.append(row['modality'])
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmModality.addItem(ho)
        cmModality.activated[str].connect(self.onActivatedProcessorModality)

        selected = []
        cur.execute("select station_name from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            if row['station_name'] is None:
                selected.append(self.stNone)
            else:
                selected.append(row['station_name'])
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmStations.addItem(ho)
        cmStations.activated[str].connect(self.onActivatedProcessorStation)

        self.selector_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_selector))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['name'])
            if minid is None:
                minid = row['pk']
            else:
                minid = min (minid,row['pk'])
            self.selector_id_text[row['name']] = row['pk']
            self.selector_id_text[row['pk']] = row['name']
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmSelectors.addItem(ho)
        zerotext = str(cmSelectors.itemText(0))
        self.selector_id_text[zerotext] = minid-1
        self.selector_id_text[minid-1] = zerotext
        self.selectedSelector = self.selector_id_text[zerotext]
        cmSelectors.activated[str].connect(self.onActivatedProcessorSelector)

        self.status_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_status_omschrijving))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['veld_omschrijving'])
            self.status_id_text[row['veld_omschrijving']] = row['nummer']
            self.status_id_text[row['nummer']] = row['veld_omschrijving']
            if minid is None:
                minid = row['nummer']
            else:
                minid = min (minid,row['nummer'])
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmStatus.addItem(ho)
        zerotext = str(cmStatus.itemText(0))
        self.status_id_text[zerotext] = minid-1
        self.status_id_text[minid-1] = zerotext
        self.selectedStatus = self.status_id_text[zerotext]
        cmStatus.activated[str].connect(self.onActivatedProcessorStatus)

        ### select all processes
        self.onActivatedProcessorSelector(str(cmSelectors.itemText(0)))


        ### close db
        self.closedb(con)

        ### build GUI
        self.layout.addWidget(cmModality,0,0)
        self.layout.addWidget(cmStations,0,1)
        self.layout.addWidget(cmSelectors,1,0)
        self.layout.addWidget(cmStatus,1,1)

        self.layout.addWidget(self.qctable,2,0,1,2)
        self.layout.addWidget(btReset,3,1)
        self.qw.setLayout(self.layout)

    def databaseTruncatePopUp(self):
        # clear GUI
        self.clearLayout(self.layout)
        qctext = "Did NOT truncate ANY tables"

        title = "Warning"
        text = "Are you sure you want to truncate?\n" \
               "This will remove all datasets and results from the IQC database,\n" \
               "which will trigger the WAD-Collector to treat all data in the DCM4CHEE database as new\n\n" \
               "This is NOT a clean reset of the WAD-Server.\n" \
               "XML input and output files are NOT removed (but new ones will generated, so expect clutter)\n" \
               "Configuration files and Selectors will remain and stay active.\n\n" \
               "Last chance: press \"OK\" to truncate or \"CANCEL\" to bail out."

        ret = QtGui.QMessageBox.warning(self, title,text,
                QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Ok)
        if ret == QtGui.QMessageBox.Ok:
            # textfield showing number of selected processes
            nrows = self.databaseTruncate()
            qctext = "Did truncate %d tables" %nrows

        self.statusLabel.setText(qctext)

    def databaseTruncate(self):
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        # truncating
        cur.execute("SET FOREIGN_KEY_CHECKS=0")
        nrows = 0
        for ta in self.truncate_tables:
            cur.execute("truncate %s" % ta)
            nrows += 1
        cur.execute("SET FOREIGN_KEY_CHECKS=1")

        con.commit()

        ### close db
        self.closedb(con)
        return nrows

    def popupPurge(self):
        # clear GUI
        self.clearLayout(self.layout)
        qctext = "Did NOT purge!"

        path = self.pathXMLInputOutput()
        if path is None:
            return
        path = os.path.join(path,'XML')
        title = "Warning"
        text = "Are you sure you want to purge?\n" \
               "This will remove the files in %s\n" \
               "and remove all datasets and results from the IQC database,\n" \
               "which will trigger the WAD-Collector to treat all data in the DCM4CHEE database as new\n\n" \
               "This is NOT a clean reset of the WAD-Server.\n" \
               "Configuration files and Selectors will remain and stay active.\n\n" \
               "Last chance: press \"OK\" to purge or \"CANCEL\" to bail out." % path

        ret = QtGui.QMessageBox.warning(self, title,text,
                QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Ok)
        if ret == QtGui.QMessageBox.Ok:
            path1 = os.path.join(path,'analysemodule_input')
            path2 = os.path.join(path,'analysemodule_output')
            cpt = sum([len(files) for r, d, files in os.walk(path1)])
            cpt += sum([len(files) for r, d, files in os.walk(path2)])
            path1 = os.path.join(path1,'*')
            path2 = os.path.join(path2,'*')
            os.popen("pkexec rm -rf %s %s" %(path1,path2))
            trows = self.databaseTruncate()
            qctext = "Did truncate %d tables and deleted %d files" % (trows,cpt)

        self.statusLabel.setText(qctext)


    def pathXMLInputOutput(self):
        if self.xmlroot: # provided in ini file
            return self.xmlroot

        path = None
        # print in dialog: only works if on db (and server) on localhost, and installed as services
        if os.path.isfile('/etc/init/WAD-Selector.conf'):
            with open('/etc/init/WAD-Selector.conf') as infile:
                rows = infile.readlines()
            for row in rows:
                if 'chdir' in row:
                    try:
                        path = os.path.dirname( (row.split(' ')[1]).split('WAD_Services')[0] )
                        break
                    except:
                        path = None
        elif os.path.isfile('/etc/init.d/WAD-Selector'):
            with open('/etc/init.d/WAD-Selector') as infile:
                rows = infile.readlines()
            for row in rows:
                if 'RUN_FROM' in row:
                    try:
                        path = os.path.dirname( (row.split('=')[1]).split('WAD_Services')[0] )
                        break
                    except:
                        path = None
        else:
            print "[removeXMLInputOutput] not on localhost or WAD services not installed as services"

        if path is None:
            print "[removeXMLInputOutput] Could not determine path of WAD input/output"
        return path


    def periodicStatus(self):
        self.runmode = self.stModeReport
        # clear GUI
        self.clearLayout(self.layout)
        self.status_id_text = {}

        # table showing results
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)
        self.qctable.horizontalHeader().sectionClicked.connect(self.qctable.sortByColumn)
        # graph to viz data
        dateaxis = DateAxis(orientation='bottom')
        self.qcgraph = pg.PlotWidget(title="Report", axisItems={'bottom': dateaxis})
        self.qcgraph.setLabel('left','scaled measurement')
        self.qcgraph.setTitle('Report')
        self.qcgraph.showGrid(y=True, alpha=0.85)
        # buttonboxes
        self.qcgroupbox = QtGui.QGroupBox(self,title="Show graphs")
        self.qcgroupbox.setLayout(QtGui.QGridLayout())

        # dropdown with period of report and QC periodicity
        cmReportPeriod = QtGui.QComboBox(self)
        items = [
            "[Period to Report]",
            "Last 3 months",
            "Last 6 months",
            "Last 12 months",
            "Previous year",
            "This year",
            ]
        for it in items:
            cmReportPeriod.addItem(it)
        self.selectedReportPeriod = cmReportPeriod.itemText(0)
        cmReportPeriod.activated[str].connect(self.onActivatedReportPeriod)

        cmQCPeriodicity = QtGui.QComboBox(self)
        items = [
            "[QC Periodicity]",
            "Monthly",
            "Bi-Weekly",
            "Weekly",
            "Daily"
            ]
        for it in items:
            cmQCPeriodicity.addItem(it)
        self.selectedReportQCPeriodicity = cmQCPeriodicity.itemText(0)
        cmQCPeriodicity.activated[str].connect(self.onActivatedReportQCPeriodicity)

        # dropdown with modalities, stationnames, selector
        cmModality = QtGui.QComboBox(self)
        cmModality.addItem("[Modality]")
        self.selectedModality = cmModality.itemText(0)

        cmStations = QtGui.QComboBox(self)
        cmStations.addItem("[StationName]")
        self.selectedStation = cmStations.itemText(0)

        zerotext = "[Selector]"
        cmSelectors = QtGui.QComboBox(self)
        cmSelectors.addItem(zerotext)

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        ### fill ComboBoxes
        selected = []
        cur.execute("select modality from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            selected.append(row['modality'])
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmModality.addItem(ho)
        cmModality.activated[str].connect(self.onActivatedProcessorModality)

        selected = []
        cur.execute("select station_name from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            if row['station_name'] is None:
                selected.append(self.stNone)
            else:
                selected.append(row['station_name'])
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmStations.addItem(ho)
        cmStations.activated[str].connect(self.onActivatedProcessorStation)

        self.selector_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_selector))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['name'])
            if minid is None:
                minid = row['pk']
            else:
                minid = min (minid,row['pk'])
            self.selector_id_text[row['name']] = row['pk']
            self.selector_id_text[row['pk']] = row['name']
        names = sorted(OrderedDict.fromkeys(selected).keys())
        for ho in names:
            cmSelectors.addItem(ho)
        zerotext = str(cmSelectors.itemText(0))
        self.selector_id_text[zerotext] = minid-1
        self.selector_id_text[minid-1] = zerotext
        self.selectedSelector = self.selector_id_text[zerotext]
        cmSelectors.activated[str].connect(self.onActivatedProcessorSelector)

        self.status_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_status_omschrijving))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['veld_omschrijving'])
            self.status_id_text[row['veld_omschrijving']] = row['nummer']
            self.status_id_text[row['nummer']] = row['veld_omschrijving']
            if minid is None:
                minid = row['nummer']
            else:
                minid = min (minid,row['nummer'])
        names = OrderedDict.fromkeys(selected).keys()
        zerotext = "[Status]"
        self.status_id_text[zerotext] = minid-1
        self.status_id_text[minid-1] = zerotext
        self.selectedStatus = self.status_id_text[zerotext]

        ### close db
        self.closedb(con)

        self.onActivatedReportPeriod(cmReportPeriod.itemText(0))
        ### build GUI
        self.layout.addWidget(cmReportPeriod,0,0)
        self.layout.addWidget(cmQCPeriodicity,0,1)

        self.layout.addWidget(cmModality,1,0)
        self.layout.addWidget(cmStations,1,1)
        self.layout.addWidget(cmSelectors,1,2)

        self.layout.addWidget(self.qctable,2,0,1,3)
        self.layout.addWidget(self.qcgraph,3,0,1,3)
        self.layout.addWidget(self.qcgroupbox,4,0,1,3)
        self.qw.setLayout(self.layout)

    def finishReport(self,reportfile,machine='',verwacht=0,count=0,bads=0,retakes=0,errors=0):
        if self.reporter is None:
            return

        abstract_data = []
        abstract_data.append( ['Machine:', '%s'%machine] )
        abstract_data.append( ['Periode:', '%s t/m %s'%(self.selectedReportDateTimestart.strftime('%Y-%m-%d'),self.selectedReportDateTimeuntil.strftime('%Y-%m-%d'))] )
        abstract_data.append( ['Verwachte #testen:', '%d'%verwacht] )
        abstract_data.append( ['Uitgevoerde #testen:','%d (%d onbruikbaar, %d herhalingen)'%(count,bads,retakes)] )
        if verwacht == 0:
            uitgevoerd = 100.
        else:
            uitgevoerd = 100.*(count-bads-retakes)/verwacht
        abstract_data.append( ['Netto score:','%.1f%% uitgevoerd' %uitgevoerd] )
        abstract_data.append( ['Gevonden problemen:','%d'%errors] )

        self.reporter.addtitlepagewithtable(abstract_data)

        self.reporter.render(reportfile)
        self.statusLabel.setText("Report written to %s." %reportfile)
        self.reporter = None

    def addToReport(self,machine='',selectors=[],scaleddata={},errors={},
                             verwacht=0,counts={},bads={},retakes={}):
        if self.reporter is None:
            self.reporter = Reporter(title='QC Status Rapportage',author='dbTool.py versie %d'%self.qcversion,url=self.host)

        for sel in sorted(selectors):
            section_title   = self.reporter.header(machine+': '+sel)
            section_explain1 = self.reporter.paragraph("Periode: %s t/m %s" %(self.selectedReportDateTimestart.strftime('%Y-%m-%d'),self.selectedReportDateTimeuntil.strftime('%Y-%m-%d')))
            section_explain2 = self.reporter.paragraph("Verwachte #testen: %d" % verwacht)
            section_explain3 = self.reporter.paragraph("Uitgevoerde #testen: %d (%d onbruikbaar, %d herhalingen)" %(counts[sel],bads[sel],retakes[sel]))
            section_explain4 = self.reporter.paragraph("Netto score: %.1f%% uitgevoerd" %(100.*(counts[sel]-bads[sel]-retakes[sel])/verwacht) )

            total_errors = 0
            table_data = []
            table_data.append(['meting','#problemen'])
            for key,val in sorted(errors[sel].iteritems()):
                table_data.append( [key,val] )
                total_errors += val

            section_tabletext = self.reporter.paragraph("In totaal zijn %d problemen gevonden." % total_errors)
            section_table = self.reporter.table( table_data )

            graph_data = []
            for name in sorted(scaleddata[sel].keys()):
                dadata = sorted(scaleddata[sel][name])
                xyarr = [ (da[0],da[1]) for da in dadata ]
                graph_data.append( (name,xyarr) )
            if len(graph_data) == 0:
                continue
            if len(graph_data)>self.reporter.maxnumgraphsinchart: # max two charts to this page
                subset = graph_data[0:min(2*self.reporter.maxnumgraphsinchart,len(graph_data))]
                section_chart = self.reporter.stackedGraphs(subset,tmin=datetime2seconds(self.selectedReportDateTimestart),tmax=datetime2seconds(self.selectedReportDateTimeuntil))
                section = [section_title, section_explain1,section_explain2,section_explain3,section_explain4,section_tabletext, section_table, section_chart]
                self.reporter.addsection(section)

                if len(subset)<len(graph_data):
                    graph_data = graph_data[len(subset):]
                    if len(graph_data)>4*self.reporter.maxnumgraphsinchart:
                        print "ERROR! Still too many items to plot!",len(graph_data)
                        return
                    else:
                        section_chart = self.reporter.stackedGraphs(graph_data,tmin=datetime2seconds(self.selectedReportDateTimestart),tmax=datetime2seconds(self.selectedReportDateTimeuntil))
                        self.reporter.addsection(section_chart)
            else:
                section_chart = self.reporter.dategraph(graph_data,tmin=datetime2seconds(self.selectedReportDateTimestart),tmax=datetime2seconds(self.selectedReportDateTimeuntil))
                section = [section_title, section_explain1,section_explain2,section_explain3,section_explain4,section_tabletext, section_table, section_chart]
                self.reporter.addsection(section)


    def currentStatus(self,mode='selector+description'):
        """
        Shows the date and outcome of the last run of each selector
        Flowchart:
        1. Get list of all selectors
        2. For each selector find last run
        3. Of each last run, determine result if NOT status == error
        :return:
        """
        self.runmode = self.stModeReport
        # clear GUI
        self.clearLayout(self.layout)
        self.status_id_text = {}

        # table showing results
        # table showing results
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)
        self.qctable.horizontalHeader().sectionClicked.connect(self.qctable.sortByColumn)
        self.qcgraph = None
        self.qcgroupbox = None

        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        ### fill ComboBoxes
        selected = []
        cur.execute("select modality from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            selected.append(row['modality'])
        modalities = sorted(OrderedDict.fromkeys(selected).keys())

        selected = []
        cur.execute("select station_name from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            if row['station_name'] is None:
                selected.append(self.stNone)
            else:
                selected.append(row['station_name'])
        stations = sorted(OrderedDict.fromkeys(selected).keys())

        self.selector_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_selector))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['name'])
            if minid is None:
                minid = row['pk']
            else:
                minid = min (minid,row['pk'])
            self.selector_id_text[row['name']] = row['pk']
            self.selector_id_text[row['pk']] = row['name']
        selectors = sorted(OrderedDict.fromkeys(selected).keys())

        self.status_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_status_omschrijving))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['veld_omschrijving'])
            self.status_id_text[row['veld_omschrijving']] = row['nummer']
            self.status_id_text[row['nummer']] = row['veld_omschrijving']
            if minid is None:
                minid = row['nummer']
            else:
                minid = min (minid,row['nummer'])
        statusi = OrderedDict.fromkeys(selected).keys()

        allentries = {}
        # 2.1 first make list linking all gp_ids to series_ids and add analysisdate and status
        cur.execute("SELECT * FROM %s" % self.table_gewenste_processen)
        rows_processen = cur.fetchall()
        for row in rows_processen:
            gpid = row['pk']
            sfk = []
            if row['series_fk']:
                sfk.append(row['series_fk'])
            elif row['study_fk']:
                cur.execute("SELECT pk FROM %s WHERE study_fk=%d" % (self.table_series,row['study_fk']))
                rows_series = cur.fetchall()
                for roww in rows_series:
                    sfk.append(roww['pk'])
            elif row['instance_fk']:
                cur.execute("SELECT series_fk FROM %s WHERE pk=%d" % (self.table_instances,row['instance_fk']))
                rows_instances = cur.fetchall()
                sfk.append(rows_instances[0]['series_fk']) # can be only one result

            for sk in sfk: # WIJKT AF VAN EERDER, OMDAT ONDERSCHEID MEERDERE SERIES PER STUDY NODIG
                allentries[str(gpid)+'_'+str(sk)]= reportentry(gpid=gpid,seriesid=sk,
                                                      status=self.status_id_text[row['status']],
                                                      analysisdate=row['creation_time'],
                                                      selector=self.selector_id_text[row['selector_fk']],
                                                      description='')

        # 2.2: Add modality and stationname to all entries
        for gpid,entry in allentries.iteritems():
            cur.execute("SELECT * FROM %s WHERE pk=%d" % (self.table_series,entry.seriesid))
            series_rows = cur.fetchall()
            if series_rows[0]['station_name'] is None:
                series_rows[0]['station_name'] = self.stNone
            if series_rows[0]['modality'] is None:
                series_rows[0]['modality'] = self.stNone
            if series_rows[0]['series_desc'] is None:
                series_rows[0]['series_desc'] = ''
            if series_rows[0]['src_aet'] is None:
                series_rows[0]['src_art'] = self.stNone
            entry.description += series_rows[0]['series_desc']
            entry.stationname=series_rows[0]['station_name']
            entry.modality=series_rows[0]['modality']
            entry.srcaet = series_rows[0]['src_aet']

            # add content_datetime to entry
            cur.execute("SELECT content_datetime FROM %s WHERE series_fk=%d" % (self.table_instances,entry.seriesid))
            instance_rows = cur.fetchall()
            seriesdate = min ( [ row['content_datetime'] for row in instance_rows] )
            entry.date = seriesdate

        ### close db
        self.closedb(con)

        self.addFilterNames(allentries,doubleid=True) # just make sure both RH and MO are different series

        # 3. make list of only the latest series given to each selector
        self.reportentries = {}
        for key,entry in allentries.iteritems():
            uid = entry.selector
            if mode == 'selector+description':
                uid+='_'+entry.description

            if uid in self.reportentries:
                if entry.date > self.reportentries[uid].date:
                    self.reportentries[uid] = entry
            else:
                self.reportentries[uid] = entry

        # 4. checkresults
        self.checkResults(self.reportentries)

        # 5. show in table
        headers = ['Selector','StationName','Description','SrcAET','AcquisitionDate','WADStatus','AnalysisDate','Result']
        self.qctable.setRowCount(len(self.reportentries))
        if len(self.reportentries) > 0:
            self.qctable.setColumnCount(len(headers))
            rowid = 0
            self.qctable.setHorizontalHeaderLabels(headers)
            for key,entry in sorted(self.reportentries.iteritems()):
                colid = 0
                vals = [entry.selector, entry.stationname,entry.description,entry.srcaet,entry.date,entry.status,entry.analysisdate,entry.analysisresult]
                for v in vals:
                    self.qctable.setItem(rowid, colid, QtGui.QTableWidgetItem(str(v)))
                    colid += 1
                rowid += 1

        self.statusLabel.setText("Done")

        ### build GUI
        self.layout.addWidget(self.qctable,0,0)
        self.qw.setLayout(self.layout)


    def reportDemo(self):
        from numpy import random
        randn = random.randn
        ## 2. for each selector make dictionary of scaled measurements vs date
        machine = 'AZUMAMM'
        selectors = ['sel_a','sel_b','sel_c']
        measurements = ['meas_a','meas_b','meas_c','meas_d','meas_e']
        nsamples = 10
        entries = []
        for sel in selectors:
            for meas in measurements:
                secs = sorted(3600*24*30*randn(nsamples)+datetime2seconds(datetime.datetime(2014,1,2,3,4,5,6)))
                vals = randn(nsamples)
                for s,v in zip (secs,vals):
                    entries.append( (sel,meas,datetime.datetime.fromtimestamp(s),v) )

        scaleddata = {}
        for en in entries:
            selector = en[0]
            measurement_name = en[1]
            date = en[2]
            value = en[3]
            if not selector in scaleddata:
                scaleddata[selector]={}
            if not measurement_name in scaleddata[selector]:
                scaleddata[selector][measurement_name]=[]
            scaleddata[selector][measurement_name].append( (datetime2seconds(date),value,2.*value) )

        errors = {}
        errors['meas_a'] = 3
        errors['meas_b'] = 1
        errors['meas_a'] = 5

        total_errors = 0
        for key,val in errors.iteritems():
            total_errors += val
        print '     : %d problems:'%(total_errors),
        for key,val in errors.iteritems():
            print key,':',val,';',
        print ''

        reporter = Reporter(author='Me')
        reporter.addtitlepage()

        #print 'Overal: found %d (%d bad, %d retakes) of %d expected tests (net %.1f%% done) showing %d errors' \
        #      %( overal_total,overal_np,overal_retakes,overal_expected,100.*(overal_total-overal_np-overal_retakes)/overal_expected,overal_errors )

        for sel in sorted(scaleddata.keys()):
            section_title   = reporter.header(machine+': '+sel)
            section_explain1 = reporter.paragraph("Periode: 2014-01-01 t/m 2014-06-30")
            section_explain2 = reporter.paragraph("Verwachte #testen: 52")
            section_explain3 = reporter.paragraph("Uitgevoerde #testen: 95 (0 onbruikbaar, 2 herhalingen)")
            section_explain4 = reporter.paragraph("Netto score: 178.8% uitgevoerd")

            section_tabletext = reporter.paragraph("In totaal zijn 9 problemen gevonden.")
            table_data = []
            table_data.append(['meting','#problemen'])
            for key,val in sorted(errors.iteritems()):
                table_data.append( [key,val] )
            section_table = reporter.table( table_data )

            graph_data = []
            for name in sorted(scaleddata[sel].keys()):
                dadata = scaleddata[sel][name]
                xyarr = [ (da[0],da[1]) for da in dadata ]
                graph_data.append( (name,xyarr) )
            section_chart = reporter.dategraph(graph_data)

            section = [section_title, section_explain1,section_explain2,section_explain3,section_explain4,section_tabletext, section_table, section_chart]
            reporter.addsection(section)


        # colortest
        section_title = reporter.header("ColorTest")
        section_table = reporter.colorTest()
        reporter.addsection([section_title,section_table])

        # stackedChart test
        for sel in sorted(scaleddata.keys()):
            section_title   = reporter.header(machine+': '+sel)
            section_explain1 = reporter.paragraph("Periode: 2014-01-01 t/m 2014-06-30")
            section_explain2 = reporter.paragraph("Verwachte #testen: 52")
            section_explain3 = reporter.paragraph("Uitgevoerde #testen: 95 (0 onbruikbaar, 2 herhalingen)")
            section_explain4 = reporter.paragraph("Netto score: 178.8% uitgevoerd")

            section_tabletext = reporter.paragraph("In totaal zijn 9 problemen gevonden.")
            table_data = []
            table_data.append(['meting','#problemen'])
            for key,val in sorted(errors.iteritems()):
                table_data.append( [key,val] )
            section_table = reporter.table( table_data )

            graph_data = []
            for name in sorted(scaleddata[sel].keys()):
                dadata = scaleddata[sel][name]
                xyarr = [ (da[0],da[1]) for da in dadata ]
                graph_data.append( (name,xyarr) )
            section_chart = reporter.stackedGraphs(graph_data)

            section = [section_title, section_explain1,section_explain2,section_explain3,section_explain4,section_tabletext, section_table, section_chart]
            reporter.addsection(section)


        reportfile = 'gfe.pdf'
        reporter.render(reportfile)
        self.statusLabel.setText("Report written to %s." %reportfile)

    def reportMIRBucky(self):
        """
        Idee: Voor iedere bucky de volgende waardes dumpen:
            Exposure_ms, DOP, MTF_Area5, Sensitivity
        """
        selectors = ['CR F10','CR F11a', 'CR F11b', 'CR Trauma1', 'CR Trauma2']
        subset = [
            'ExposureTime (ms)', 'ImageAreaDoseProduct', 'Sensitivity', 'AreaMTF5' ]
        subset_table = []
        subset_wall = []
        for ss in subset:
            subset_table.append(ss+'_Table')
            subset_wall.append(ss+'_Wall')

        self.runmode = self.stModeReport
        # clear GUI
        self.clearLayout(self.layout)
        self.status_id_text = {}

        # table showing results
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)
        self.qctable.horizontalHeader().sectionClicked.connect(self.qctable.sortByColumn)
        self.qcgraph = None


        # make list of selectable selectors
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)
        self.selector_id_text = {}
        selected = []
        cur.execute("select * from %s" % (self.table_selector))
        rows = cur.fetchall()
        minid = None
        for row in rows:
            selected.append(row['name'])
            if minid is None:
                minid = row['pk']
            else:
                minid = min (minid,row['pk'])
            self.selector_id_text[row['name']] = row['pk']
            self.selector_id_text[row['pk']] = row['name']
        zerotext = '[Selector]'
        self.selector_id_text[zerotext] = minid-1
        self.selector_id_text[minid-1] = zerotext

        ### close db
        self.closedb(con)
        """
        Workflow:
                2. Find processes of correct modality, stationname, selector and in given reporting period
                3. Determine expected number of processes
        """
        ## 2. Find processes of correct modality, stationname, selector and in given reporting period

        for sel in selectors:
            if sel in self.selector_id_text:
                print "SEL",sel
                # select given selector
                self.selectedSelector=self.selector_id_text[sel]
                # select correct time period
                self.selectedReportPeriod = 'Previous year'
                self.determineReportPeriod()
                # clear irrelevant search fields
                self.selectedModality ='[Modality]'
                self.selectedStation = '[StationName]'
                self.selectedStatus = None

                self.onActivatedProcessorChanged()
                for lab,sss in zip(['Table','Wall'],[subset_table,subset_wall]):
                    self.getSubSelection(self.reportentries,sss)

                    fname = sel+'_'+lab+'.tsv'
                    copyline = 'selector\tdate'
                    for ss in sorted(sss):
                        copyline += '\t%s' % ss
                    copyline += '\n'
                    with open(fname, 'w') as myfile:
                        myfile.write(copyline)
                        for en in self.reportentries.values():
                            if en.analysisresult:
                                copyline = '%s\t%s' %(en.selector,str(en.date))
                                for key,val in sorted(en.analysisresult.iteritems()):
                                    copyline += '\t%s' % str(val)
                                copyline += '\n'
                                myfile.write(copyline)

    def getSubSelection(self,entries,subselection):
        """
        Check all results of given entries, identified by gewensteprocessen_id.
        If any of the tables resultaten_char, resultaten_floating is outside specs, put in report
        """
        result = {}
        ###1 open db iqc
        con,cur,dbversion = self.connectdb(host=self.host)

        for en in entries.values():
            if en.status == 'Error':
                en.analysisresult = self.stResultError
                continue
            en.analysisresult = {}
            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (self.table_char,en.gpid))
            rows = cur.fetchall()
            for row in rows:
                if row['omschrijving'] in subselection:
                    en.analysisresult[row['omschrijving']] = row['waarde']

            cur.execute("SELECT * FROM %s WHERE gewenste_processen_fk=%d" % (self.table_float,en.gpid))
            rows = cur.fetchall()
            for row in rows:
                if row['omschrijving'] in subselection:
                    en.analysisresult[row['omschrijving']] = row['waarde']

        ### close db
        self.closedb(con)
        return
#----------------------------------------------------------------------
#----------------------------------------------------------------------
class DateAxis(pg.AxisItem):
    def tickStrings(self, values, scale, spacing):
        strns = []
        try:
            rng = max(values)-min(values)
        except:
            return strns

        string = '%d-%b-%Y'
        label1 = '%d-%b,%Y -'
        label2 = ' %d-%b,%Y'

        for x in values:
            try:
                strns.append(time.strftime(string, time.localtime(x)))
            except ValueError:  ## Windows can't handle dates before 1970
                strns.append('')
        try:
            label = time.strftime(label1, time.localtime(min(values)))+time.strftime(label2, time.localtime(max(values)))
        except ValueError:
            label = ''
        #self.setLabel(text=label)
        return strns

def process_args(argv):
    parser = argparse.ArgumentParser(description='PyQt4 argstest',
                                   add_help=False)

    # reimplement help (which we disabled above) so that -help works rather
    # than --help; done to be consistent with the style of args Qt wants
    parser.add_argument("-h", "-help", action='help',
                      help="show this help message and exit")

    # AS: now my own
    parser.add_argument("-f", "--file", help="No gui, just run for given file")
    return parser.parse_args(argv[1:])

class ignore(argparse.Action):
  # we create an action that does nothing, so the Qt args do nothing
  def __call__(self, parser, namespace, values, option_string=None):
    pass

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv) # QApplication eats argv in constructor

    # We can get a QStringList out of QApplication of those arguments it
    # didn't decide were reserved by Qt.
    argv2 = app.arguments()

    # now we need to turn them back into something that optparse/argparse
    # can understand, since a QStringList is not what it wants
    argv3 = []
    for i in argv2:
        argv3.append(str(i))

    # now we can pass this to optparse/argparse
    args=process_args(argv3)

    form = dbTool()
    if(args.file != None):
        form.batchFile(args.file)
    else:
        form.show()
        app.exec_()

