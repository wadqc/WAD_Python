#!/usr/bin/python
""" Changelog:
20140521: initial version
"""
import argparse
from collections import OrderedDict
from pyqtgraph.Qt import QtGui,QtCore
import MySQLdb as mdb
import sys
from ConfigParser import SafeConfigParser

class dbTool(QtGui.QMainWindow):
    qcversion = 20140521

    host = None
    user = None
    pswd = None

    iqcdb = 'iqc'
    table_gewenste_processen = 'gewenste_processen'
    table_studies = 'study'
    table_series = 'series'
    table_instances = 'instance'
    table_selector = 'selector'
    table_status_omschrijving = 'status_omschrijving'
    status_reset = 0

    qctable = None
    qctext = None
    selectedModality = None
    selectedStation = None
    selectedSelector = None
    selectedStatus = None

    status_id_text = {}
    selector_id_text = {}

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
        print "ERROR: configfile %s corrupt. Should look like:\n[iqc]\nhost = localhost\n" \
              "user = root\npswd = ***\n\nwith valid values for host and user.\n" \
              "if pswd is blank then a dialog will ask for a valid password" % ininame
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
        if self.pswd is '':
            self.pswd = self.pswdFromDialog()

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

        self.setGeometry(1, 1, 512, 512)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        exitAction = QtGui.QAction("&Quit", self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(QtGui.qApp.quit)
        fileMenu.addAction(exitAction)

        processorMenu = menubar.addMenu("&Processor")
        procresetAction = QtGui.QAction("Reset processor status for selection", self)
        procresetAction.setShortcut('Ctrl+P')
        procresetAction.triggered.connect(self.processorReset)
        processorMenu.addAction(procresetAction)

    def clearLayout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clearLayout(item.layout())

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
        except mdb.Error, e:
            print "[connectdb] Error %d: %s" % (e.args[0],e.args[1])
            sys.exit(1)
        return con,cur

    def closedb(self,con):
        if con:
            con.close()

    def processorShowSelected(self,selected):
        ###1 open db iqc
        con,cur = self.connectdb(host=self.host)

        selrows = []
        headers = []
        for pkid in selected:
            cur.execute("select * from %s where pk=%d" % (self.table_gewenste_processen,pkid))
            rows = cur.fetchall()
            selrows.append(rows[0])  # should be only 1
            if len(headers) == 0:
                headers = rows[0].keys()

        self.qctable.setRowCount(len(selrows))
        if len(selected) > 0:
            self.qctable.setColumnCount(len(selrows[0]))
            rowid = 0
            self.qctable.setHorizontalHeaderLabels(headers)
            for row in selrows:
                colid = 0
                for hd in headers:
                    if hd == 'status':
                        valstr = self.status_id_text[row[hd]]
                    elif hd == 'selector_fk':
                        valstr = self.selector_id_text[row[hd]]
                    else:
                        valstr = str(row[hd])
                    self.qctable.setItem(rowid, colid, QtGui.QTableWidgetItem(valstr))
                    colid += 1
                rowid += 1
        ### close db
        self.closedb(con)

        self.qctext.setText("Selected %d processes." %len(selrows))

    def onActivatedProcessorSelector(self,text):
        print "[onActivatedProcessorSelector]",text,self.selector_id_text[str(text)]
        self.selectedSelector = self.selector_id_text[str(text)]

        if self.qctable is None or self.qctext is None:
            return
        self.onActivatedProcessorChanged()

    def onActivatedProcessorStatus(self,text):
        print "[onActivatedProcessorStatus]",str(text),self.status_id_text[str(text)]
        self.selectedStatus = self.status_id_text[str(text)]

        if self.qctable is None or self.qctext is None:
            return
        self.onActivatedProcessorChanged()

    def onActivatedProcessorStation(self,text):
        print "[onActivatedProcessorStation]",text
        station = str(text)
        self.selectedStation = station
        if self.qctable is None or self.qctext is None:
            return

        self.onActivatedProcessorChanged()

    def onActivatedProcessorModality(self,text):
        print "[onActivatedProcessorModality]",text
        mod = str(text)
        self.selectedModality = mod
        if self.qctable is None or self.qctext is None:
            return

        self.onActivatedProcessorChanged()

    def onActivatedProcessorChanged(self):
        ###1 open db iqc
        con,cur = self.connectdb(host=self.host)

        # 1-2: build list with (gp_is.series ids)
        cur.execute("SELECT * FROM %s" % self.table_gewenste_processen)
        rows = cur.fetchall()

        print self.selectedSelector,self.selectedStatus,self.selectedModality,self.selectedStation
        gp_series_fks = []
        for row in rows:
            if self.selectedSelector==self.selector_id_text['[Selector]'] or self.selectedSelector == row['selector_fk']:
                if self.selectedStatus==self.status_id_text['[Status]'] or self.selectedStatus == row['status']:
                    gpid = row['pk']
                    sfk = []
                    if row['series_fk']:
                       sfk.append(row['series_fk'])
                    elif row['study_fk']:
                        cur.execute("SELECT pk FROM %s WHERE study_fk=%d" % (self.table_series,row['study_fk']))
                        rows = cur.fetchall()
                        for row in rows:
                            sfk.append(row['pk'])
                    elif row['instance_fk']:
                        cur.execute("SELECT series_fk FROM %s WHERE pk=%d" % (self.table_instances,row['instance_fk']))
                        rows = cur.fetchall()
                        sfk.append(rows[0]['series_fk']) # can be only one result

                    if len(sfk) >0:
                        gp_series_fks.append([gpid,sfk])

        # 3: check modality and stationname in all series per gp_id; if match, add to reset_list
        reset_list = []
        for gpid,skfs in gp_series_fks:
            for skf in skfs:
                cur.execute("SELECT * FROM %s WHERE pk=%d" % (self.table_series,skf))
                rows = cur.fetchall()
                if self.selectedModality=="[Modality]" or self.selectedModality in rows[0]['modality']: # can be only one result
                    if self.selectedStation=="[StationName]" or self.selectedStation in rows[0]['station_name']: # can be only one result
                        reset_list.append(gpid)
                        break
        ### close db
        self.closedb(con)


        print "selected pks:",reset_list
        self.processorShowSelected(reset_list)


    def resetListedProcessor(self):
        if self.qctable is None:
            return

        if self.qctable.rowCount() == 0:
            return

        if self.qctable.columnCount() == 0:
            return

        colid = -1
        for i in range(self.qctable.columnCount()):
            if self.qctable.horizontalHeaderItem(i).text() == 'pk':
                colid = i
                break

        if colid <0:
            print "ERROR: could not find column pk"
            return

        reset_list = []
        for rowid in range(self.qctable.rowCount()):
            reset_list.append( int(self.qctable.item(rowid,colid).text()) )

        if len(reset_list) == 0:
            return

        ###1 open db iqc
        con,cur = self.connectdb(host=self.host)

        # 4: for all gp_id in reset_list, set status to "reset"
        nrows = 0
        for gpid in reset_list:
            cur.execute("UPDATE %s set status=%d where %s.pk=%d" % (self.table_gewenste_processen,self.status_reset,self.table_gewenste_processen,gpid))
            nrows += cur.rowcount

        self.qctext.setText("Did reset of %d processes" %nrows)
        con.commit()

        ### close db
        self.closedb(con)

    def processorReset(self):
        # clear GUI
        self.clearLayout(self.layout)

        # table showing selection
        self.qctable = QtGui.QTableWidget()
        self.qctable.setEditTriggers(QtGui.QTableWidget.NoEditTriggers)

        # textfield showing number of selected processes
        self.qctext = QtGui.QLabel()

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
        con,cur = self.connectdb(host=self.host)

        ### fill ComboBoxes
        selected = []
        cur.execute("select modality from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            selected.append(row['modality'])
        names = OrderedDict.fromkeys(selected).keys()
        for ho in names:
            cmModality.addItem(ho)
        cmModality.activated[str].connect(self.onActivatedProcessorModality)

        selected = []
        cur.execute("select station_name from %s" % (self.table_series))
        rows = cur.fetchall()
        for row in rows:
            selected.append(row['station_name'])
        names = OrderedDict.fromkeys(selected).keys()
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
        names = OrderedDict.fromkeys(selected).keys()
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
        for ho in names:
            cmStatus.addItem(ho)
        zerotext = str(cmStatus.itemText(0))
        self.status_id_text[zerotext] = minid-1
        self.status_id_text[minid-1] = zerotext
        self.selectedStatus = self.status_id_text[zerotext]
        cmStatus.activated[str].connect(self.onActivatedProcessorStatus)

        ### select all processes
        selected = []
        cur.execute("select pk from %s" % (self.table_gewenste_processen))
        rows = cur.fetchall()
        for row in rows:
            selected.append(row['pk'])

        ### close db
        self.closedb(con)

        self.processorShowSelected(selected)
        ### build GUI
        self.layout.addWidget(cmModality,0,0)
        self.layout.addWidget(cmStations,0,1)
        self.layout.addWidget(cmSelectors,1,0)
        self.layout.addWidget(cmStatus,1,1)

        self.layout.addWidget(self.qctable,2,0,1,2)
        self.layout.addWidget(self.qctext,3,0)
        self.layout.addWidget(btReset,3,1)
        self.qw.setLayout(self.layout)

#----------------------------------------------------------------------
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

