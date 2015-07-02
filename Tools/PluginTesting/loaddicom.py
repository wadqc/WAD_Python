import os
import dicom


class _LoadDicom():
    def __init__(self, initialdir):

        
        self.result = {}
        self.series = {}
        self.UIDs = []
        
        self.initialdir = initialdir if initialdir and os.path.exists(initialdir) else None
        self.path = ''
        
        self.loadSeries()

    def loadSeries(self):
        filelist = []
        new_keys = []
        path = self.initialdir
        if path and os.path.exists(path):
            self.initialdir = path


            abslist = []

            for dirname, dirnames, filenames in os.walk(path):
                for filename in filenames:
                    abslist.append(os.path.join(dirname,filename))

#            abslist = [os.path.join(path, entry) for entry in os.listdir(path)]
            filelist = [entry for entry in abslist if os.path.isfile(entry)]
            

        for i, filename in enumerate(filelist):

            #print i, filename

            try:
                dc = dicom.read_file(filename)
            except dicom.filereader.InvalidDicomError:
                continue   
            key = self.getKey(dc)
        
        


            uid = dc.get('SOPInstanceUID')
            if uid in self.UIDs:
                ##Skip duplicates
                continue
            self.UIDs.append(uid)


            try:
                self.series[key]["dc"].append(dc)
                self.series[key]["filelst"].append((dc.InstanceNumber,filename))
 
    
                #self.series[key]["thumb"].append(photo)


            except KeyError:
                new_entry = {"path": path,
                             "file": filename,
                             "filelst":[(dc.InstanceNumber,filename)],
                             "id": dc.get('SeriesInstanceUID'),
                             "scroll": None,
                             "dc": [dc]}
                self.series[key] = new_entry
                new_keys.append(key)

            

        
        if new_keys:
            for key in new_keys:
                try:
                    self.sortDCs(key, 'ImageIndex')
                except TypeError:
                    self.sortDCs(key, 'InstanceNumber')
            self.sortKeys(new_keys)


            

    
    def sortBox(self):
        selected = None
        sels = self.box.getcurselection()
        if sels:
            selected = sels[0]            
        keys = list(self.box.get())
        self.sortKeys(keys)
        self.box.setlist(tuple(keys))
        if selected:
            self.box.setvalue(selected)
    
    def sortDCs(self, key, to):
        self.series[key]["dc"].sort(key=lambda x: int(x.get(to)))
        
    def sortKeys(self, keys):
        keys.sort(key=lambda x: self.series[x]["id"])
            
    def getKey(self, dc):


        patient = dc.get('PatientsName')
        study = dc.get('StudyDescription')
        series = dc.get('SeriesDescription')
        snum = dc.get('SeriesNumber')
        
        try:
            dt = dc.get('SeriesDate')
            tm = dc.get('SeriesTime')
            date = "%s-%s-%s"%(dt[6:], dt[4:6], dt[:4])
            time = "%s:%s:%s"%(tm[:2], tm[2:4], tm[4:6])
        
        except:
            date = ''
            time = ''
            #print ('Warning could not get series date/time')
        
        return "%s - %s - %s - %s - %s - %s"%(date, time, patient, study, series, snum)

def LoadDicom(initialdir=None):
    dc, path = [], ''
    w = _LoadDicom(initialdir)

    #for key in w.series.keys():
    #    print w.series[key]['filelst']
    
    return w.series





    #return w.series
        
if __name__ == "__main__":
    t = LoadDicom('/home/dickrsch/KLIFIO/ExampleData/NM/QCData/EARL')
    for key in t.keys():
        print t[key]["filelst"]
