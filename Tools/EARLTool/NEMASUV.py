
## python suv.py -f ~/KLIFIO/ExampleData/NM/QCData/EARL/WB_CTACHEMIBODY/ -ord zxy -s 20 -A 70.77 -Ta 1710 -Ar 0.85 -Tr 1603 -V 9287



import numpy as np
import pylab as plt
from scipy import ndimage
from operator import itemgetter
from numpy import ma
import copy
import dicom
#import pydicom_series
import sys
import argparse
import os
import time
import sys
import NEMAIQ as niq
import NEMAlib as niqlib
from NEMAIQ import simple_nema_phantom
from NEMAIQ import sphere
from datetime import datetime
import time





class uniform_phantom(simple_nema_phantom):
    def  __init__(self,data,order=['x','y','z'],phantomdim=[20,10],phantomvolume=9500,dose=75,dosecaltime='000000',remdose=0,remdosecaltime='000000'):


        print order
        self.xind = order.index('x')
        self.yind = order.index('y')
        self.zind = order.index('z')
        
        axes = [self.xind,self.yind,self.zind]
        print data
        pixelsize,array,info = niq.loaddata(data)

        print info
        
        self.pixel_size = pixelsize
        print "pixelsize: ", pixelsize

        self.phantomvolume = float(phantomvolume)
        self.dose = float(dose)
        self.dosecaltime = str(dosecaltime)
        self.remdose = float(remdose)
        self.remdosecaltime = str(remdosecaltime)

        
        self.rows = info.info["0028","0010"].value
        self.cols = info.info["0028","0011"].value

        self.vendor = info.info["0008","0070"].value
        
        if self.vendor == 'Philips Medical Systems':
            self.suvfactor = info.info["7053","1009"].value
        else:
            self.suvfactor = 1.0

        self.weight = info.info["0010","1030"].value

        print "Acquisition date",info.info["0008","0022"].value #Date
        acqtime = info.info["0008","0032"].value.split('.')[0] #Acq. Time and we remove the microseconds

        
        
        print "Acquisition time:  ", acqtime
        print "Activity calibrated: %s at %s "%(self.dose,self.dosecaltime)
        print "Remainder activity: %s at %s "%(self.remdose,self.remdosecaltime)

        self.fmt = '%H%M%S'
        self.fmt2 = '%H%M'

        

        t1 = datetime.strptime(acqtime, self.fmt) #acquisition time dcmhdr
        t2 = datetime.strptime(self.dosecaltime, self.fmt) #act. cal time
        t3 = datetime.strptime(self.remdosecaltime, self.fmt)#remaining activity cal time

        # Activity concentration at acquisitiontime:


        deltat0 = float((t1 - t2).seconds/60.) #minutes
        deltatr = float((t2 - t3).seconds/60.) #minutes


        self.currconc = niqlib.calc_act_con(self.phantomvolume,self.dose,deltat0,self.remdose,deltatr)

        print "Net activity at acquisition : %s "%self.currconc
        


#float(indata.info["0028","0030"].value[0])

        self.bgradius = 4

        



        self.data = np.array(array)

        print type(self.data)
        print type(self.data[0])

        tmpshape  = np.shape(self.data)

        print tmpshape
        
        self.ddimz = tmpshape[self.zind]
        self.ddimy = tmpshape[self.yind]
        self.ddimx = tmpshape[self.xind]

        self.maxdata = self.data.max()
        self.volume = np.zeros(np.shape(data))

#        self.mode = mode
        self.z = self.z_init()


        self.thrdata = copy.deepcopy(self.data)

        lowind = (self.thrdata < np.max(self.thrdata[self.z,:,:])*0.2)
        self.thrdata[lowind] = 0.0
        highind = (self.thrdata >= np.max(self.thrdata[self.z,:,:])*0.2)
        self.thrdata[highind] = 10

        self.bgdata = copy.deepcopy(self.data)

        self.roughcom = ndimage.measurements.center_of_mass(self.data[self.z,:,:])

        print "rough com location:", self.roughcom

        self.cor = self.roughcom
        self.bgspheres = {}

        self.set_bg_spheres()






    def guess_point(self,profile):
        return simple_nema_phantom.guess_point(self,profile)

    def set_bg_spheres(self):
        return simple_nema_phantom.set_bg_spheres(self,nrangles=25,radlength=10)

    def get_bg_stats(self,data=None):
        return simple_nema_phantom.get_bg_stats(self,data=None)

    def check_overlap(self,newsphere):
        a = newsphere.volume
        for key in self.bgspheres.keys():
            b = self.bgspheres[key].volume
            if niq.compare_coordsets(a,b) == True:
                return None
        return newsphere

    def plot_bgspheres(self):
        return simple_nema_phantom.plot_bgspheres(self)

    def plot_dataspheres(self,sliceid,data=None,showfig=False):
        if data == None:
            data = self.data

        tmpdat = data[sliceid,:,:]
        tmpbg = self.plot_bgspheres()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.imshow(tmpdat)
        ax1.imshow(tmpbg,alpha=0.5)

        if showfig==True:
            plt.show()
        return fig



    def plot_img_hist(self, sid, bins=256,showfig=False):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        tmpdat = self.data[sid,:,:]
        
        plt.hist(tmpdat.flatten(), bins, range=(tmpdat.min(),tmpdat.max()), fc='k', ec='k')
        if showfig==True:
            plt.show()
        return fig


    def find_border(self,line):
        for i in range(len(line)):
            if line[i] < max(line)/2.0:
                return i


    def slice_uniformity(self,sliceid,data=None,showfig=False):

        if data == None:
            data = self.data

        thrav = np.average(data[sliceid,:,:])
        
        
        tmpdat = data[sliceid,:,:]
        
        tmpmaskdat = np.zeros(np.shape(tmpdat))
        tmpmaskdat[(data[sliceid,:,:]>thrav*.50)] = 1

        hline = []
        vline = []
        
        for i in range(self.ddimx):
            try:
                hline.append(tmpdat[self.cor[1],self.cor[0]+i])
                vline.append(tmpdat[self.cor[1]+i,self.cor[0]])
            except:
                pass
                

        tmpradius = 0.8*(min(self.find_border(vline),self.find_border(hline)))
        
        
        x_grid, y_grid = np.meshgrid(np.arange(self.ddimx), np.arange(self.ddimy))
        disk = ((x_grid - self.cor[0])**2 + (y_grid - self.cor[1])**2) <= tmpradius**2  # array of booleans with the disk shape

        points_in_disk = tmpdat[disk]
        points_in_circle2 = np.ma.masked_array(tmpdat, ~disk)

        
        

        fig = plt.figure()
        plt.imshow(points_in_circle2)
        ax1 = fig.add_subplot(111)

        if showfig==True:
            plt.show()
        return fig, np.average(points_in_disk), np.std(points_in_disk)



    def z_init(self):
        return simple_nema_phantom.z_init(self)        

    def z_slice(self):
        return simple_nema_phantom.z_slice(self)        

    def z_profile(self,data):
        return simple_nema_phantom.z_profile(self)        
    

    def visualize_slices(self, sliceid):
        
        tmpdat = self.data[sliceid,:,:]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.imshow(tmpdat)

        return fig

        #plt.show()


    def plot_zprofile(self,showfig=False):
        tmpdat = self.z_slice()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(tmpdat)
        if showfig==True:
            plt.show()
        return fig


def update_wad(filename, mode, axorder, sporder, **kwargs):

    output = {}


    #Float results

    volume = kwargs['volume']
    output['volume'] = volume
    output['length'] = kwargs['dimensions'][0]
    output['radius'] = kwargs['dimensions'][1]
    output['activity'] = kwargs['bgac']

    output['caltime'] = kwargs['bgactime']
    output['remaindoer'] = kwargs['bgremac']
    output['remtime'] = kwargs['bgremactime']

    output['deltat0'] = 999999
    output['deltatr'] = 999999
    output['scanto'] = 999999

    output['maxsuv'] =999999
    output['avgsuv'] =999999
    output['minsuv'] =999999


    
    th = niqlib.thf18
    fmt = '%H%M%S' #assumed format for all times

    # convert all 'string' times to datetime format so we can calculate something

    scantime = datetime.strptime(kwargs['scantime'],fmt)
    ctimespac = datetime.strptime(kwargs['ctimespac'],fmt) 
    bgactime = datetime.strptime(kwargs['bgactime'],fmt) 
    spremactime = datetime.strptime(kwargs['spremactime'],fmt) 
    bgremactime = datetime.strptime(kwargs['bgremactime'],fmt) 

    actbgsctime = niqlib.decay(float(kwargs['bgac']),float((scantime - bgactime).seconds/60.),th) - niqlib.decay(float(kwargs['bgremac']),float((scantime - bgremactime).seconds/60.),th) 

    output['Act BG (scantime)'] = actbgsctime

    output['Scantime'] = scantime.time()



    print "Writing results to wadserver" 


    print "1 - loading data"
    phantom = uniform_phantom(data=filename, order=axorder,phantomdim=kwargs['dimensions'],phantomvolume=kwargs['volume'],dose=kwargs['bgac'],dosecaltime=kwargs['bgactime'],remdose=kwargs['bgremac'],remdosecaltime=kwargs['bgremactime'])

    print phantom.z
    output['z_slice'] = phantom.z

    zslicefig = plt.figure()
    ax1 = zslicefig.add_subplot(111)

    plt.xtitle = "Z slice"
    plt.ytitle = "Counts in slice"
    plt.plot(phantom.z_slice())

    #plt.show()

    output['EARL zprofile'] = zslicefig
    
    print "2 - info from header"
    """We need to calculate the activity concentration in Bq/ml from the
volume of the phantom and the administered activity.
     """

    # suv =  c(t) / (A(t) / bodyweight)
    # 
    # In the case of suv phantom:
    # suv = c(t) / (A(t) / 9.293), where A in MBq and c in kBq/cc
    # 
    # dicom[7053,1000] converts to SUV

    # For Philips specific:
    # There are 2 ways to calculate the suv
    # 1 - multiply pixel value with tag 7053,1000 and make sure the correct weight was used (or correct)
    # 2 - multiply pix. val. with tag 7053, 1009 (rescale slope) which converts to bq/cc and correct for decay
    

    suv = lambda x:1000*x*float(phantom.suvfactor)/phantom.currconc #option 2

 
    print "3 - show result: "

    sid  = phantom.z
    output['EARL histogram']=phantom.plot_img_hist(sid, bins=256,showfig=False)





    for k in [sid - 6, sid - 4, sid, sid + 4, sid + 6]:
        tmpfig, tmpunif,tmpstd = phantom.slice_uniformity(k,showfig=False)
        print "uniformity: ", k, tmpunif, tmpstd, actbgsctime, suv(tmpunif)
        output['SUV slice %s'%str(k - sid)] =  suv(tmpunif)
        output['COV slice %s'%str(k - sid)] =  tmpstd/tmpunif


    output['EARL halfslice'] =  phantom.plot_dataspheres(sid) #phantom.visualize_slices(sid)

    

    print "4 - sphere statistics: "
    

    bgresults  = phantom.get_bg_stats()

    meansuvlist = []
    maxsuvlist = []

    for key in bgresults.keys():
        print "BG Sphere: ", key 
        tmpresult = bgresults[key]
        meansuvlist.append([key,tmpresult['avg'],tmpresult['avg']/phantom.currconc])
        maxsuvlist.append([key,tmpresult['max'],tmpresult['max']/phantom.currconc])
        for subkey in tmpresult.keys():
            print '\t', subkey, tmpresult[subkey]


    print "SUVmean: ",meansuvlist
    print "SUVmax: ",maxsuvlist

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(zip(*meansuvlist)[2])
    #plt.show()
    


    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(zip(*maxsuvlist)[2])
    output['EARL suvcurve'] = fig

    #plt.show()
    

    



    output['procstat'] = 'processed'
    output['status'] = 'processed'



    return output
    


def main():
    parser = argparse.ArgumentParser(description='Commandline arguments to run reconstructor')
  
    parser.add_argument('-f',metavar='file',nargs='?',type=str)
    parser.add_argument('-ord',metavar='order',nargs='?',type=list)
    parser.add_argument('-s',metavar='slice',nargs='?',type=int)

    parser.add_argument('-A',metavar='Activity (MBq)',nargs='?',type=float)
    parser.add_argument('-Ta',metavar='Calibration time',nargs='?',type=str)
    parser.add_argument('-Ar',metavar='Remaining activity',nargs='?',type=float)
    parser.add_argument('-Tr',metavar='Calib time',nargs='?',type=str)

    parser.add_argument('-V',metavar='Volume',nargs='?',type=float)
    parser.add_argument('-len',metavar='Phantom length',nargs='?',type=float)
    parser.add_argument('-rad',metavar='Phantom radius',nargs='?',type=float)

    args = vars(parser.parse_args())
    filename = args['f'] #dicom directory
    order = args['ord']
    sid = args['s']

    activity = args['A']
    caltime = args['Ta']
    remainder = args['Ar']
    remtime = args['Tr']

    length = args['len']
    rad = args['rad']
    volume = args['V']

    dimensions = [length,rad]
    print dimensions

    print "1 - loading data"
    phantom = uniform_phantom(data=filename, order=order,phantomdim=dimensions,phantomvolume=volume,dose=activity,dosecaltime=caltime,remdose=remainder,remdosecaltime=remtime)

    print phantom.z
    
    print "2 - info from header"
    """We need to calculate the activity concentration in Bq/ml from the
volume of the phantom and the administered activity.
 
 """

    print "z_slice bollen: "

    print "2 - registering data"
 
    print "3 - show result: "
    phantom.plot_zprofile()
    phantom.visualize_slices(sid)
    phantom.plot_dataspheres(sid)


    print "4 - sphere statistics: "

    bgresults  = phantom.get_bg_stats()

    meansuvlist = []
    maxsuvlist = []

    for key in bgresults.keys():
        print "BG Sphere: ", key 
        tmpresult = bgresults[key]
        meansuvlist.append([key,tmpresult['avg'],tmpresult['avg']/phantom.currconc])
        maxsuvlist.append([key,tmpresult['max'],tmpresult['max']/phantom.currconc])
        for subkey in tmpresult.keys():
            print '\t', subkey, tmpresult[subkey]


    print "SUVmean: ",meansuvlist
    print "SUVmax: ",maxsuvlist

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(zip(*meansuvlist)[2])
    plt.show()
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(zip(*maxsuvlist)[2])
    plt.show()
    


if __name__ == "__main__":
    sys.exit(main())
