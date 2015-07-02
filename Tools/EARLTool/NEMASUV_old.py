
## python suv.py -f ~/KLIFIO/ExampleData/NM/QCData/EARL/WB_CTACHEMIBODY/ -ord zxy -s 20 -A 70.77 -Ta 1710 -Ar 0.85 -Tr 1603 -V 9287



import numpy as np
import pylab as plt
from scipy import ndimage
from operator import itemgetter
from numpy import ma
import copy
import dicom
import pydicom_series
import sys
import argparse
import os
import time
import sys
import NEMAIQ as niq
from NEMAIQ import simple_nema_phantom
from NEMAIQ import sphere
from datetime import datetime
import time

thf18 = 119 #half life time in minutes



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
        
        print "info - date",info.info["0008","0022"].value #Date
        acqtime = info.info["0008","0032"].value.split('.')[0] #Acq. Time and we remove the microseconds

        
        
        print "acquisition time:  ", acqtime
        print "acquisition time:  ", self.dosecaltime
        print "acquisition time:  ", self.remdosecaltime

        self.fmt = '%H%M%S'
        self.fmt2 = '%H%M'

        

        t1 = datetime.strptime(acqtime, self.fmt) #acquisition time dcmhdr
        t2 = datetime.strptime(self.dosecaltime, self.fmt) #act. cal time
        t3 = datetime.strptime(self.remdosecaltime, self.fmt)#remaining activity cal time

        # Activity concentration at acquisitiontime:


        deltat0 = float((t1 - t2).seconds/60.) #minutes
        deltatr = float((t2 - t3).seconds/60.) #minutes


        self.currconc = self.calc_act_con(self.phantomvolume,self.dose,deltat0,self.remdose,deltatr)

        print self.currconc
        


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


    def calc_act_con(self,volume,activity,deltat0, remainder,deltatr):
        
        print "volume: (ml)", volume
        print "activity: (MBq)", activity
        print "deltat0: ", deltat0
        print "remainder: ", remainder
        print "deltatr: ", deltatr 
 
        Ar = self.decay(remainder,deltatr,thf18) 
        A0 = self.decay((activity - Ar),deltat0,thf18) 

        print "Decay corrected Ar, A: ", Ar,A0

        accon = 1e6*(A0 - Ar)/volume
        print "Activity concentration: Bq/ml",accon
        return accon


    def decay(self,A0,t,th):
        return A0*np.power(0.5,(t/th))


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

    def plot_dataspheres(self,sliceid,data=None):
        if data == None:
            data = self.data

        tmpdat = data[sliceid,:,:]
        tmpbg = self.plot_bgspheres()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.imshow(tmpdat)
        ax1.imshow(tmpbg,alpha=0.5)

        plt.show()


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

        plt.show()


    def plot_zprofile(self):
        tmpdat = self.z_slice()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(tmpdat)
        plt.show()



def update_wad(filename, mode, axorder, sporder,bgac, dimensions, volume, caltime, remainder, remtime ):

    output = {}


    #Float results

    output['volume'] = volume
    output['length'] = dimensions[0]
    output['radius'] = dimensions[1]
    output['activity'] = bgac
    output['deltat0'] = 999999
    output['caltime'] = caltime
    output['remaindoer'] = remainder
    output['remtime'] = remtime
    output['maxsuv'] =999999
    output['avgsuv'] =999999
    output['minsuv'] =999999
    output['procstat'] = 'processed'
    
    print "Writing results to wadserver" 


    print "1 - loading data"
    phantom = uniform_phantom(data=filename, order=axorder,phantomdim=dimensions,phantomvolume=volume,dose=bgac,dosecaltime=caltime,remdose=remainder,remdosecaltime=remtime)

    print phantom.z
    output['z_slice'] = phantom.z

    zslicefig = plt.figure()
    ax1 = zslicefig.add_subplot(111)

    plt.xtitle = "Z slice"
    plt.ytitle = "Counts in slice"
    plt.plot(phantom.z_slice())

    plt.show()

    output['EARL zprofile'] = zslicefig
    output['EARL halfslice'] = zslicefig
    output['EARL suvcurve'] = zslicefig
    
    print "2 - info from header"
    """We need to calculate the activity concentration in Bq/ml from the
volume of the phantom and the administered activity.
 
 """


    print "2 - registering data"
 
    print "3 - show result: "
    phantom.plot_zprofile()

    sid  = phantom.z
    
    phantom.visualize_slices(sid)
    phantom.plot_dataspheres(sid)


    print "4 - sphere statistics: "

    #floatlist = ['procstat','volume','length','radius','activity','deltat0','remaindoer','deltatr','scanto','maxsuv','avgsuv','minsuv']
    #imlist = ['zprofile','halfslice','suvcurve']

    
    

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
