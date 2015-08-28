"""
Dit script analyseert automatisch opnames gemaakt met het NemaIQ fantoom in het kader van de EARL registratie.

De output is een dicom-file met resultaten die naar de DCM4CHEE database wordt gegenereerd.

Het script parsenemaiq.py vertaalt de resultaten naar de waddb

Er is een simple_nema_phantom class die de daadwerkelijke data en een model van het nema fantoom registreert en vervolgens per bol statistische informatie, contrasten met bgrois en recovery coefficienten berekend.

Vanaf de command-line wordt de code als volgt gerund:

python NemaIQ.py -f ../NemaIQtool/Voorbeelddata/EARLIQ -m max  -ord zxy -ort acw -spa 10 -bga 1

Hierbij zijn:

m - gaat het om hotspots = max of coldspots = min
ord - de volgorde van de assen
ort - de volgorde van de bollen, tegen de klok in acw of met de klok mee cw
spa - de activiteitsconcentratie in de bollen
bga - de activiteitsconcentratie in de achtergrond

DCM4CHEE - locatie DCM4CHEE server

TO DO:

- partiele voxels meenemen om een correcter volume mee te nemen.

- testen op andere datasets

2013 DD
"""

import numpy as np
import pylab as plt
from scipy import ndimage
from operator import itemgetter
from numpy import ma
import copy


import sys
import argparse
import os
import time




from datetime import datetime
import NEMAlib as niqlib

import sys
sys.path.append("../../")
from pyWADLib import pydicom_series



def plot_rclist(rclist,order='az',showfig=False):
    RCspecsA50 = [[0.71, 0.83], [0.67, 0.79], [0.59, 0.73], [0.53, 0.68], [0.38, 0.52], [0.25, 0.35]]
    RCspecsMax = [[0.88, 1.08], [0.85, 1.05], [0.77, 1.01], [0.75, 0.94], [0.51, 0.74], [0.29, 0.46]]

    data = [(key,rclist[key]) for key in rclist.keys()]
    X,Y = zip(*data)

    if order == 'za':
        X = sorted(X,reverse=True)


    for elem in range(len(RCspecsA50)):
        l = plt.axvline(x=elem+1, ymin=RCspecsA50[elem][0], ymax= RCspecsA50[elem][1], linewidth=5, color='r')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    plt.plot( X , Y )
    plt.xtitle = "Sphere index"
    plt.ytitle = "RC coefficient"
    plt.grid(True)

    if showfig==True:
        plt.show()

    return fig


def loaddata_old(dicomdata):
    indata = pydicom_series.read_files(dicomdata,showProgress=True)[0]  


    pixelsize =  float(indata.info["0028","0030"].value[0])

    array = indata.get_pixel_array()
    return pixelsize, array, indata


def loaddata(dicomdata):
    indata = pydicom_series.read_files(dicomdata,showProgress=True)[0]  


    pixelsize =  float(indata.info["0028","0030"].value[0])
    
    rows = indata.info["0028","0010"].value
    cols = indata.info["0028","0011"].value

    newshape = (rows,cols)
    print "shape of dataslices: ", newshape
    
    bitdepth = indata.info["0028","0100"].value #bits allocated
    if bitdepth == 16:
        bv = np.int16
    elif bitdepth == 8:
        bv = np.int8    
    elif bitdepth == 32:
        bv = np.int32
    elif bitdepth == 64:
        bv = np.int64
    
    
    datlist = []

    sequence = indata.__dict__['_datasets']
    
    for i in range(len(sequence)):
       tmpkeys = sequence[i].keys()      
       datlist.append(np.reshape(np.fromstring(sequence[i][dicom.tag.Tag("7fe0","0010")].value,dtype=bv),newshape))
       

    array = np.array(datlist)
    print "loaddata complete"
    
    return pixelsize, array, indata


def compare_coordsets(a,b):
    i = 0
    while i < len(a):
        if a[i] in b:
            return True
        i+=1
    return False

class progressBar:
    def __init__(self, minValue = 0, maxValue = 100, totalWidth=75):
        """ Initializes the progress bar. """
        self.progBar = ""   # This holds the progress bar string
        self.oldprogBar = ""
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.updateAmount(0)  # Build progress bar string

    def appendAmount(self, append):
        """ Increases the current amount of the value of append and 
        updates the progress bar to new ammount. """
        self.updateAmount(self.amount + append)
    
    def updatePercentage(self, newPercentage):
		""" Updates the progress bar to the new percentage. """
		self.updateAmount((newPercentage * float(self.max)) / 100.0)

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
        values set at initialization; if it is over or under, it takes the
        min or max value as a default. """
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1), ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ' '.join([self.progBar, percentString])
    
    def draw(self):
        """ Draws the progress bar if it has changed from it's previous value.  """
        if self.progBar != self.oldprogBar:
            self.oldprogBar = self.progBar
            sys.stdout.write(self.progBar + '\r')
            sys.stdout.flush()      # force updating of screen

    def __str__(self):
        """ Returns the current progress bar. """
        return str(self.progBar)


class sphere():
    def __init__(self,datsize,spherecom):
        self.com_coord = (spherecom[0],spherecom[1])
        self.z = spherecom[2]
        self.radius = spherecom[3]
        self.rcntr=[0,0]
        self.substeps = 2
        self.volume = self.volume_init()
        self.matsize = datsize
        self.angle = 0

    def rotate_point(self,point,angle,rcntr):
 
        
        tmppoint = np.subtract([point[0],point[1]],rcntr)
        rotpoint = (np.cos(angle - self.angle)*tmppoint[0]-np.sin(angle - self.angle)*tmppoint[1], np.sin(angle-self.angle)*tmppoint[0]+np.cos(angle - self.angle)*tmppoint[1])

        self.angle=angle

        return np.add(rcntr,rotpoint)

    def rotate_volume(self,angle,rcntr):
        
        self.com_coord = self.rotate_point(self.com_coord,angle,rcntr)
        self.volume = self.volume_init()

    def translate_volume(self,vec):
        tmpvec = np.add(np.array(self.com_coord),vec)
        self.com_coord = tmpvec
        self.volume = self.volume_init()

    def distance(self,voxel):
        x,y,z = voxel

        dist2 = x**2 + y**2 + z**2

        return dist2

    def inorout(self,voxel):

        dist2 = self.distance(voxel)
        if dist2 <= self.radius**2:
            flag = 1
        else:
            flag = 0
        
        return flag

    def get_weight(self,voxel):
        
        out = []
        x,y,z = voxel
        dist2 = self.distance(voxel)
          
        if dist2 < self.radius**2:
            weight = 1.
        elif dist2 > (self.radius+1)**2:
            weight = 0.
        else:
            weight=0
            step = 1./(2**int(self.substeps))
            #print "Voxel substeps: ", self.substeps, step

            xsub = np.arange(x-.5,x+.5,step)
            ysub = np.arange(y-.5,y+.5,step)
            zsub = np.arange(z-.5,z+.5,step)
            
            subvoxels = zip(xsub,ysub,zsub)

            for elem in subvoxels:
                if self.distance(elem) <= self.radius**2:
                    weight += 1.0/(2**int(self.substeps))


        return [x+int(self.com_coord[0]),y+int(self.com_coord[1]), z+self.z,weight]


    def volume_init(self,subpixlevel=1):
        
        '''
        This function returns a weighted list of pixels that make up the volume
        of the sphere. (0,0,0) is the center of the sphere from which distances are measured.
        '''

        out = []

        surein = self.radius - 1 # alle pixels die geheel binnen de cirkel liggen (straal cirkel - verste hoek van midden van pixel
        sureout = self.radius + 1

        for i in range(-int(self.radius + 1),int(self.radius) + 1):
            for j in range(-int(self.radius + 1),int(self.radius) + 1):
                for k in range(-int(self.radius + 1),int(self.radius) + 1):

                    tmpvoxel = (i,j,k)
                    out.append(self.get_weight(tmpvoxel))

        return out

    def vol2array(self):
        outarray = np.zeros(self.matsize)
        for elem in self.volume:
            try:
                outarray[elem[2],elem[0],elem[1]] = 1
            except:
                print "warning, point not in volume"
                pass
        return outarray

    def sphere_mask(self,data):
        out = []
        for elem in self.volume:
            try:
                out.append(elem[3]*data[elem[2],elem[0],elem[1]])
            except:
                print "requested datapoint outside volume"
                pass
        return out
            


class simple_nema_phantom():
    def __init__(self,data,mode,order=['x','y','z'],orientation='cw'):
     
        print order
        self.xind = order.index('x')
        self.yind = order.index('y')
        self.zind = order.index('z')
        
        axes = [self.xind,self.yind,self.zind]

        self.spheres = {}
        self.bgspheres = {}



        print data
        pixelsize,array,info = loaddata_old(data)
        self.pixel_size = pixelsize
        print "pixelsize: ", pixelsize

        self.bgradius = 37/self.pixel_size/2.0

        self.data = array 
    

        tmpshape  = np.shape(self.data)
        self.ddimz = tmpshape[self.zind]
        self.ddimy = tmpshape[self.xind]
        self.ddimx = tmpshape[self.yind]

        

        self.maxdata = self.data.max()
        self.volume = np.zeros(np.shape(data))
    
        self.mode = mode
        self.z = self.z_init(mode=self.mode)
        


        self.thrdata = copy.deepcopy(self.data)
        lowind = (self.thrdata < np.max(self.thrdata[self.z,:,:])*0.2)
        self.thrdata[lowind] = 0.0
        highind = (self.thrdata >= np.max(self.thrdata[self.z,:,:])*0.2)
        self.thrdata[highind] = 10


        self.bgdata = copy.deepcopy(self.data)
#        bgthrp = (self.bgdata > np.max(self.data)/10.0)
#        self.bgdata[bgthrp] = 1.0
        bgthrm = (self.bgdata <= np.max(self.data)/20.0)
        self.bgdata[bgthrm] = 0.0


        self.angle = 0
        self.searchrange = 5
        self.cor = [self.ddimx/2,self.ddimy/2]

        self.roughcom = ndimage.measurements.center_of_mass(self.data[self.z,:,:])

        print "rough com location:", self.roughcom

        self.sphror = orientation
        self.sphere_com = self.nema_coords(pixelsize,0,self.cor,orientation=self.sphror)

        if orientation=='cw':
            self.maxkey = 1
        elif orientation=='acw':
            self.maxkey = 6

        for i in range(1,7):
            self.spheres[i] = sphere(tmpshape,self.sphere_com[i])



    def check_overlap(self,newsphere):
        
        a = newsphere.volume
        

        for key in self.spheres.keys():
            b = self.spheres[key].volume
            if compare_coordsets(a,b) == True:
                return None

        for key in self.bgspheres.keys():
            b = self.bgspheres[key].volume
            if compare_coordsets(a,b) == True:
                return None

        return newsphere


    def set_bg_spheres(self,nrangles=25,radlength=250/4.0):
        

        tmpshape = np.shape(self.data)

        mag = lambda x: np.sqrt(np.dot(x,x.T))

        vec = (0,radlength)
        steps = 100
        X = []
        Y = []
        bgid = 0

        for alpha in range(nrangles):
            tmpangle = alpha*2*np.pi/nrangles 
            tmpvec = (np.cos(tmpangle)*vec[0]-np.sin(tmpangle)*vec[1]/steps, np.sin(tmpangle)*vec[0]+np.cos(tmpangle)*vec[1]/steps)

            tmpprofile = []


            for i in range(steps):
                
                tmppoint = np.add(self.cor,(tmpvec[0]*i,tmpvec[1]*i))
                

                try:
                    tmpprofile.append([self.bgdata[self.z,int(tmppoint[0]),int(tmppoint[1])],int(tmppoint[1]),int(tmppoint[0])])
                except:
                    print "error!"
                    print "tmppoint", tmppoint, self.z, np.shape(self.bgdata)                    
 
            newpoint = self.guess_point(tmpprofile)
            tmpcoords = tmpprofile[newpoint]
#            print alpha,newpoint, tmpcoords
            tmpx = tmpcoords[2]
            tmpy = tmpcoords[1]
            X.append(tmpx)
            Y.append(tmpy)

            tmpcom = [tmpx,tmpy,self.z,self.bgradius]
            tmpbgsphere = sphere(tmpshape,tmpcom)
            
            
            addsphere = self.check_overlap(tmpbgsphere)
            if addsphere is not None:
                self.bgspheres[bgid] = tmpbgsphere
                bgid+=1


    def guess_point(self,profile):
        status = True
        j = len(profile)
        while status == True and j>0:
            j-=1
            if profile[j][0] > 0:
                status = False
        epsilon = 3
        radius = self.bgradius
        offset = radius + epsilon
        return int(j - offset)


       

    def z_slice(self):
        out = []
        for i in range(self.ddimz):
            if self.zind == 0:
                out.append(np.sum(self.data[i,:,:]))
            elif self.zind == 1:
                out.append(np.sum(self.data[:,i,:]))
            elif self.zind == 2:
                out.append(np.sum(self.data[:,:,i]))
        return out

    def z_init(self,mode='min'):
        
        tmpY = self.z_slice()
        tmpX = np.arange(0,len(tmpY))

        av = int(np.average(tmpX,weights=tmpY))
        width = int(len(tmpX)/4)
        
        

        if mode == 'min':
            out = tmpY.index(np.min(tmpY[av-width:av+width]))
        elif mode == 'max':
            out = tmpY.index(np.max(tmpY[av-width:av+width]))
        return out
            
    def nema_coords(self,pixel_size, rotation=0, com=[0,0],orientation='cw'):

        angle = 2.0*np.pi/6.0
        longradius = (114.4/2.0)/pixel_size #  (114.4/2.0)/pixel_size #in pixels

        vec = (0,longradius) 

        coords = {}
        bgcoords = {}

        for i in range(1,7):

            tmpangle = (i-1)*angle + rotation
            tmpvec = (np.cos(tmpangle)*vec[0]-np.sin(tmpangle)*vec[1], np.sin(tmpangle)*vec[0]+np.cos(tmpangle)*vec[1])


            if orientation == 'cw':
                coords[i] = [com[1]+tmpvec[1],com[0]+tmpvec[0],self.z,self.rad_sphere(i)/pixel_size]    #in pixels not mm
            elif orientation == 'acw':
                coords[i] = [com[1]+tmpvec[1],com[0]+tmpvec[0],self.z,self.rad_sphere(7-i)/pixel_size]    #in pixels not mm



        return coords








    def rad_sphere(self,index):
        if index == 1:
            out = 37
        elif index == 2:
            out = 28
        elif index == 3:
            out = 22
        elif index == 4:
            out = 17
        elif index == 5:
            out = 13
        elif index == 6:
            out = 10
        else:
            out = -1
        return out/2.0



    def register_phantom(self):

        if self.mode == 'min':            
            initval = 100000000

        if self.mode == 'max':
            initval = 0

        regpar = []
        regval = initval

        
        prog = progressBar(maxValue = 4*self.searchrange*self.searchrange)


        for i in range(int(self.roughcom[0])-self.searchrange, int(self.roughcom[0])+self.searchrange):
            for j in range(int(self.roughcom[1])-self.searchrange,int(self.roughcom[1])+self.searchrange):

                prog.appendAmount(1)
                prog.draw()

                

                tmpcor = (i,j)

                self.translate_spheres(tmpcor)

                tmpangle = self.find_angle()
                self.rotate_spheres(tmpangle,self.cor)

                tmpval = self.get_stats()

                summedval = 0
                for key in tmpval.keys():
                    summedval+=tmpval[key]['avg']

                if self.mode == 'min':
                   if  summedval < regval:
                       regval = summedval
                       regpars = [tmpcor,tmpangle,tmpval]

                elif self.mode == 'max':
                   if  summedval > regval:
                       regval = summedval
                       regpars = [tmpcor,tmpangle,tmpval]

                       
        self.set_bg_spheres()
                
        self.translate_spheres(regpars[0])
        self.rotate_spheres(regpars[1],self.cor)

        return regpars
        

    def data_com(self,z):
        return ndimage.measurements.center_of_mass(self.data[z,:,:])

    def rotate_spheres(self,angle,rotp):
        self.angle = angle
        
        for key in self.spheres.keys():
            self.spheres[key].rotate_volume(angle,rotp)

    def rotate_bg_spheres(self,angle,rotp):
        self.angle = angle
        
        for key in self.bgspheres.keys():
            self.bgspheres[key].rotate_volume(angle,rotp)


    def translate_spheres(self,vec):

        for key in self.spheres.keys():
            tmpcom = self.spheres[key].com_coord            
            tmpdiff = np.subtract(self.cor,vec)
            tmpdiff = -tmpdiff
            self.spheres[key].translate_volume(tmpdiff)
        self.cor = vec


    def translate_bgspheres(self,vec):

        for key in self.bgspheres.keys():
            tmpcom = self.bgspheres[key].com_coord            
            tmpdiff = np.subtract(self.cor,vec)
            tmpdiff = -tmpdiff
            self.bgspheres[key].translate_volume(tmpdiff)
        self.cor = vec


    def get_stats(self,data=None):
        if data == None:
            data = self.data
        out = {}
        for key in self.spheres.keys():
            tmpdat = self.spheres[key].sphere_mask(data)
            out[key] = {'volume':len(tmpdat),'sum':np.sum(tmpdat),'mean':np.mean(tmpdat),'std':np.std(tmpdat),'min':np.min(tmpdat),'max':np.max(tmpdat),'avg':np.average(tmpdat)}

        return out

    def get_bg_stats(self,data=None):
        if data == None:
            data = self.data

        bgout = {}
        for key in self.bgspheres.keys():
            tmpdat = self.bgspheres[key].sphere_mask(data)
            bgout[key] = {'volume':len(tmpdat),'sum':np.sum(tmpdat),'mean':np.mean(tmpdat),'std':np.std(tmpdat),'min':np.min(tmpdat),'max':np.max(tmpdat),'avg':np.average(tmpdat)}
        return bgout

    def plot_spheres(self):
        outarray = np.zeros((self.ddimx,self.ddimy))
        for i in range(1,7):
            tmpsphere = (self.spheres[i])
            for elem in tmpsphere.volume:
                try:
                    outarray[elem[0],elem[1]] = elem[3]
                except:
                    print "warning, point not in volume"
                    pass

        return outarray

    def plot_bgspheres(self):
        outarray = np.zeros((self.ddimx,self.ddimy))
        for i in range(len(self.bgspheres.keys())):
            tmpsphere = (self.bgspheres[i])
            for elem in tmpsphere.volume:
                try:
                    outarray[elem[0],elem[1]] = elem[3]
                except:
                    print "warning, point not in volume"
                    pass
        return outarray


    def plot_dataspheres(self,sliceid,data=None,showfig=False):
        if data == None:
            data = self.data

        tmpdat = data[sliceid,:,:]
        tmpsph = self.plot_spheres()
        tmpbg = self.plot_bgspheres()
        

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.imshow(tmpdat)
        ax1.imshow(tmpsph,alpha=0.4)
        ax1.imshow(tmpbg,alpha=0.4)

        if showfig == True:
            plt.show()
        return fig


    def rot_profile(self,cor=None):
        X = np.arange(0,2*np.pi,0.05) #absolute angles
        Y = []
        if cor == None:
            cor = self.cor
        for x in X:
            
            self.rotate_spheres(x,cor)
            tmpdat = self.get_stats(self.thrdata)
            
            summed_mean = np.sum([tmpdat[key]['sum'] for key in tmpdat.keys()])
            
            Y.append(summed_mean)
        return X,Y


    def find_angle(self,cor=None):
        if cor == None:
            cor = self.cor
        tmpX,tmpY=self.rot_profile(cor)
        
        if self.mode == 'min':
            out = tmpX[tmpY.index(np.min(tmpY))]
        elif self.mode == 'max':
            out = tmpX[tmpY.index(np.max(tmpY))]
        return out


def update_wad(filename, mode, axorder, sporder, **kwargs):

    #sphereac = 999, ctimespac = '120000', bgac = 999, bgactime = '120000', spremac = 999, spremactime = '120000', bgremac = 999, bgremactime = '120000',scantime='120000'



    output = {}


    #Float values

    fmt = '%H%M%S' #assumed format for all times



    output['Sphere activity'] = kwargs['sphereac']
    output['Calibration time sphere'] = kwargs['ctimespac']



    output['Background activity'] = kwargs['bgac']
    output['Background act time'] = kwargs['bgactime']



    output['Sphere rem act'] = kwargs['spremac']
    output['Sphere rem act time'] = kwargs['spremactime']



    output['Bg rem act'] = kwargs['bgremac']
    output['Bg rem act time'] = kwargs['bgremactime']




    ''' With all the activities known we can simply calculate
the activity concentration of the spheres @ scantime and the activity concentration of the background @ scantime
    '''
    th = niqlib.thf18

    # convert all 'string' times to datetime format so we can calculate something

    scantime = datetime.strptime(kwargs['scantime'],fmt)
    ctimespac = datetime.strptime(kwargs['ctimespac'],fmt) 
    bgactime = datetime.strptime(kwargs['bgactime'],fmt) 
    spremactime = datetime.strptime(kwargs['spremactime'],fmt) 
    bgremactime = datetime.strptime(kwargs['bgremactime'],fmt) 



    actspsctime = niqlib.decay(float(kwargs['sphereac']),float((scantime - ctimespac).seconds/60.),th) - niqlib.decay(float(kwargs['spremac']),float((scantime - spremactime).seconds/60.),th) 

    actbgsctime = niqlib.decay(float(kwargs['bgac']),float((scantime - bgactime).seconds/60.),th) - niqlib.decay(float(kwargs['bgremac']),float((scantime - bgremactime).seconds/60.),th) 

    output['Act Spheres (scantime)'] = actspsctime
    output['Act BG (scantime)'] = actbgsctime


    output['Scantime'] = scantime.time()

    print axorder 

    print "1 - loading data"
    iqphantom = simple_nema_phantom(data=filename,mode=mode,order=axorder,orientation=sporder)

    print "z_slice bollen: ",iqphantom.z

    zslicefig = plt.figure()
    ax1 = zslicefig.add_subplot(111)

    plt.xtitle = "Z slice"
    plt.ytitle = "Counts in slice"
    #plt.plot(iqphantom.z_slice())
    #plt.show()


    # images

    output['EARL z_slice'] = iqphantom.z
    output['EARL zprofile'] = zslicefig



    print "2 - registering data"

    print "COR and angle before registration:", iqphantom.cor, iqphantom.angle

    #regpars = iqphantom.register_phantom()
    #print "registration pars:", regpars[0], regpars[1]
   # output['regpars'] = regpars

    print "COR and angle after registration:", iqphantom.cor, iqphantom.angle

    output['cor'] = iqphantom.cor
    output['angle'] = iqphantom.angle

    print "3 - show result: "
    print "zslice", iqphantom.z

    output['EARL halfslice'] = iqphantom.plot_dataspheres(iqphantom.z)

    

    print "4 - sphere statistics: "
    results  = iqphantom.get_stats()
    for key in results.keys():
        print "Sphere: ", key 
        tmpresult = results[key]
        for subkey in tmpresult.keys():
            print '\t', subkey, tmpresult[subkey]
    bgresults  = iqphantom.get_bg_stats()
    for key in bgresults.keys():
        print "BG Sphere: ", key 
        tmpresult = bgresults[key]
        for subkey in tmpresult.keys():
            print '\t', subkey, tmpresult[subkey]
        

    bgaverage = np.average([bgresults[key]['avg'] for key in bgresults.keys()])
    print bgaverage

    meancontrastlist = {}
    maxcontrastlist = {}
    for key in results.keys():
        meancontrastlist[key] = (results[key]['avg'] - bgaverage)/bgaverage
        maxcontrastlist[key] = (results[key]['max'] - bgaverage)/bgaverage

    output['maxcontrast'] = maxcontrastlist
    output['meancontrast'] = meancontrastlist


    print "5 - contrasts: "
    for key in meancontrastlist.keys():
        print "mean:", key, meancontrastlist[key]
        print "max:", key, maxcontrastlist[key]


    print "6 - rc list: "
    
    inputcontrast = (actspsctime - actbgsctime)/ actbgsctime
    meanrclist = {}
    maxrclist = {}

    for key in results.keys():
        meanrclist[key] = meancontrastlist[key]/inputcontrast
        maxrclist[key] = maxcontrastlist[key]/inputcontrast

    for key in meanrclist.keys():
        print "mean: ",key, meanrclist[key]
        print "max: ", key, maxrclist[key]

    output['maxrc'] = maxrclist
    output['meanrc'] = meanrclist


    output['RC-1'] = meanrclist[1]
    output['RC-2'] = meanrclist[2]
    output['RC-3'] = meanrclist[3]
    output['RC-4'] = meanrclist[4]
    output['RC-5'] = meanrclist[5]
    output['RC-6'] = meanrclist[6]

    output['EARL rc mean curve'] =  plot_rclist(meanrclist,'za')
    output['EARL rc max curve'] =  plot_rclist(maxrclist,'za')


    output['status'] = 'processed'

    print "Writing results to wadserver" 



    return output
    
   



def main():
    parser = argparse.ArgumentParser(description='Commandline arguments to run reconstructor')
  
    parser.add_argument('-f',metavar='file',nargs='?',type=str)
    parser.add_argument('-m',metavar='mode',nargs='?',type=str)
    parser.add_argument('-ord',metavar='order',nargs='?',type=list)
    parser.add_argument('-ort',metavar='orientation',nargs='?',type=str)

    parser.add_argument('-bga',metavar='bg activity',nargs='?',type=float)
    parser.add_argument('-spa',metavar='sphere activity',nargs='?',type=float)

    parser.add_argument('-Tbga',metavar='Cal time bg activity',nargs='?',type=float)
    parser.add_argument('-Tspa',metavar='Cal time sphere activity',nargs='?',type=float)

    parser.add_argument('-s',metavar='slice',nargs='?',type=int)


    args = vars(parser.parse_args())
    filename = args['f'] #dicom directory

    m = args['m'] #mode can be max or min
    order = args['ord']
    ort = args['ort']
    bga = args['bga'] #background activity
    spa = args['spa'] #sphere activity

    order = args['ord']
    sid = args['s']

    #activity = args['A']
    #caltime = args['Ta']
    #remainder = args['Ar']
    #remtime = args['Tr']

    #length = args['len']
    #rad = args['rad']
    #volume = args['V']

    #dimensions = [length,rad]
    #print dimensions




    update_wad(filename = filename, mode = m, axorder = order, sporder = ort,bgac = bga, sphereac = spa, ctimespac = '120000', bgactime = '120000', spremac = 999, spremactime = '120000', bgremac = 999, bgremactime = '120000',scantime='120000')


if __name__ == "__main__":
    series_lst = [['./dicom.dcm']]
    plugin_data = PluginData(series_lst)
    sys.exit(main())
