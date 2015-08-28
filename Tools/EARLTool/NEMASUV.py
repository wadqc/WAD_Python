
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
        

#---------


def loadData(path):
#    dirlist = [filename for filename in os.listdir(path) if filename.endswith(".IMA")]    
    dirlist = [filename for filename in os.listdir(path) ]    

    first = True
    info = {}
    slice_locations=[]
    slice_index = []
    for filename in dirlist:
        ### Load data
        print filename

        dc=dicom.read_file(os.path.join(path,filename))
        
        if first:
            first = False
            print dc
            info["rescale_slope"] = float(dc.get('RescaleSlope'))
            info["pixel_spacing"] = float(dc.get('PixelSpacing')[0])
            info["SliceThickness"] = float(dc.get('SliceThickness'))
            info["ActualFrameDuration"] = float(dc.get('ActualFrameDuration'))
            info["AcquisitionDate"] = dc.get('AcquisitionDate')
            print "wtp????"
            wtp = 100 #float(dc.get("Radiopharmaceutical_Information_Sequence")[0].get("Radionuclide_Total_Dose"))
            print wtp
            info["RadionuclideTotalDose"] = wtp
            nr_of_slices =  111 #(int(dc.get('NumberOfSlices')))

            rows = (int(dc.get('Rows')))
            cols = (int(dc.get('Columns')))
            if rows != cols:
                print "Slices niet vierkant!"
            data = numpy.empty([rows, rows, nr_of_slices],dtype=numpy.int32)

        slice_index.append(int(dc.get('ImageIndex'))-1)
        slice_locations.append(float(dc.get('SliceLocation')))
        data[:,:,slice_index[-1]] = dc.pixel_array
        
    ### Determine boundaries
    ind = numpy.argsort(slice_index)
    info["slice_locations"] = numpy.array(slice_locations)[ind]

    return [data, info]

       axial_center = numpy.median(info["slice_locations"])
        ax_low = axial_center - const["AxialRange"]/2
        ax_high = axial_center + const["AxialRange"]/2
        use_slices = [i for i, loc in enumerate(info["slice_locations"]) if (ax_low <= loc <= ax_high)]

        pict_size = data.shape[0]
        xmin = int(pict_size*(1.0 - const["CropHeight"])/2.0)
        xmax = int(pict_size*(1.0 + const["CropHeight"])/2.0)

        r_size = info["pixel_spacing"]
        z_size = abs(info["slice_locations"][0]-info["slice_locations"][1])

        center_slice = use_slices[int(len(use_slices)/2)]
        vmax = float(numpy.max(data[xmin:xmax,xmin:xmax,center_slice ]))
        vmin = 0.0

        max_factor = 2.0
        def_factor = 1.2

     rad_pixel_bound = const["Diameter"]/info["pixel_spacing"]/2.0
        square_size = max(1, int(const["SquareSize"]/info["pixel_spacing"]))


        const = {"Path": self.path.get(),
                 "AxialRange": float(self.cnf["AxialRange"]),
                 "Diameter": float(self.cnf["Diameter"]),
                 "SquareSize": float(self.cnf["SquareSize"]),
                 "Activity": float(self.QC_activity.get()),
                 "GeGaVolume": float(self.cnf["GeGaVolume"]),
                 "CropHeight": float(self.cnf["CropHeight"])}

def calculate(data, info, const, use_slices, center_slice, xmin, rad_pixel_bound, square_size):
    allcounts = []
    center_of_mass = []
    squares_pos = []
    vertices = [[],[]]
    for index in use_slices:
        ### Determine center for each slice and fill with squares up to the radial boundary
        center_of_mass.append( get_center_of_mass(data[:,:,index])[::-1] )
        squares_pos.append( getSquares(data.shape[0], center_of_mass[-1], rad_pixel_bound, square_size) )

        ### List the total nr of counts in each square, in each slice
        allcounts.append([])
        for pos in squares_pos[-1]:
            allcounts[-1].append(numpy.sum(data[pos[0]:pos[0]+square_size,pos[1]:pos[1]+square_size,index]))

        ### Create bounding box in image 1 and 2. Indicates the area used for calculation.
        for i, verts in enumerate(vertices):
            start = (index, center_of_mass[-1][i]-rad_pixel_bound-xmin)
            end = (start[0]+1,start[1])
            verts+=[start,end]

        if index == center_slice:
            center_slice_index = len(allcounts)-1

    for verts in vertices:
        verts += [(v[0],v[1]+2.0*rad_pixel_bound) for v in reversed(verts)]
        verts.insert(0, verts[0])
        verts.append(verts[0])
    codes = [Path.LINETO for v in vertices[0]]
    codes[0] = Path.MOVETO
    codes[-1] = Path.CLOSEPOLY

    return [allcounts, squares_pos, center_of_mass, vertices, codes]

def getNU(counts):
        av = float(numpy.average(counts))
        mx = float(numpy.max(counts))
        mn = float(numpy.min(counts))
        return [(mx-av)/av*100,(mn-av)/av*100]
    
def getResults(info, const, allcounts, square_size):#, use_slices, center_slice_index, r_size, z_size, square_size):
    print "\n--- Results ---"
    results = {}

    
    ### Average slice counts
    #Vpix = info["SliceThickness"] * r_size**2 * 0.001
    Time = info["ActualFrameDuration"]*1e-3
    '''
    avcounts = numpy.average(allcounts[center_slice_index])/square_size**2

    print "RescaleSlope : ", info["rescale_slope"]
    print "Slice-Thickness, Slice-Spacing, Pixel-Spacing: ", info["SliceThickness"], z_size, r_size
    act = 1000 * const["Activity"] / const["GeGaVolume"]
    print "C/s / mL         : ", avcounts/(Time*Vpix)
    print "Activity / mL    : ", act
    print "Ratio            : ", avcounts/(Time*Vpix)/act
    '''
    
    ### Slice uniformity
    NU_slice = []
    CV_slice = []
    for slicecounts in allcounts:
        NU_slice.append(getNU(slicecounts))
        CV_slice.append(coefficient_of_variation(slicecounts)*100)
    showdata = numpy.transpose(numpy.append(numpy.transpose(NU_slice),[CV_slice], axis=0))
    results["SlcNU_max"] = numpy.max(numpy.transpose(NU_slice)[0])
    results["SlcNU_min"] = numpy.max(numpy.transpose(NU_slice)[1])
    results["SlcCV"] = numpy.max(CV_slice)

    ### Volume uniformity
    flattened = [counts for slicecounts in allcounts for counts in slicecounts]
    NU_volume = getNU(flattened)
    CV_volume = coefficient_of_variation(flattened)*100
    print "Volume: ", NU_volume, CV_volume
    results["VolNU_max"] = NU_volume[0]
    results["VolNU_min"] = NU_volume[1]
    results["VolCV"] = CV_volume
    
    ### System uniformity
    av_counts = [numpy.average(slicecounts) for slicecounts in allcounts]
    NU_system = getNU(av_counts)
    CV_system = coefficient_of_variation(av_counts)*100
    print "System: ", NU_system, CV_system
    results["SysNU_max"] = NU_system[0]
    results["SysNU_min"] = NU_system[1]
    results["SysCV"] = CV_system
    
    ### SUV = counts/pixel * rescale_slope / (activity*10^6 / volume*10^3) = counts * SUV_factor
    SUV_factor = 0.001 * info["rescale_slope"] * const["GeGaVolume"] / (const["Activity"] * square_size**2)
    SUV = float(numpy.average(flattened))*SUV_factor
    print "SUV average      : ", SUV
    results["SUV_av"] = SUV
    
    ### Measured activity (relying on rescale_slope).
    decays = SUV * const["Activity"] * Time
    results["Decays"] = decays
    print "Expected decays  : ", const["Activity"] * Time
    print "Decays (1E6)     : ", decays
    print "---------------\n"
    
    return [results, showdata]

def get_center_of_mass(image):
    ### Copy from scipy.ndimage.measurements.center_of_mass (copied to reduce compiled size, prevents importing scipy)
    normalizer = numpy.sum(image)
    grids = numpy.ogrid[[slice(0, i) for i in image.shape]]
    results = [numpy.sum(image * grids[dir].astype(float)) / normalizer for dir in range(image.ndim)]
    return tuple(results)     

def coefficient_of_variation(lijst):
    return numpy.std(lijst)/numpy.average(lijst)
    
def getSquares(pict_size, center_of_mass, radius, square_size):
    ### Returns a list of square locations (lower-left corner) within the bounded radius.
    #Only check within (minx,miny) and (maxx,maxy). Integer multiples of the square_size, -1.
    minx = int((center_of_mass[0]-radius)/square_size)*square_size
    miny = int((center_of_mass[1]-radius)/square_size)*square_size
    maxx = int((center_of_mass[0]+radius)/square_size)*square_size + square_size
    maxy = int((center_of_mass[1]+radius)/square_size)*square_size + square_size

    #Don't check within center square: Certain to fall within circle.
    max_minx = int(center_of_mass[0]-(0.7*radius))+1
    min_maxx = int(center_of_mass[0]+(0.7*radius))
    max_miny = int(center_of_mass[1]-(0.7*radius))+1
    min_maxy = int(center_of_mass[1]+(0.7*radius))

    radsq = numpy.square(radius)
    mat = numpy.zeros((maxx-minx, maxy-miny), dtype=numpy.int8)
    #For every pixel within (minx,miny)and (maxx,maxy), check whether it lies in the circle.
    for x in range(minx,maxx):
        for y in range(miny,maxy):
            if max_minx < x < min_maxx and max_miny < y < min_maxy:
                mat[x-minx, y-miny] = 1
            else:
                dist = (float(x)-center_of_mass[0])**2 + (float(y)-center_of_mass[1])**2
                if dist < radsq:
                    mat[x-minx,y-miny]=1
            
    squares_pos=[]
    s = square_size - 1
    for x in range(0,maxx-minx,square_size):
        for y in range(0,maxy-miny,square_size):
            if mat[x,y] and mat[x+s,y] and mat[x,y+s] and mat[x+s,y+s]:
                squares_pos.append([minx+x,miny+y])
                
    return squares_pos


#---------







    def plot_img_hist(self, sid, bins=256,showfig=False):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        tmpdat = self.data[sid,:,:]
        
        plt.hist(tmpdat.flatten(), bins, range=(tmpdat.min(),tmpdat.max()), fc='k', ec='k')
        if showfig==True:
            plt.show()
        return fig

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
