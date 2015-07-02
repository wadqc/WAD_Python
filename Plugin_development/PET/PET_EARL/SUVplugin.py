#!/usr/bin/env python 
print '...'
import numpy
import pylab
from matplotlib.patches import Rectangle, Circle, PathPatch
from matplotlib.path import Path
import os, os.path
import sys
import dicom

import time
import datetime

#CONFIGFILE = ".\\config.txt"
CONFIGFILE = "config.txt"

class SUVanalyzer:
    def __init__(self,data,info):


        config_variables = ["DefaultPath",
                            "ImageSize","CropHeight","GraphHeight",
                            "AxialRange","Diameter","SquareSize",
                            "GermaniumHalflife","GeGaCalibrationDate",
                            "GeGaCalibrationActivity","GeGaVolume",
                            "RequiredDecays", "RDMargin",
                            "NU_limit","CV_limit","SUV_limit","Counts_limit"]
        self.cnf = self.loadConfig(config_variables)


        self.QC_date = datetime.datetime.strptime(info['AcquisitionDate'],"%Y%m%d").strftime("%d-%m-%Y")
        #now = datetime.datetime.now()
        #today_string = "-".join([str(now.day),str(now.month),str(now.year)])
        #self.QC_date = '01-05-2012' #today_string


        self.QC_activity = ''
        self.QC_duration = ''

        self.Dicom_activity = ''
        self.Dicom_duration = ''
        

        self.results = {}
        self.showdata = []


        self.data = data
        self.info = info
        self.calculated = {}

        self.getActivity()


        print self.QC_date
        print self.QC_activity
        print self.QC_duration

        print '*'*20
        print '*'*20
        print self.calculated.keys()
        print '*'*20
        print '*'*20

        print self.calculated['center_slice']

        pylab.imshow(self.data[:,:,self.calculated['center_slice']])
        pylab.show()
        
        pylab.imshow(self.data[:,200,:])
        pylab.show()

        pylab.imshow(self.data[200,:,:])
        pylab.show()

        print '*'*20
        print '*'*20
        #print self.results
        print '*'*20
        print '*'*20

        


        xsize = 300
        ysize = 240
        fig = pylab.figure(figsize=(xsize/100.0,ysize/100.0), dpi=100)
        ax = fig.add_subplot(111)
        ax.grid(True)
        ax.plot(self.showdata) #COV curve
        ax.plot(self.results["SUV_curve"]) #COV curve
        pylab.show()


    def getActivity(self):
        halflife = float(self.cnf["GermaniumHalflife"]) * 86400
        cal_activity = float(self.cnf["GeGaCalibrationActivity"])
        cal_date = self.cnf["GeGaCalibrationDate"].split("-")
        qc_date = self.QC_date.split("-")
        activity = calcActivity(cal_date, qc_date, halflife, cal_activity)
        self.QC_activity = "%.2f"%activity
        
        settime = float(self.cnf["RequiredDecays"])/activity/60.0
        self.QC_duration = "%.2f"%settime
        
        [data, info] = [self.data,self.info]
        
        qcd = info["AcquisitionDate"]
        self.QC_date = qcd[6:]+"-"+qcd[4:6]+"-"+qcd[:4]


        dicom_activity = info["RadionuclideTotalDose"]
        dicom_activity = dicom_activity * 1E-6
        dicom_duration = info["ActualFrameDuration"]
        dicom_duration = (dicom_duration * 1E-3)/60                

        '''
        self.dicom_activity_label.config(text="")
        self.dicom_duration_label.config(text="")
        if abs(dicom_activity - float(self.QC_activity.get())) > 0.02:
            self.dicom_activity_label.config(text="%.2f MBq volgens dicom header"%dicom_activity)
        if abs(dicom_duration - float(self.QC_duration.get())) > 0.02:
            self.dicom_duration_label.config(text="%.2f min volgens dicom header"%dicom_duration)
        '''

        const = {"AxialRange": float(self.cnf["AxialRange"]),
                 "Diameter": float(self.cnf["Diameter"]),
                 "SquareSize": float(self.cnf["SquareSize"]),
                 "Activity": float(self.QC_activity),
                 "GeGaVolume": float(self.cnf["GeGaVolume"]),
                 "CropHeight": float(self.cnf["CropHeight"])}
        
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

        t1=time.time()
        [allcounts, squares_pos, center_of_mass, vertices, codes] = calculate(data, info, const, use_slices, center_slice, xmin, rad_pixel_bound, square_size)
        print "Calculating %.4f sec."%(time.time()-t1)


        
        [self.results, self.showdata] = getResults(info, const, allcounts, square_size)
   
        ### self.calculated holds information needed to draw the patches in the images
        self.calculated["done"] = True
        self.calculated["use_slices"]= use_slices
        self.calculated["center_slice"] = center_slice
        self.calculated["center_of_mass"]= center_of_mass
        self.calculated["xmin"] = xmin
        self.calculated["rad_pixel_bound"] = rad_pixel_bound
        self.calculated["squares_pos"] = squares_pos
        self.calculated["square_size"] = square_size


            
                        
        
    def loadConfig(self, config_variables):
        cnf = {}
        error = None
        if not os.path.exists(CONFIGFILE):
            open(CONFIGFILE, 'w').close()
            
        fp = open(CONFIGFILE, "r+")
        for line in fp:
            line = line.strip().split()
            if len(line) and line[0] in config_variables:
                if len(line) == 1 or line[1][0] == "#":
                    error=line[0]
                else:
                    cnf[line[0]] = line[1]
        
        for var in config_variables:
            if var not in cnf:
                fp.write("\n"+var+"\n")
                error = var

        fp.close()        
                
        if error != None:
            print "\n"+20*"*"
            print "Geen waarde voor " + error + " in'" + CONFIGFILE + "'!"
            print 20*"*"+"\n"
            sys.exit()
        return cnf

    
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

    

    SUV_factor = 0.001 * info["rescale_slope"] * const["GeGaVolume"] / (const["Activity"] * square_size**2)
    print "eeeee",numpy.shape(allcounts)
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

    SUV = float(numpy.average(flattened))*SUV_factor
    print "SUV average      : ", SUV
    results["SUV_av"] = SUV
    print numpy.shape(flattened)

   
 
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


def getDicom(path, string):
  #  dirlist = [filename for filename in os.listdir(path) if filename.endswith(".IMA")]
    dirlist = [filename for filename in os.listdir(path) ]
    if len(dirlist):
        dc=dicom.read_file(os.path.join(path, dirlist[0]))
        return dc.get(string)

def calcActivity(cal_date, qc_date, halflife, cal_activity):
        sec = []
        for date in [cal_date, qc_date]:
            pydate = datetime.datetime(int(date[2]), int(date[1]), int(date[0]), 12, 0)
            sec.append(time.mktime(pydate.timetuple()))
        decay = (sec[1]-sec[0])/halflife
        return cal_activity * 0.5**decay

def suv_main(data,results,**kwargs):

    
    first = True
    info = {}
    slice_locations=[]
    slice_index = []

    filelst = data.getAllInstances()
    for dcmfile in filelst:
        ### Load data
    
        dc=dcmfile
        
        if first:
            first = False
            #print dc
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

    print info
    SUVanalyzer(data,info)





if __name__ == "__main__":
    app = SUVanalyzer('/data/SoftwareTools/PythonScripts/PETQC/TEST2')

