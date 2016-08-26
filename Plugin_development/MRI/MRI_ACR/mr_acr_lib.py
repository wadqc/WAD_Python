# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# PyWAD is an open-source set of plugins for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# This package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN) and Arnold Schilham (UMCU) 
#
#
# Date: 
# Version: 0.1
# Authors: 
# Changelog:
#
#
# Description of this plugin:
#
# 
#
from pyWADLib import pydicom_series as dicomseries

import pylab as plt
from scipy.optimize import curve_fit
from dicom import tag

def sigmoid(x, x0, k):
    y = 1 / (1 + np.exp(-k*(x-x0)))
    #print "hallo"
    return y


def _2d_ellipse(shape,point,r1,r2=None,show='False'):
    if r2 == None:
        r2 = r1
    tmp = np.zeros(shape)
    for i in range(shape[0]):
        for j in range(shape[1]):
            dist = np.sqrt((float(i - point[0])/float(r1))**2 + (float(j - point[1])/float(r2))**2)
            if dist <= 1:
                tmp[i,j] = 1
    if show=='True':
        matshow(tmp)
    return ma.make_mask(tmp)

class ACR_analyzer:
    def __init__(self,path):

        self.dcmseries = dicomseries.read_files(path,True,True)

        self.tagdict = {'Protocol Name': ['0018', '1030'], 'Scan sequence': ['0018', '0020'], 'Acq Mat': ['0018', '1310'], 'Instance No': ['0020', '0013'], 'TR': ['0018', '0080'], 'TE': ['0018', '0081'], 'Image No': ['0051', '1008'], 
                        'Series No': ['0020', '0011'], 'PixelSize':['0028','0030']}

        self.datadict = {}
        self.fill_datadict()

        

    def qc_protocol(self):
        qc_protocol = {}
        return



    def print_datadict_keys(self):
        for key in self.datadict.keys():
            for subkey in self.datadict[key].keys():
                for subsubkey in self.datadict[key][subkey].keys():
                    print key,subkey,subsubkey


                

    def show_me_image(self,protocolname,seriesno,instanceno,show='False'):
        array = self.datadict[protocolname][seriesno][instanceno][0]
        if show == 'True':
            plt.matshow(array)
        return array

    def fill_datadict(self):
        for study in self.dcmseries:
            sequence  = study.__dict__['_datasets'] 
            tmpshape = study.__dict__['_shape']
            if len(tmpshape) == 2:
                nrofinstances = 1
            elif len(tmpshape) == 3:
                nrofinstances = tmpshape[0]
            for instancenr in range(nrofinstances):
                print "Instance -----", instancenr
                tmpinstance = sequence[instancenr]
                tmparray = tmpinstance.pixel_array
                
                tmpprotocol = str(tmpinstance[self.tagdict['Protocol Name']].value)
                tmpseriesno = str(tmpinstance[self.tagdict['Series No']].value)
                tmpdcminstance = str(tmpinstance[self.tagdict['Instance No']].value)
                
                tmpdata = [tmparray,self.metadata(tmpinstance)]

                if tmpprotocol not in self.datadict.keys():
                    self.datadict[tmpprotocol] = {tmpseriesno:{tmpdcminstance:tmpdata}}
                else:
                    if tmpseriesno not in self.datadict[tmpprotocol].keys():
                        self.datadict[tmpprotocol][tmpseriesno] = {tmpdcminstance:tmpdata}
                    else:
                        if tmpdcminstance not in self.datadict[tmpprotocol][tmpseriesno].keys():
                            self.datadict[tmpprotocol][tmpseriesno][tmpdcminstance] = tmpdata                        

    def metadata(self,dcminstance):
        outdict = {}
        for key in self.tagdict.keys():
            print key
            outdict[key] = dcminstance[tag.Tag(self.tagdict[key])]
            
    
    def slice_position_accuracy(self,tmparray):

        threshold = np.mean(tmparray[200:210,100:110])
        
        y_offset = 100
        x_offset = 200
        cross_profile = np.mean(tmparray[y_offset:y_offset+10,x_offset:x_offset+100],0)        

        points = []

        # determine where the half threshold is being crossed

        sign = 1
        for i in range(len(cross_profile)):
            if np.abs(cross_profile[i] - threshold/2) < 500:
                points.append(i)
                
        

        print points
        # define left and right profiles 
        midpoint = x_offset + (np.max(points) + np.min(points))/2
        print midpoint
        width_profile = 5
        mid_off = 2

        left_profile = np.mean(tmparray[80:150,midpoint-width_profile-mid_off:midpoint-mid_off],1)
        right_profile = np.mean(tmparray[80:150,midpoint+mid_off:mid_off+midpoint + width_profile],1) 

        # Both profiles have identical shapes but are slightly shifted due to the diff. in lengths
        # to calculate the difference in lengths we fit each curve to a sigmoid curve and from
        # the difference we determine the difference in lengths

        
        p0 = [45,0.5] #initial estimate
        
        lnorm = left_profile/np.max(left_profile)
        rnorm = right_profile/np.max(right_profile)

        x = np.arange(len(lnorm))

        lopt, lcov = curve_fit(sigmoid, x, lnorm,p0)
        ropt, rcov = curve_fit(sigmoid, x, rnorm,p0)

        return np.abs(lopt[0] - ropt[0])


    def geometric_accuracy(self,tmparray):
        return



    def slice_thickness_accuracy(self,tmparray):

        threshold = np.mean(tmparray[200:210,100:110])
        
        y_offset = 100
        x_offset = 200
        cross_profile = np.mean(tmparray[y_offset:y_offset+10,x_offset:x_offset+100],0)        

        points = []

        # determine where the half threshold is being crossed

        sign = 1
        for i in range(len(cross_profile)):
            if np.abs(cross_profile[i] - threshold/2) < 500:
                points.append(i)
                
        

        print points
        # define left and right profiles 
        midpoint = x_offset + (np.max(points) + np.min(points))/2
        print midpoint
        width_profile = 5
        mid_off = 2

        top_profile = np.mean(tmparray[80:150,midpoint-width_profile-mid_off:midpoint-mid_off],1)
        bottom_profile = np.mean(tmparray[80:150,midpoint+mid_off:mid_off+midpoint + width_profile],1) 

        # Both profiles have identical shapes but are slightly shifted due to the diff. in lengths
        # to calculate the difference in lengths we fit each curve to a sigmoid curve and from
        # the difference we determine the difference in lengths

        
        p0 = [45,0.5] #initial estimate
        
        lnorm = left_profile/np.max(left_profile)
        rnorm = right_profile/np.max(right_profile)

        x = np.arange(len(lnorm))

        lopt, lcov = curve_fit(sigmoid, x, lnorm,p0)
        ropt, rcov = curve_fit(sigmoid, x, rnorm,p0)

        return np.abs(lopt[0] - ropt[0])


    def ghosting(self,data=None,showv='False'):

        if data == None:
            data = self.show_me_image('t1_se_tra','3','7',show='False')


        mymaskdict = {}
        mymaskdict['center']= _2d_ellipse(shape(data),[256,256],164.,164.,show=showv)
        mymaskdict['top']= _2d_ellipse(shape(data),[35,250],20.,100.,show=showv)
        mymaskdict['right']= _2d_ellipse(shape(data),[250,480],100.,20.,show=showv)
        mymaskdict['bottom']= _2d_ellipse(shape(data),[475,250],20.,100.,show=showv)
        mymaskdict['left']= _2d_ellipse(shape(data),[250,30],100.,20.,show=showv)

        return self.calculate_ghosting(data,mymaskdict)

    def calculate_ghosting(self,data,roimaskdict, show='False'):
        results = {}
        for key in roimaskdict.keys():
            marr = ma.masked_array(data,mask=1-roimaskdict[key])
            tmpsize = ma.count(marr)
            tmpmean = ma.mean(marr)
            tmpstd = ma.std(marr)
            results[key] = [tmpmean,tmpstd,tmpsize]
            if show=='True':
                plt.matshow(marr)
            
            
        ghosting = ((results['top'][0]+results['bottom'][0])- (results['left'][0]+results['right'][0]))/results['center'][0]
        results['ghosting'] = ghosting

        return results
