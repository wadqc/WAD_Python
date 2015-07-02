'''
This code borrows extensively from pylinac software:

https://pypi.python.org/pypi/pylinac

PyTG142 (pytg) provides TG-142 quality assurance (QA) tools to Python programmers in the field of therapy medical
physics. This package's intention is to create a free and open-source alternative to TG-142 QA commercial solutions.
'''


from __future__ import print_function
import numpy as np
import os, os.path as osp
from dicom import tag

if not 'MPLCONFIGDIR' in os.environ:
    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.


import scipy.ndimage.measurements as meas
import dicom
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt



import common_functions as ptg

__version__='03062015'


"""Default constants"""
nominal_HU = [-1000, -200, -100, -35, 120, 340, 990]
nominal_GEO = (50.0)
tolerances = {'HU':38, 'UN':40, 'SR':None, 'GEO':1}


field_of_view = ['small', 'large']
FOV_thresh = 0  # 300 the field of view size threshold between small and large
small_FOV_settings = {'size_thresh':100000, 'pix_thresh':500,
                      'GEO_box_bounds':[220, 220, 180+110, 180+110], 'GEO_min_size':25,
                      'HU_radius':50, 'HU_angles':[0, -2.62, -1.575, -0.525, 0.53, 1.56, -3.67], 'HU_boxsize':8,
                      'CON_angles':[0.4, -0.-5, 2.03], 'CON_radius':100, 'CON_ROIsize':20,
#                      'SR_radius':[34.1, 45.5], 'SR_theta':[(-3.77, -3.47), (-4.16, -3.91), (-4.57, -4.35), (-4.89, -4.69), (-5.2, -5.02), (-5.5, -5.35), (-5.78, -5.67), (-6.08, -5.97)],
                      'UN_boxsize':30, 'UN_angles':[np.pi/2, np.pi*3/2, np.pi, 0], 'UN_radius':60,
                      'SR_radius':[94.1, 100.5], 'SR_theta':[(-3.77, -3.47), (-4.16, -3.91), (-4.57, -4.35), (-4.89, -4.69), (-5.2, -5.02), (-5.5, -5.35), (-5.78, -5.67), (-6.08, -5.97)]}#,


class CATPHAN(object):
    """
    A class for loading and analyzing Cone-Beam CT DICOM files of a CatPhan 500. Analyzes: Uniformity,
    Spatial Resolution, Contrast, & HU Linearity.
    #TODO: Low Contrast
    """
    def __init__(self, using_pyqa=False):

        self.slices_of_interest = {'HU':14, 'UN':35, 'SR':20 } #, 'LOCON':23} The Catphan 503 does not have  a low contrast module
        self.images = []
        self.slice_centers = {}
        self.SOI_cleaned = {}
        self.FOV = ''
        self.UN_roi_bounds = []
        self.HU_roi_bounds = []
        self.MTF = {'LP/mm':[0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.41, 1.59], # the line pair/mm values for the CatPhan SR slice
                    'LP points': np.zeros(8), # the point values for the MTF function
                    'interp_x': np.zeros(139), # an interpolated (0.01) array of the line pair values; used for finding MTF_50, etc.
                    'interp_y': np.zeros(139),  # an interpolated (0.01) array of the actual MTF
                    'ROI points': [], # point values for where the SR ROIs were taken; used for plotting
                    'MTF_50': 0}  # the MTF 50% point in lp/mm
        self.GEO_dist = {} # dictionary of distances of geometric nodes; will be 4 values
        self.GEO_diff = {} # dictionary of differences of geometric distances from nominal
        self.results = {} # dictionary of ROI value results
        self.results_diffs = {} # dictionary of differences between nominal values and results
        self.UN_val = np.zeros(5)
        self.UN_diff = np.zeros(5)
        self._roll_found = False
        self.using_pyqa = using_pyqa


    def _load_files(self, file_list):
        """
        Load CT DICOM files given a list of images
        :return:
        """
        # intialize some lists
        images = []
        im_order = []
        # check that enough slices were loaded
        if len(file_list) < self.slices_of_interest['SR']:
            if self.using_pyqa:
                return #TODO: post message box
            else:
                raise ValueError, "Not enough slices were selected. Select all the files or change the SR slice value."
        # load dicom files from list names
        maxlist = []
        for item in file_list:
            dcm = item
            _tmparray = dcm._get_pixel_array()
            rescaleintercept = int(dcm[tag.Tag("0028","1052")].value)
            rescaleslope = int(dcm[tag.Tag("0028","1053")].value)

            tmparray =  _tmparray*rescaleslope + rescaleintercept

            images.append(tmparray)
            
            tmporder = np.round(dcm.ImagePositionPatient[-1]/dcm.SliceThickness)
            im_order.append(tmporder)

            maxlist.append([np.max(tmparray),tmporder])

        self.improps = {'FOV': dcm.DataCollectionDiameter,
                        'width': np.size(images[0],0),
                        'height': np.size(images[0],1),
                        'FOV size':field_of_view[0] if dcm.DataCollectionDiameter < FOV_thresh else field_of_view[1],
                        'Rescale Slope':dcm.RescaleSlope,
                        'Rescale Intercept':dcm.RescaleIntercept,
                        'mm/Pixel':dcm.PixelSpacing[0]}

        self._sort_images(images, im_order)
        


    def _sort_images(self, image_list, order_list):
        """
        Sort the CBCT files into the correct order based on DICOM tags
        """
        # preallocate
        #self.images = np.zeros((self.improps['width'], self.improps['height'], len(image_list)))

        print('len imglst',len(image_list))

        self.images = np.zeros((len(image_list),self.improps['width'], self.improps['height']))

        # shift order list from -xxx:+xxx to 0:2xxx
        order_list += np.round(abs(min(order_list)))

        print ('ordl',sorted(order_list))
        
        for idx, item in enumerate(order_list.astype(int)):  #TODO: possible list comp

            print("index",idx,item)
            self.images[item,:,:] = image_list[idx]
           
        

    def find_phan_centers(self,params):
        """
        Find the center of the phantom for the slices of analysis
        :return:
        """

        self.slices_of_interest =  ptg.find_slices_of_interest(self.images,params)
        print('haleluja!',self.slices_of_interest)
        for key in self.slices_of_interest:
        
            print ("curr SOI",key)
            #plt.imshow(self.images[self.slices_of_interest[key],:,:])
            #plt.title("%s"%key)
            #plt.show()

            SOI = self.images[self.slices_of_interest[key],:,:]
            SOI = np.where(SOI > np.percentile(SOI,48), 1, 0)  #  convert slice to binary based on threshold
            #SOI = morph.binary_fill_holes(SOI)  # fill in air pockets to make one solid ROI
            SOI_labeled, no_roi = meas.label(SOI)  # identify the ROIs
            if no_roi < 1 or no_roi is None:
                raise ValueError, "Unable to locate the CatPhan"
            hist,bins = np.histogram(SOI_labeled,bins=no_roi) # hist will give the size of each label
            SOI = np.where(SOI_labeled == np.argmax(hist), 1, 0)  # remove all ROIs except the largest one (the CatPhan)
            self.slice_centers[key] = meas.center_of_mass(SOI)
            self.SOI_cleaned[key] = SOI

    def find_roll(self):
        """
        Determine the roll of the phantom based on the phantom air bubbles on the HU slice.
        :return:
        """
        SOI = self.SOI_cleaned['HU']
        invSOI = ptg.invert(SOI)
        labels, no_roi = meas.label(invSOI)
        roi_sizes = [meas.sum(invSOI, labels, index=item) for item in range(1,no_roi+1)]
        air_bubbles = [idx+1 for idx, item in enumerate(roi_sizes) if item < np.median(roi_sizes)*1.5 and item > np.median(roi_sizes)*(1/1.5)]
        if len(air_bubbles) != 2:
            self.phan_roll = 0

            #raise RuntimeWarning, "Roll unable to be determined; assuming 0"
        else:
            air_bubble_CofM = meas.center_of_mass(invSOI, labels, air_bubbles)
            y_dist = air_bubble_CofM[0][0] - air_bubble_CofM[1][0]
            x_dist = air_bubble_CofM[0][1] - air_bubble_CofM[1][1]
            angle = np.arctan2(y_dist, x_dist) * 180/np.pi
            if angle < 0:
                roll = abs(angle) - 90
            else:
                roll = angle - 90
            self.phan_roll = roll

            self._roll_found = True

    def find_HU(self):
        """
        Determine HU values from HU slice. Averages 3 slices.
        :return:
        """
        # preallocation
        HU_val = np.zeros((len(nominal_HU)))
        HU_diff = np.zeros((len(nominal_HU)))

        # average 3 slices around the nominal slice
        HU_slice = np.mean(self.images[self.slices_of_interest['HU']-1:self.slices_of_interest['HU']+1,:,:],0)

        # For each HU ROI...
        for idx, HU in enumerate(nominal_HU):
            # calculate the box bounds of the ROI
            self.HU_roi_bounds.append([self.slice_centers['HU'][0] + small_FOV_settings['HU_radius']*np.cos(np.radians(self.phan_roll)+small_FOV_settings['HU_angles'][idx]),
                          self.slice_centers['HU'][1] + small_FOV_settings['HU_radius']*np.sin(np.radians(self.phan_roll) + small_FOV_settings['HU_angles'][idx])])
            # sample the ROI
            HU_roi = HU_slice[self.HU_roi_bounds[idx][0] - small_FOV_settings['HU_boxsize']/2:self.HU_roi_bounds[idx][0] + small_FOV_settings['HU_boxsize']/2,
                     self.HU_roi_bounds[idx][1] - small_FOV_settings['HU_boxsize']/2:self.HU_roi_bounds[idx][1] + small_FOV_settings['HU_boxsize']/2]
            # Calculate the mean HU value and the difference from nominal
            HU_val[idx] = np.mean(HU_roi*self.improps['Rescale Slope'] + self.improps['Rescale Intercept'])
            HU_diff[idx] = HU_val[idx] - HU

        self.results['HU'] = HU_val
        self.results_diffs['HU'] = HU_diff


        if all(self.results_diffs['HU'] < tolerances['HU']):
            self.HU_pass = True
        else:
            self.HU_pass = False



    def find_SR(self, plot_ROIs=False):
        """
        Determine the Spatial Resolution from the Line-Pair slice. Averages 3 slices.
        :return:
        """
        # take max from 3 slices around the nominal slice
        SR_slice = np.max(self.images[self.slices_of_interest['SR'] - 1:self.slices_of_interest['SR'] + 1,:,:], 0)
        # pilSRslice = Image.fromarray(SR_slice)
        if plot_ROIs:
            fig, ax = plt.subplots(1,8)

        # create y,x coordinates
        values = []
        SRboxes = []
        for idx, (left, right) in enumerate(small_FOV_settings['SR_theta']):
            # for in_out in small_FOV_settings['SR_radius']:
            #TODO: could probably be optimized
            SRboxes.append([[self.slice_centers['SR'][0] + small_FOV_settings['SR_radius'][0] * np.sin(np.radians(self.phan_roll) + left), # inner left y-coordinate
                             self.slice_centers['SR'][0] + small_FOV_settings['SR_radius'][1] * np.sin(np.radians(self.phan_roll) + left),  # outer left y-coordinate
                             self.slice_centers['SR'][0] + small_FOV_settings['SR_radius'][1] * np.sin(np.radians(self.phan_roll) + right),  # outer right y-coordinate
                             self.slice_centers['SR'][0] + small_FOV_settings['SR_radius'][0] * np.sin(np.radians(self.phan_roll) + right),  # inner right y-coordinate
                            #repeat first point so that box completes when plotted
                            self.slice_centers['SR'][0] + small_FOV_settings['SR_radius'][0] *  np.sin(np.radians(self.phan_roll) + left)],  # inner left y-coordinate

                            [self.slice_centers['SR'][1] + small_FOV_settings['SR_radius'][0] * np.cos(np.radians(self.phan_roll) + left),  # inner left x-coordinate
                             self.slice_centers['SR'][1] + small_FOV_settings['SR_radius'][1] * np.cos(np.radians(self.phan_roll) + left),  # outer left x-coordinate
                             self.slice_centers['SR'][1] + small_FOV_settings['SR_radius'][1] * np.cos(np.radians(self.phan_roll) + right),  # outer right x-coordinate
                             self.slice_centers['SR'][1] + small_FOV_settings['SR_radius'][0] * np.cos(np.radians(self.phan_roll) + right), # inner right x-coordinate

                             self.slice_centers['SR'][1] + small_FOV_settings['SR_radius'][0] * np.cos(np.radians(self.phan_roll) + left)]])  # inner left x-coordinate


            #TODO: this could use some serious optimization
            big_mask = self.sector_mask(SR_slice.shape, (self.slice_centers['SR']), small_FOV_settings['SR_radius'][1], (left, right))
            little_mask = self.sector_mask(SR_slice.shape, (self.slice_centers['SR']), small_FOV_settings['SR_radius'][0], (left, right))
            masked_img = np.where(big_mask == True, SR_slice, 0)
            masked_img2 = np.where(little_mask == True, 0, masked_img)
            final_mask = np.where(masked_img2 == 0, np.NaN, masked_img2)
            if plot_ROIs:
                ax[idx].imshow(final_mask)
            values.append(np.percentile(final_mask[~np.isnan(final_mask)],[5,95]))

            self.MTF['LP points'][idx] = (values[idx][1] - values[idx][0]) / (values[0][1] - values[0][0])
            # ttt = np.where(big_mask == True, little_mask == False, SR_slice, 0)
            # ax[idx].imshow(masked_img2)

        self.MTF['ROI points'] = SRboxes
        self.MTF['interp_x'] = np.arange(self.MTF['LP/mm'][0], self.MTF['LP/mm'][-1], 0.01)
        self.MTF['interp_y'] = np.interp(self.MTF['interp_x'],self.MTF['LP/mm'],self.MTF['LP points'])

        self.MTF['MTF_50'] = self.MTF['interp_x'][np.argmin(np.abs(self.MTF['interp_y'] - 0.5))]

        pass

    def sector_mask(self, shape, centre, radius, angle_range):
        """
        Return a boolean mask for a circular sector. The start/stop angles in
        `angle_range` should be given in clockwise order.
        """

        x, y = np.ogrid[:shape[0], :shape[1]]
        cx, cy = centre
        # tmin, tmax = np.deg2rad(angle_range)
        tmin, tmax = angle_range

        # ensure stop angle > start angle
        if tmax < tmin:
            tmax += 2 * np.pi

        # convert cartesian --> polar coordinates
        r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
        theta = np.arctan2(x - cx, y - cy) - tmin

        # wrap angles between 0 and 2*pi
        theta %= (2 * np.pi)

        # circular mask
        circmask = r2 <= radius * radius

        # angular mask
        anglemask = theta <= (tmax - tmin)

        return circmask * anglemask

    def find_GEO(self):
        """
        Determine the geometric distortion. Averages 3 slices around HU slice.
        :return:
        """

        # average 3 slices around the nominal slice

        GEO_slice = np.mean(self.images[
self.slices_of_interest['HU'] - 1:self.slices_of_interest['HU'] + 1,
small_FOV_settings['GEO_box_bounds'][0]:small_FOV_settings['GEO_box_bounds'][2],
small_FOV_settings['GEO_box_bounds'][1]:small_FOV_settings['GEO_box_bounds'][3]], 0)

        
        #plt.imshow(GEO_slice)

        # construct black & white image
        GEO_bw1 = np.where(GEO_slice > np.median(GEO_slice) * 1.3, 1, 0)
        GEO_bw2 = np.where(GEO_slice < np.median(GEO_slice) * 0.7, 1, 0)
        GEO_bw = GEO_bw1 + GEO_bw2

        plt.imshow(GEO_bw)
        plt.title("GEO_slice")
        #plt.show()

        # find centers of geometric points
        labels, no_roi = meas.label(GEO_bw)
        print (no_roi)

        #if no_roi != 4:
        #    raise ValueError, "Unable to locate the geometric nodes. May need to adjust the geometric ROI slice. "


        geo_CofM = meas.center_of_mass(GEO_bw, labels, index=[1,2,3,4])

        self.geo_CofM = geo_CofM

        # distance calculations (result is in mm; nominal distance is 50mm (5cm)
        self.results['v1'] = ptg.dist_2points(geo_CofM[0], geo_CofM[1]) * self.improps['mm/Pixel']
        self.results_diffs['v1'] = self.results['v1'] - nominal_GEO
        self.results['h1'] = ptg.dist_2points(geo_CofM[0], geo_CofM[2]) * self.improps['mm/Pixel']
        self.results_diffs['h1'] = self.results['h1'] - nominal_GEO
        self.results['v2'] = ptg.dist_2points(geo_CofM[3], geo_CofM[1]) * self.improps['mm/Pixel']
        self.results_diffs['v2'] = self.results['v2'] - nominal_GEO
        self.results['h2'] = ptg.dist_2points(geo_CofM[3], geo_CofM[2]) * self.improps['mm/Pixel']
        self.results_diffs['h2'] = self.results['h2'] - nominal_GEO

        res_diffs = np.array((self.results_diffs['v1'], self.results_diffs['h1'], self.results_diffs['v2'],
                     self.results_diffs['h2']))

        if all(res_diffs < tolerances['GEO']) and all(res_diffs > -tolerances['GEO']):
            self.GEO_pass = True
        else:
            self.GEO_pass = False

        pass


    def find_UNIF(self):
        """
        Determine Uniformity from Uniformity slice. Averages 3 slices
        :return:
        """


        # mean from 3 slices around nominal slice
        UN_slice = np.mean(self.images[self.slices_of_interest['UN'] - 1:self.slices_of_interest['UN'] + 1,:,:], 0)

        for idx in range(5):
            # calculate ROI centers
            if idx == 0:
                self.UN_roi_bounds.append(self.slice_centers['UN'])
            else:
                self.UN_roi_bounds.append([self.slice_centers['HU'][0] + small_FOV_settings['UN_radius'] * np.cos(np.radians(self.phan_roll) + small_FOV_settings['UN_angles'][idx-1]),
                             self.slice_centers['HU'][1] + small_FOV_settings['UN_radius'] * np.sin(np.radians(self.phan_roll) + small_FOV_settings['UN_angles'][idx-1])])
            # sample the ROI
            UN_roi = UN_slice[
                     self.UN_roi_bounds[idx][0] - small_FOV_settings['UN_boxsize'] / 2:self.UN_roi_bounds[idx][0] + small_FOV_settings['UN_boxsize'] / 2,
                     self.UN_roi_bounds[idx][1] - small_FOV_settings['UN_boxsize'] / 2:self.UN_roi_bounds[idx][1] + small_FOV_settings['UN_boxsize'] / 2]
            # Calculate the mean HU value and the difference from nominal
            self.UN_val[idx] = np.mean(UN_roi * self.improps['Rescale Slope'] + self.improps['Rescale Intercept'])
            self.UN_diff[idx] = self.UN_val[idx]

        if all(self.UN_diff < tolerances['UN']) and all(self.UN_diff > -tolerances['UN']):
            self.UNIF_pass = True
        else:
            self.UNIF_pass = False


    def draw_ROIs(self,results, show=True):
        """
        Draw the ROI boxes for a given analysis (e.g. The boxes where HU was calculated from).
        :return:
        """

        fig, ((UN_ax, HU_ax), (SR_ax, LOCON_ax)) = plt.subplots(2,2)
        # Uniformity objects
        UN_ax.imshow(self.images[self.slices_of_interest['UN'],:,:])
        for idx in range(len(self.UN_roi_bounds)):
            UN_ax.add_patch(Rectangle((self.UN_roi_bounds[idx][1]-small_FOV_settings['UN_boxsize']/2,
                                       self.UN_roi_bounds[idx][0]-small_FOV_settings['UN_boxsize']/2),
                                      small_FOV_settings['UN_boxsize'], small_FOV_settings['UN_boxsize'],
                                      fill=False, edgecolor='black'))
        UN_ax.autoscale(tight=True)
        UN_ax.set_title('Uniformity Slice')


        # HU objects
        HU_ax.imshow(self.images[self.slices_of_interest['HU'],:,:])
        for idx in range(len(self.HU_roi_bounds)):
            HU_ax.add_patch(Rectangle((self.HU_roi_bounds[idx][1] - small_FOV_settings['HU_boxsize'] / 2,
                                       self.HU_roi_bounds[idx][0] - small_FOV_settings['HU_boxsize'] / 2),
                                      small_FOV_settings['HU_boxsize'], small_FOV_settings['HU_boxsize'],
                                      fill=False, edgecolor='black'))
        HU_ax.autoscale(tight=True)
        HU_ax.set_title('HU & Geometric Slice')

        # GEO objects
        for idx in [0,3]: # for 2 corner geo points...
            for idx2 in [1,2]:  # connect to the other 2 geo points
                HU_ax.plot([self.geo_CofM[idx][1] + small_FOV_settings['GEO_box_bounds'][1],
                           self.geo_CofM[idx2][1] + small_FOV_settings['GEO_box_bounds'][0]],
                           [self.geo_CofM[idx][0] + small_FOV_settings['GEO_box_bounds'][1],
                           self.geo_CofM[idx2][0] + small_FOV_settings['GEO_box_bounds'][0]],'black')



        # SR objects
        SR_ax.imshow(self.images[self.slices_of_interest['SR'],:,:])
        for pnt in self.MTF['ROI points']:
            SR_ax.plot(pnt[1],pnt[0])
        SR_ax.autoscale(tight=True)
        SR_ax.set_title('Spatial Resolution Slice')

        #TODO: Low contrast
        '''
        LOCON_ax.imshow(self.images[:,:,self.slices_of_interest['LOCON']])
        LOCON_ax.set_title('Low Contrast (In Development)')
        '''

        fig.savefig('Overview')
        results.addObject('Overview','Overview.png')

    #if show:
    #    plt.show()

    def return_results(self, results, using_pyqa=False):
        """
        """
        #TODO: make a bit prettier
        print('HU Regions: ', self.results['HU'])
        for i in range(1,7):
            results.addFloat('HU Regions %s'%i, self.results['HU'][i])

        print('HU Passed?: ', self.HU_pass)
        if self.HU_pass == True:
            results.addChar('HU Passed', 'True')
        else:
            results.addChar('HU Passed', 'False')

        print('Uniformity: ', self.UN_val)
        for i in range(len(self.UN_val)):
            results.addFloat('Uniformity %s'%i, self.UN_val[i])

        print('Uniformity Passed?: ', self.UNIF_pass)
        if self.UNIF_pass == True:
            results.addChar('Uniformity Passed', 'True')
        else:
            results.addChar('Uniformity Passed', 'False')
                
        print('MTF 50% (lp/mm): ', self.MTF['MTF_50'])
        results.addFloat('MTF 50% (lp/mm): ', self.MTF['MTF_50'])


        results.addFloat('Geometric distances v1',self.results['v1'])
        results.addFloat('Geometric distances v2',self.results['v2'])
        results.addFloat('Geometric distances h1',self.results['h1'])
        results.addFloat('Geometric distances h2',self.results['h2'])
        
        print('Geometry Passed?: ', self.GEO_pass)
        if self.GEO_pass == True:
            results.addChar('Geometry Passed','True')
        else:
            results.addChar('Geometry Passed','False')

    def analyze(self,params):
        """
        One-method analysis of catphan files.
        :return:
        """
        self.find_phan_centers(params)
        self.find_roll()

        self.find_HU()
        self.find_GEO()
        self.find_SR()
        self.find_UNIF()


def catphan500(data,results,**kwargs):


     paramdict = kwargs.get('params', None) #read out all the parameter tags from the config_xml



     catphan = CATPHAN()
     filelist = data.getAllInstances()

     catphan._load_files(filelist)
     catphan.analyze(paramdict)
     catphan.return_results(results)
     catphan.draw_ROIs(results)
    



