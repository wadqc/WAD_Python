#!/usr/bin/python
# -*- coding: latin-1 -*-

import os

import numpy as np


from pylinac import CBCT

__version__='22122015'


class WADCatphan(CBCT):
    """Inherit CBCT class from pylinac and add some methods to evaluate it with WAD QC software"""

    def save_numeric_results_to_wad(self, results):
        """Save numeric results"""

        # HU
        for name in self.hu.roi_names:
            region = self.hu.rois[name]
            results.addFloat('HU region %s measured' % name, region.pixel_value)
            results.addFloat('std HU region %s measured' % name, region.std)
            results.addFloat('HU region %s  nominal' % name, region.nominal_val)

        # Geometric distances
        for key in self.geometry.lines.keys():
            line = self.geometry.lines[key]
            results.addFloat('Line length %s (mm)' % key, line.length_mm)

        for key in self.geometry.diagonal_lines.keys():
            line = self.geometry.diagonal_lines[key]
            results.addFloat('Line length %s (mm)' % key, line.length_mm)

        results.addFloat('Geometric Line Average (mm)', self.geometry.avg_line_length)

        # Slice thickness
        results.addFloat('Slice thickness (mm)', self.thickness.avg_slice_thickness)

        # Uniformity
        if self.settings.uniformity_distance_mm is not None:
            univalues = []
            for name in self.uniformity.rois:
                region = self.uniformity.rois[name]
                results.addFloat('Uniformity region %s measured' % name, region.pixel_value)
                results.addFloat('std Uniformity region %s measured' % name, region.std)
                univalues.append(region.pixel_value)

            results.addFloat('Uniformity region average', np.average(univalues))
            results.addFloat('std uniformity region average', np.nanstd(univalues))

        # MTF
        if self.settings.spatial_resolution_distance_mm is not None:
            results.addFloat('MTF 80% (lp/mm)', self.spatialres.mtf(80))
            results.addFloat('MTF 50% (lp/mm)', self.spatialres.mtf(percent=50))
            results.addFloat('MTF 10% (lp/mm)', self.spatialres.mtf(percent=10))

        for i in range(len(self.spatialres.line_pair_frequency)):
            results.addFloat('MTF ' + str(self.spatialres.line_pair_frequency[i]) + ' (lp/mm)',
                             self.spatialres.line_pair_mtfs[i])

        #Expected number of peaks
        results.addFloat('Expected number of peaks', self.settings.nr_peaks_to_be_found)

        #Found number of peaks
        results.addFloat('Found number of peaks', len(self.spatialres.profile_peaks_idxs))

        #Expected number of valleys
        exp_valleys = self.settings.nr_peaks_to_be_found - self.settings.nr_line_pairs_to_analyze
        results.addFloat('Expected number of valleys', exp_valleys)

        #Found number of valleys
        results.addFloat('Found number of valleys', len(self.spatialres.profile_valleys_idxs))

        #Low contrast
        for name in self.lowcontrast.rois:
            results.addFloat('Contrast target ' + name + ' mm', self.lowcontrast.rois[name].contrast)
            results.addFloat('Contrast constant target ' + name + ' mm', self.lowcontrast.rois[name].contrast_constant)
            results.addFloat('Pixel value target ' + name + ' mm', self.lowcontrast.rois[name].pixel_value)
            results.addChar('Center target ' + name + ' mm', '(x=' + str(self.lowcontrast.rois[name].center.x) +
                            ', y=' + str(self.lowcontrast.rois[name].center.y) + ')')
            #nog meer toevoegen?

        for name in self.lowcontrast.bg_rois:
            results.addFloat('Pixel value ' + name, self.lowcontrast.bg_rois[name].pixel_value)
            results.addChar('Center target ' + name, '(x=' + str(self.lowcontrast.bg_rois[name].center.x) + ', y=' +
                            str(self.lowcontrast.bg_rois[name].center.y) + ')')
            #nog meer toevoegen?

    def save_images_to_wad(self, results):
        """Save figures/images"""
        self.save_analyzed_subimage('hu.png', subimage='hu', show=False)
        results.addObject('HU linearity image','hu.png')

        if self.settings.uniformity_distance_mm is not None:
            self.save_analyzed_subimage('un.png', subimage='un', show=False)
            results.addObject('HU uniformity image','un.png')

        if self.settings.spatial_resolution_distance_mm is not None:
            self.save_analyzed_subimage('sp.png', subimage='sp', show=False)
            results.addObject('Spatial Resolution image','sp.png')

        if self.settings.low_contrast_distance_mm is not None:
            self.save_analyzed_subimage('lc.png', subimage='lc', show=False)
            results.addObject('Low Contrast image','lc.png')

        self.save_analyzed_subimage('mtf.png', subimage='mtf', show=False)
        results.addObject('RMTF plot','mtf.png')

        self.save_analyzed_subimage('lin.png', subimage='lin', show=False)
        results.addObject('HU linearity values','lin.png')

        self.save_analyzed_subimage('prof.png', subimage='prof', show=False)
        results.addObject('HU uniformity profiles','prof.png')

        self.save_analyzed_subimage('CircleProfile.png', subimage='circle_profile', show=False)
        results.addObject('Original circle profile', 'CircleProfile.png')

        self.save_analyzed_subimage('SpacedProfileFoundPeaksValleys.png', subimage='spaced_profile_peaks_valleys',
                                    show=False)
        results.addObject('Spaced circle profile with the found peaks and valleys','SpacedProfileFoundPeaksValleys.png')


def analyze_catphan(data,results,**kwargs):

    paramdict = kwargs.get('params', None) #read out all the parameter tags from the config_xml

#is dit nodig voor allemaal?
    low_contrast_dist_mm = None
    angle_first_target_lc = None
    spatial_resolution_dist_mm = None
    nr_line_pairs_to_analyze = None
    start_angle_line_pairs_rad = None
    uniformity_dist_mm = None
    air_angle_deg = None
    pmp_angle_deg = None
    ldpe_angle_deg = None
    poly_angle_deg = None
    acrylic_angle_deg = None
    delrin_angle_deg = None
    teflon_angle_deg = None

    #params HU module
    if paramdict.find('hu_air_angle_deg') is not None:
       air_angle_deg = int(paramdict.find('hu_air_angle_deg').text)
    if paramdict.find('hu_pmp_angle_deg') is not None:
       pmp_angle_deg = int(paramdict.find('hu_pmp_angle_deg').text)
    if paramdict.find('hu_ldpe_angle_deg') is not None:
       ldpe_angle_deg = int(paramdict.find('hu_ldpe_angle_deg').text)
    if paramdict.find('hu_poly_angle_deg') is not None:
       poly_angle_deg = int(paramdict.find('hu_poly_angle_deg').text)
    if paramdict.find('hu_acrylic_angle_deg') is not None:
       acrylic_angle_deg = int(paramdict.find('hu_acrylic_angle_deg').text)
    if paramdict.find('hu_delrin_angle_deg') is not None:
       delrin_angle_deg = int(paramdict.find('hu_delrin_angle_deg').text)
    if paramdict.find('hu_teflon_angle_deg') is not None:
       teflon_angle_deg = int(paramdict.find('hu_teflon_angle_deg').text)

    #params Spatial Resolution module
    if paramdict.find('spatial_resolution_dist_mm') is not None:
        spatial_resolution_dist_mm = int(paramdict.find('spatial_resolution_dist_mm').text)
    if paramdict.find('nr_line_pairs_to_analyze') is not None:
        nr_line_pairs_to_analyze = int(paramdict.find('nr_line_pairs_to_analyze').text)
    if paramdict.find('start_angle_line_pairs_rad') is not None:
        start_angle_line_pairs_rad = float(paramdict.find('start_angle_line_pairs_rad').text)

    #params Low Contrast module
    if paramdict.find('low_contrast_dist_mm') is not None:
        low_contrast_dist_mm = int(paramdict.find('low_contrast_dist_mm').text)
    if paramdict.find('angle_first_target_low_contrast') is not None:
        angle_first_target_lc = float(paramdict.find('angle_first_target_low_contrast').text)

    #params Uniformity module
    if paramdict.find('uniformity_dist_mm') is not None:
        uniformity_dist_mm = int(paramdict.find('uniformity_dist_mm').text)

    dicomfolder = os.path.dirname(data.series_filelist[0][1])

    catphan = WADCatphan(folderpath=dicomfolder, lc_dist_mm=low_contrast_dist_mm, sr_dist_mm=spatial_resolution_dist_mm,
                         un_dist_mm=uniformity_dist_mm, nr_line_pairs_to_analyze=nr_line_pairs_to_analyze,
                         start_angle_line_pairs_rad=start_angle_line_pairs_rad, air_angle_deg=air_angle_deg,
                         pmp_angle_deg=pmp_angle_deg, ldpe_angle_deg=ldpe_angle_deg, poly_angle_deg=poly_angle_deg,
                         acrylic_angle_deg=acrylic_angle_deg, delrin_angle_deg=delrin_angle_deg,
                         teflon_angle_deg=teflon_angle_deg, lc_angle_first_target=angle_first_target_lc)
    catphan.analyze()
    catphan.save_numeric_results_to_wad(results)
    catphan.save_images_to_wad(results)
