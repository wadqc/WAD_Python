#!/usr/bin/python
# -*- coding: latin-1 -*-

import os

import numpy as np


from pylinac import CBCT

__version__='20160229'


class WADCatphan(CBCT):
    """Inherit CBCT class from pylinac and add some methods to evaluate it with WAD QC software"""

    def save_numeric_results_to_wad(self, results):
        """Save numeric results"""

        # HU
        for name in self.hu.roi_names:
            region = self.hu.rois[name]
            results.addFloat('HU region %s measured' % name, region.pixel_value, level=1)
            results.addFloat('std HU region %s measured' % name, region.std, level=2)

        # Uniformity
        univalues = []
        for name in self.uniformity.rois:
            region = self.uniformity.rois[name]
            results.addFloat('Uniformity region %s measured' % name, region.pixel_value, level=1)
            results.addFloat('std Uniformity region %s measured' % name, region.std, level=2)
            univalues.append(region.pixel_value)

        results.addFloat('Uniformity region average', np.average(univalues), level=1)
        results.addFloat('std uniformity region average', np.nanstd(univalues), level=2)

        # MTF
        results.addFloat('MTF 80% (lp/mm)', self.spatialres.mtf(80), level=1)

        # Geometric distances
        for key in self.geometry.lines.keys():
            line = self.geometry.lines[key]
            results.addFloat('Line length %sl' % key, line.length_mm, level=1)
        results.addFloat('Geometric Line Average (mm)', self.geometry.avg_line_length, level=1)

        # Slice thickness
        results.addFloat('Slice thickness', self.thickness.avg_slice_thickness, level=1)

        # Low contrast
        results.addFloat('Low contrast', self.lowcontrast.rois_visible, level=1)

    def save_images_to_wad(self, results):
        """Save figures/images"""
        self.save_analyzed_subimage('hu.png', subimage='hu', show=False)
        results.addObject('HU linearity image','hu.png')

        self.save_analyzed_subimage('un.png', subimage='un', show=False)
        results.addObject('HU uniformity image','un.png')

        self.save_analyzed_subimage('sp.png', subimage='sp', show=False)
        results.addObject('Spatial Resolution image','sp.png')

        self.save_analyzed_subimage('lc.png', subimage='lc', show=False)
        results.addObject('Low Contrast image','lc.png')

        self.save_analyzed_subimage('mtf.png', subimage='mtf', show=False)
        results.addObject('RMTF plot','mtf.png')

        self.save_analyzed_subimage('lin.png', subimage='lin', show=False)
        results.addObject('HU linearity values','lin.png')

        self.save_analyzed_subimage('prof.png', subimage='prof', show=False)
        results.addObject('HU uniformity profiles','prof.png')


def catphan504(data,results,**kwargs):
    dicomfolder = os.path.dirname(data.series_filelist[0][1])
    catphan = WADCatphan(dicomfolder)
    catphan.analyze()
    catphan.save_numeric_results_to_wad(results)
    catphan.save_images_to_wad(results)


