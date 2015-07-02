// C++ library for WAD image quality control (IQC) analysis modules.
// For use with the NVKF WAD IQC framework for automatic analysis of
// DICOM objects.
//
// Copyright (C) 2013, 2014  E.J. Rijkhorst / e.rijkhorst@vumc.nl
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef WADRADIOLOGYPHANTOM_H
#define WADRADIOLOGYPHANTOM_H

// Custom includes
#include "wadPhantom.h"

// ITK includes
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkStatisticsImageFilter.h>

namespace wad
{
/**
 * \class RadiologyPhantom
 * \brief Derived class implementing a radiology phantom
 *
 * \note Only 2D radiology phantoms are currently implemented.
 */
	class RadiologyPhantom : public Phantom
	{
	public:

    private:
	    typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractFilterType;
		typedef itk::Image< unsigned char, 2 > PngImageType;
		typedef itk::ImageDuplicator< ImageType > DuplicatorType;
		typedef itk::ImageFileWriter< PngImageType > PngWriterType;
		typedef itk::ImageRegionConstIterator< ImageType > RegionConstIteratorType;
		typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;
		typedef itk::InvertIntensityImageFilter < ImageType > InvertIntensityFilterType;
		typedef itk::MinimumMaximumImageCalculator <ImageType> MinMaxCalculatorType;
		typedef itk::RescaleIntensityImageFilter< ImageType, PngImageType > RescalePngFilterType;
	    typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;

	    /// Compute mean at image centre in a square of size roiSize
        void ComputeMeanAtImageCentre(const double roiSize);

        /// Compute minimum and maximum image intensity
        void ComputeMinMaxIntensity();
        
        /// Compute statistics at ROI
        void ComputeRoiStatistics();
        
        /// Copy input image to result image
        void CopyImage();

		/// Create an inverse of the input image
		void CreateInverseImage();
		
		void DoAnalyze();

		std::string DoGetResult( const std::string quantity, const std::string unit, const std::string description );

		void DoInitialize( const ImageType *image );

        void DoWriteResultImage();

        /// Set the ROI
        void SetRoi( const IndexType index, const SizeType size );
        
        /// Show the ROI in the result image
        void ShowRoi();

		/// Inverse of input image
		ImageType::Pointer m_inverseImage;

        /// Maximum intensity
        PixelType m_max;

        /// Minimum intensity
        PixelType m_min;
		
        /// Region of interest
        RegionType m_roi;

        /// Maximum of roi
        PixelType m_roiMax;
        
        /// Mean of roi
        PixelType m_roiMean;

        /// Minimum of roi
        PixelType m_roiMin;

        /// Standard deviation of roi
        PixelType m_roiSigma;

	}; // end of class
    
} // end of namespace

#endif
