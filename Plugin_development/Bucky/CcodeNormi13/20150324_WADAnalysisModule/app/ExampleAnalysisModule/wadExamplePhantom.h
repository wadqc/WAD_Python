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

#ifndef WADEXAMPLEPHANTOM_H
#define WADEXAMPLEPHANTOM_H

// Custom includes
#include "wadPhantom.h"

// ITK includes
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkStatisticsImageFilter.h>

namespace wad
{
/**
 * \class ExamplePhantom
 * \brief Derived class implementing a example phantom
 */
	class ExamplePhantom : public Phantom
	{
	public:

    private:
	    typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractFilterType;
		typedef itk::Image< unsigned char, 2 > PngImageType;
		typedef itk::ImageFileWriter< PngImageType > PngWriterType;
		typedef itk::RescaleIntensityImageFilter< ImageType, PngImageType > RescalePngFilterType;
	    typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;

	    /// Compute mean at image centre in a square of size roiSize
        void ComputeMeanAtImageCentre(const double roiSize);
        
        /// Compute statistics at ROI
        void ComputeRoiStatistics();
		
		void DoAnalyze();

		std::string DoGetResult( const std::string quantity, const std::string unit, const std::string description );

		void DoInitialize( const ImageType *image );

        void DoWriteResultImage();

        /// Set the ROI
        void SetRoi( const IndexType index, const SizeType size );
        
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
