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

#ifndef WADBUCKYPHANTOM_H
#define WADBUCKYPHANTOM_H

// Custom includes
#include "wadPhantom.h"

// ITK includes
#include <itkExtractImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionExclusionConstIteratorWithIndex.h>
#include <itkImageRegionExclusionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkMeanProjectionImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkStatisticsImageFilter.h>

#include <itkBilateralImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkZeroCrossingBasedEdgeDetectionImageFilter.h>

namespace wad
{
/**
 * \class BuckyPhantom
 * \brief Derived class implementing a radiology phantom
 *
 * Currently the NORMI13 phantom for Bucky X-ray images is implemented.
 */	
	class BuckyPhantom : public Phantom
	{
	public:
        /// Set the scale factor using the SID and thickness
		void SetScaleFactor( const double thickness, const double SID );

    protected:
        
    private:
        struct InhomogeneityStruct
        {
            // 0=top left, 1=top right, 2=centre, 3=bottom left, 4=bottom right
            double mean[5];
            double sigma[5];
            double min[5];
            double max[5];
        };

        typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractFilterType;
        typedef itk::Image< unsigned char,  2 > PngImageType;
        typedef itk::ImageDuplicator< ImageType > DuplicatorType;
        typedef itk::ImageFileWriter< ImageType > WriterType;
        typedef itk::ImageFileWriter< PngImageType > PngWriterType;
        typedef double Pixel1DType; // must be float or double!
        typedef itk::Image< Pixel1DType, 1 > Image1DType;
        typedef itk::ImageRegionConstIterator< Image1DType > RegionConstIterator1DType;	
        typedef itk::ImageRegionConstIterator< ImageType > RegionConstIteratorType;	
        typedef itk::ImageRegionExclusionConstIteratorWithIndex< ImageType > RegionExclusionConstIteratorType;	
        typedef itk::ImageRegionExclusionIteratorWithIndex< ImageType > RegionExclusionIteratorType;	
        typedef itk::ImageRegionIterator< Image1DType > RegionIterator1DType;	
        typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;	
        typedef itk::InvertIntensityImageFilter < ImageType > InvertIntensityFilterType;
        typedef itk::MeanProjectionImageFilter < ImageType, Image1DType > MeanProjectionFilterType;
        typedef itk::MinimumMaximumImageCalculator <Image1DType> MinMaxCalculator1DType;
        typedef itk::MinimumMaximumImageCalculator <ImageType> MinMaxCalculatorType;
        typedef itk::RescaleIntensityImageFilter< ImageType, PngImageType > RescalePngFilterType;
        typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;
        
        typedef itk::BilateralImageFilter< ImageType, ImageType > BilateralFilterType;
        
        typedef itk::Image<double, Dimension> RealImageType;
        typedef itk::CastImageFilter< ImageType, RealImageType> CastToRealFilterType;
        typedef itk::RescaleIntensityImageFilter<RealImageType, ImageType > RescaleFilter;
        typedef itk::CannyEdgeDetectionImageFilter<RealImageType, RealImageType> CannyEdgeDetectionFilter;
        typedef itk::SobelEdgeDetectionImageFilter <RealImageType, RealImageType> SobelEdgeDetectionFilterType;
        typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter< RealImageType, RealImageType> ZeroCrossingBasedEdgeDetectionFilterType;
        
        // Orientation constants
        enum ORIENTATION {UNKNOWN, NORTH, EAST, SOUTH, WEST};

		/// Check if the ROI contains the high contrast pattern
        int CheckForHighContrastPattern( ORIENTATION assumedOrientation );

        /// Compute inhomogeneity at centre and four corners of phantom
		void ComputeInhomogeneity();
        
        /// Compute mean at image centre
        void ComputeMeanAtImageCentre();

        /// Compute minimum and maximum image intensity
        void ComputeMinMaxIntensity();
        
        /// Compute statistics at ROI
        void ComputeRoiStatistics();
        
        /// Copy input image to result image
        void CopyImage();

		/// Create an inverse of the input image
		void CreateInverseImage();
	
        void DoAnalyze();

        /// ITK approach for analysis
        void AnalyzeITK();

		std::string DoGetResult( const std::string quantity, const std::string unit, const std::string description );

		void DoInitialize( const ImageType *image );

        void DoWriteResultImage();

        /// Find inner edges of diaphragm by searching along horizontal and vertical centre lines
        void FindDiaphragmEdges();

        /// Find the orientation of the phantom in the image
        void FindOrientation();

        /// Set the ROI
        void SetRoi( const IndexType index, const SizeType size );
        
        /// Set the phantom thickness (in mm)
	    void SetThickness( const double thickness );

        /// Show the diaphragm edges on the result image
        void ShowDiaphragmEdges();

        /// Show the ROI in the result image
        void ShowRoi();

        /// Diaphragm region of interest
        RegionType m_diaphragmRoi;

        /// Physical diaphragm centre in mm
		double m_diaphragmCentre[2];

		/// Physical diaphragm size in mm
		double m_diaphragmSize[2];
		
		/// Inhomogeneity
		InhomogeneityStruct m_inhomogeneity;

		/// Inverse of input image
		ImageType::Pointer m_inverseImage;

        /// Maximum intensity
        PixelType m_max;

        /// Minimum intensity
        PixelType m_min;
        
		/// Orientation of phantom in image
		ORIENTATION m_orientation;
		
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

		/// Scale factor
		double m_scaleFactor;
		
		/// Thickness in mm
		double m_thickness;
        
	}; // end of class
    
} // end of namespace

#endif
