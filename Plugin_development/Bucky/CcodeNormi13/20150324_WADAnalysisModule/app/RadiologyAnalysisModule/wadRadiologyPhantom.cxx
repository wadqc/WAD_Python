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

#include "wadRadiologyPhantom.h"

namespace wad
{
    void RadiologyPhantom::ComputeMeanAtImageCentre(const double roiSize)
    {
        if( m_debug ) std::cout << "RadiologyPhantom::ComputeMeanAtImageCentre" << std::endl;
        IndexType index;
        SizeType  size;
        index[0] = IndexValueType(m_imageSize[0]/2 - roiSize/2.0/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 - roiSize/2.0/m_imageSpacing[1]);
        size[0]  = SizeValueType(roiSize/m_imageSpacing[0]);
        size[1]  = SizeValueType(roiSize/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
    }
   
    void RadiologyPhantom::ComputeMinMaxIntensity()
    {
        if( m_debug ) std::cout << "RadiologyPhantom::ComputeMinMaxIntensity" << std::endl;
        MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();
        minMaxCalculator->SetImage( m_image );
        minMaxCalculator->Compute();
        m_min = minMaxCalculator->GetMinimum();
        m_max = minMaxCalculator->GetMaximum();
        if( m_debug )
        {
            std::cout << m_min << std::endl;
            std::cout << m_max << std::endl;
        }
    }

    void RadiologyPhantom::ComputeRoiStatistics()
    {
        if( m_debug ) std::cout << "RadiologyPhantom::ComputeRoiStatistics" << std::endl;
        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetExtractionRegion( m_roi );
        extractFilter->SetInput( m_image );
        StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
        statisticsFilter->SetInput( extractFilter->GetOutput() );
        try
        {
            statisticsFilter->Update();
        }
        catch (itk::ExceptionObject & e)
        {
            std::cerr << e << std::endl;
        }
        m_roiMean = statisticsFilter->GetMean();
        m_roiSigma = statisticsFilter->GetSigma();
        m_roiMin = statisticsFilter->GetMinimum();
        m_roiMax = statisticsFilter->GetMaximum();
        if( m_debug )
        {
            std::cout << m_roiMean << std::endl;
            std::cout << m_roiSigma << std::endl;
            std::cout << m_roiMin << std::endl;
            std::cout << m_roiMax << std::endl;
        }
    }

    void RadiologyPhantom::CreateInverseImage()
    {
        if( m_debug ) std::cout << "RadiologyPhantom::CreateInverseImage" << std::endl;
		// Compute the min and max intensity
		this->ComputeMinMaxIntensity();
		// Invert the image
		InvertIntensityFilterType::Pointer invertIntensityFilter = InvertIntensityFilterType::New();
		invertIntensityFilter->SetInput( m_image );
		invertIntensityFilter->SetMaximum( m_max );
		try
		{
			invertIntensityFilter->Update();
		}
		catch (itk::ExceptionObject & e)
		{
			std::cerr << e << std::endl;
			//return EXIT_FAILURE;
		}
		m_inverseImage = invertIntensityFilter->GetOutput();
	}

	void RadiologyPhantom::DoAnalyze()
	{
        if( m_debug ) std::cout << "RadiologyPhantom::DoAnalyze" << std::endl;
		this->ComputeMeanAtImageCentre(25.0);
		this->ShowRoi();
	}

    std::string RadiologyPhantom::DoGetResult( const std::string quantity, const std::string unit, const std::string description )
	{
        if( m_debug ) std::cout << "RadiologyPhantom::DoGetResult" << std::endl;
		// Results are identified according to quantity, unit and description
		std::string result;
        if ( quantity == ""        && unit == ""           && description == "Analysis status"           ) result = "Completed";
        if ( quantity == "mean"    && unit == "grey value" && description == "Mean at phantom center"    ) result = to_string( m_roiMean );
        if ( quantity == "sigma"   && unit == "grey value" && description == "Sigma at phantom center"   ) result = to_string( m_roiSigma );
        if ( quantity == "minimum" && unit == "grey value" && description == "Minimum at phantom center" ) result = to_string( m_roiMin );
        if ( quantity == "maximum" && unit == "grey value" && description == "Maximum at phantom center" ) result = to_string( m_roiMax );
		return result;
	}

    void RadiologyPhantom::DoInitialize( const ImageType *image )
    {
        if( m_debug ) std::cout << "RadiologyPhantom::DoInitialize" << std::endl;
		// Create an inverse of the input image
		this->CreateInverseImage();
	}
	
    void RadiologyPhantom::DoWriteResultImage()
    {
		if( m_debug ) std::cout << "RadiologyPhantom::DoWriteResultImage" << std::endl;
		std::string outputFilename = "result.png";
		// Rescale to char values
		RescalePngFilterType::Pointer filter = RescalePngFilterType::New();
		filter->SetOutputMinimum(   0 );
		filter->SetOutputMaximum( 255 );
		filter->SetInput( m_resultImage );

		// Write png image
		PngWriterType::Pointer writer = PngWriterType::New();
		writer->SetFileName(outputFilename);
		writer->SetInput(filter->GetOutput() );
        if( m_feedback )
        {
            std::cout << "Writing " << outputFilename << std::endl;
        }
		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject & e)
		{
			std::cerr << e << std::endl;
		}
	}

    void RadiologyPhantom::SetRoi( const IndexType index, const SizeType size )
    {
        if( m_debug ) std::cout << "RadiologyPhantom::SetRoi" << std::endl;
        m_roi.SetIndex( index );
        m_roi.SetSize( size );

        if( m_debug )
        {
            std::cout << "Index: " << m_roi.GetIndex() << std::endl;
            std::cout << "Size:  " << m_roi.GetSize() << std::endl;
        }
    }

    void RadiologyPhantom::ShowRoi()
    {
		if( m_debug ) std::cout << "RadiologyPhantom::ShowRoi" << std::endl;
		// Show roi in result image by inverse intensity
		RegionConstIteratorType iter1( m_inverseImage, m_roi );
		RegionIteratorType iter2( m_resultImage, m_roi );
		PixelType offset = 0.1*(m_max-m_min); // adding 10% of contrast range to make rois more visible in result image
		while(!iter2.IsAtEnd())
		{
			iter2.Set( iter1.Get() + offset );
			++iter1;
			++iter2;
		}

		// Shrink ROI
		IndexType index = m_roi.GetIndex();
		SizeType size = m_roi.GetSize();
		unsigned pix = 10;
		index[0] += pix;
		index[1] += pix;
		size[0] -= 2*pix;
		size[1] -= 2*pix;
        this->SetRoi( index, size );

		// Show shrunken roi in result image
		RegionConstIteratorType iter3( m_image, m_roi );
		RegionIteratorType iter4( m_resultImage, m_roi );
		while(!iter4.IsAtEnd())
		{
			
			iter4.Set( iter3.Get() );
			++iter3;
			++iter4;
		}
	}

} // end namespace
