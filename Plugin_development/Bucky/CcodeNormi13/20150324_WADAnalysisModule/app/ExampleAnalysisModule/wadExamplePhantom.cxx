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

#include "wadExamplePhantom.h"

namespace wad
{
    void ExamplePhantom::ComputeMeanAtImageCentre(const double roiSize)
    {
        if( m_debug ) std::cout << "ExamplePhantom::ComputeMeanAtImageCentre" << std::endl;
        IndexType index;
        SizeType  size;
        index[0] = IndexValueType(m_imageSize[0]/2 - roiSize/2.0/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 - roiSize/2.0/m_imageSpacing[1]);
        size[0]  = SizeValueType(roiSize/m_imageSpacing[0]);
        size[1]  = SizeValueType(roiSize/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
    }

    void ExamplePhantom::ComputeRoiStatistics()
    {
        if( m_debug ) std::cout << "ExamplePhantom::ComputeRoiStatistics" << std::endl;
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

	void ExamplePhantom::DoAnalyze()
	{
        if( m_debug ) std::cout << "ExamplePhantom::DoAnalyze" << std::endl;
		this->ComputeMeanAtImageCentre(25.0);
	}

    std::string ExamplePhantom::DoGetResult( const std::string quantity, const std::string unit, const std::string description )
	{
        if( m_debug ) std::cout << "ExamplePhantom::DoGetResult" << std::endl;
		// Results are identified according to quantity, unit and description
		std::string result;
        if ( quantity == ""        && unit == ""           && description == "Analysis status"           ) result = "Completed";
        if ( quantity == "mean"    && unit == "grey value" && description == "Mean at phantom center"    ) result = to_string( m_roiMean );
		return result;
	}

    void ExamplePhantom::DoInitialize( const ImageType *image )
    {
        if( m_debug ) std::cout << "ExamplePhantom::DoInitialize" << std::endl;
		// Create an inverse of the input image
//		this->CreateInverseImage();
	}
	
    void ExamplePhantom::DoWriteResultImage()
    {
		if( m_debug ) std::cout << "ExamplePhantom::DoWriteResultImage" << std::endl;
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

    void ExamplePhantom::SetRoi( const IndexType index, const SizeType size )
    {
        if( m_debug ) std::cout << "ExamplePhantom::SetRoi" << std::endl;
        m_roi.SetIndex( index );
        m_roi.SetSize( size );

        if( m_debug )
        {
            std::cout << "Index: " << m_roi.GetIndex() << std::endl;
            std::cout << "Size:  " << m_roi.GetSize() << std::endl;
        }
    }

} // end namespace
