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

#include "wadBuckyPhantom.h"

namespace wad
{	
	int BuckyPhantom::CheckForHighContrastPattern( ORIENTATION assumedOrientation )
	{
        if(m_debug) std::cout << "CheckForHighContrastPattern" << std::endl;

		// Set the return value to zero to start with (i.e. the ROI does not
		// contain the high contrast pattern)
		int returnValue = 0;

		// Extract the current roi and calculate the mean by projecting
		// along the dimension with smallest size
        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetExtractionRegion( m_roi );
        extractFilter->SetInput( m_image );
		MeanProjectionFilterType::Pointer meanProjectionFilter = MeanProjectionFilterType::New();
		meanProjectionFilter->SetInput( extractFilter->GetOutput() );
		switch( assumedOrientation )
		{
		case NORTH:
		case SOUTH:
			meanProjectionFilter->SetProjectionDimension(1);
		break;
		case EAST:
		case WEST:
			meanProjectionFilter->SetProjectionDimension(0);
		break;
		}
        try
        {
            meanProjectionFilter->Update();
        }
        catch (itk::ExceptionObject & e)
        {
            std::cerr << e << std::endl;
        }
		Image1DType::Pointer image1D = meanProjectionFilter->GetOutput();
		
		// Calculate the min and max value of the projected image
		MinMaxCalculator1DType::Pointer minMaxCalculator = MinMaxCalculator1DType::New();
		minMaxCalculator->SetImage( image1D );
        minMaxCalculator->Compute();
        Pixel1DType min = minMaxCalculator->GetMinimum();
        Pixel1DType max = minMaxCalculator->GetMaximum();

		// Calculate the step values
		Pixel1DType stepFact[7];
		// NORMI13 TERGOOI
         stepFact[0] = 0.00;
         stepFact[1] = 0.24;
         stepFact[2] = 0.42;
         stepFact[3] = 0.59;
         stepFact[4] = 0.75;
         stepFact[5] = 0.90;
         stepFact[6] = 1.00;
		
		// DIGI13 MCL
//		stepFact[0] = 0.00;
//		stepFact[1] = 0.14;
//		stepFact[2] = 0.34;
//		stepFact[3] = 0.60;
//		stepFact[4] = 0.81;
//		stepFact[5] = 0.93;
//		stepFact[6] = 1.00;
		
		Pixel1DType stepVal[7];
		for(unsigned i=0; i<7; ++i)
		{
			stepVal[i] = min + stepFact[i] * (max - min);
			if( m_debug) std::cout << "Step value " << i << ": " << stepVal[i] << std::endl;
		}
		
//DEBUG		
{
		if(assumedOrientation==NORTH)
		{
			std::ofstream myfile;
			myfile.open ("profile1.dat");
			RegionConstIterator1DType iter(image1D, image1D->GetLargestPossibleRegion());
			myfile << "# profile data" << std::endl;
			while( !iter.IsAtEnd() )
			{
				myfile << iter.Get() << std::endl;
				++iter;
			}
			myfile.close();
		}
}

		// Segment into 7 intensity levels
		Pixel1DType minLevel;
		Pixel1DType maxLevel;
		RegionIterator1DType iter(image1D, image1D->GetLargestPossibleRegion());
		iter.GoToBegin();
		while( !iter.IsAtEnd() )
		{
			for( unsigned i=0; i<7; ++i)
			{
				Pixel1DType currentStep;
				if(i==0)
				{
					currentStep = stepVal[1] - stepVal[0];
				}
				else
				{
					currentStep = stepVal[i] - stepVal[i-1];
				}
				minLevel = stepVal[i] - 0.5 * currentStep;
				maxLevel = stepVal[i] + 0.5 * currentStep;
				Pixel1DType val = iter.Get();
				if( val >= minLevel && val <= maxLevel) iter.Set( i ); 
			}
			++iter;
		}
		// Set the expected step sign according to the assumed orientation
		Pixel1DType expectedStepSign;
		switch( assumedOrientation )
		{
		case NORTH:
		case EAST:
			expectedStepSign = -1.0;
		break;
		case SOUTH:
		case WEST:
			expectedStepSign = 1.0;
		break;
		}
		// Count the number of intensity steps and compute the average step length.
		// Note that the total step length of 5 steps is computed (the first and
		// last one are ignored since there the derivative is not in the direction
		// of the expected step sign.
		iter.GoToBegin();
		Pixel1DType val1 = iter.Get();
		++iter;
		unsigned nrOfSteps = 0;
		double totalStepLength = 0;
		unsigned step = 0;
		while( !iter.IsAtEnd() )
		{
			Pixel1DType val2 = iter.Get();
			// Compare derivative to the expected step sign
			if( (val2 - val1) == expectedStepSign )
			{
				++nrOfSteps;
				if(nrOfSteps==1) totalStepLength = step;
				if(nrOfSteps==6) totalStepLength = step - totalStepLength;
			}
			val1 = val2;
			++step;
			++iter;
		}
		totalStepLength *= m_imageSpacing[0];
		totalStepLength /= m_scaleFactor;
		
		if( m_debug)
		{
			std::cout << "Nr of steps found " << nrOfSteps << std::endl;
			std::cout << "Total step length " << totalStepLength << std::endl;
		}
		
		// If this is the high contrast object, the number of steps should be 6
		// and the total step length should be close to 5 * 20 = 100 mm
		// and the total step length should be close to 5 * 18 = 90 mm
		if( nrOfSteps==6 )
		{
			if( totalStepLength > 90.0 && totalStepLength < 110.0 ) // NORMI13 TERGOOI
//			if( totalStepLength > 80.0 && totalStepLength < 100.0 ) // DIGI13 MCL
			{
				returnValue = 1;
			}
		}

//DEBUG		
{
		if(assumedOrientation==NORTH)
		{
			std::ofstream myfile;
			myfile.open ("profile2.dat");
			RegionConstIterator1DType iter(image1D, image1D->GetLargestPossibleRegion());
			myfile << "# profile data" << std::endl;
			while( !iter.IsAtEnd() )
			{
				myfile << iter.Get() << std::endl;
				++iter;
			}
			myfile.close();
		}
}

		return returnValue;
	}
	
    void BuckyPhantom::ComputeInhomogeneity()
    {
        if(m_debug) std::cout << "ComputeInhomogeneity" << std::endl;

		// The centres of the four corner squares for inhomogeneity measurements
		// are located (105, 105) mm away from the phantom's centre.
		// The square size is 30x30 mm, the roi size is 10x10 mm

        IndexType index;
        SizeType  size;
		unsigned i=0;
		
		// Location and size of high contrast pattern in image
		double l[2]; // top left corner of ROI
		double s[2]; // size of ROI
		double roiSize = 10.0; // ROI is 10 mm square
		s[0] = roiSize * m_scaleFactor;
		s[1] = roiSize * m_scaleFactor;
        size[0]  = SizeValueType(s[0]/m_imageSpacing[0]);
        size[1]  = SizeValueType(s[1]/m_imageSpacing[1]);

		std::cout << "s[0]: " << s[0] << std::endl;
		std::cout << "s[1]: " << s[1] << std::endl;
		std::cout << "size[0]: " << size[0] << std::endl;
		std::cout << "size[1]: " << size[1] << std::endl;
		
		// Top left
		l[0] = (-105.0 - roiSize/2.0) * m_scaleFactor;
		l[1] = (-105.0 - roiSize/2.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		this->ShowRoi();
		m_inhomogeneity.mean[i]  = m_roiMean;
		m_inhomogeneity.sigma[i] = m_roiSigma;
		m_inhomogeneity.min[i]   = m_roiMin;
		m_inhomogeneity.max[i]   = m_roiMax;
		++i;
		
		// Top right
		l[0] = (+105.0 - roiSize/2.0) * m_scaleFactor;
		l[1] = (-105.0 - roiSize/2.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		this->ShowRoi();
		m_inhomogeneity.mean[i]  = m_roiMean;
		m_inhomogeneity.sigma[i] = m_roiSigma;
		m_inhomogeneity.min[i]   = m_roiMin;
		m_inhomogeneity.max[i]   = m_roiMax;
		++i;

		// Centre
		l[0] = ( 0.0 - roiSize/2.0) * m_scaleFactor;
		l[1] = ( 0.0 - roiSize/2.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		this->ShowRoi();
		m_inhomogeneity.mean[i]  = m_roiMean;
		m_inhomogeneity.sigma[i] = m_roiSigma;
		m_inhomogeneity.min[i]   = m_roiMin;
		m_inhomogeneity.max[i]   = m_roiMax;
		++i;

		// Bottom left
		l[0] = (-105.0 - roiSize/2.0) * m_scaleFactor;
		l[1] = (+105.0 - roiSize/2.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		this->ShowRoi();
		m_inhomogeneity.mean[i]  = m_roiMean;
		m_inhomogeneity.sigma[i] = m_roiSigma;
		m_inhomogeneity.min[i]   = m_roiMin;
		m_inhomogeneity.max[i]   = m_roiMax;
		++i;

		// Bottom right
		l[0] = (+105.0 - roiSize/2.0) * m_scaleFactor;
		l[1] = (+105.0 - roiSize/2.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		this->ShowRoi();
		m_inhomogeneity.mean[i]  = m_roiMean;
		m_inhomogeneity.sigma[i] = m_roiSigma;
		m_inhomogeneity.min[i]   = m_roiMin;
		m_inhomogeneity.max[i]   = m_roiMax;
		++i;
	}

    void BuckyPhantom::ComputeMeanAtImageCentre()
    {
        if( m_debug ) std::cout << "ComputeMeanAtImageCentre" << std::endl;
        // Assuming the diaphragm is set on the 26 cm grid line
        // Compute mean intensity at centre of image in a 15x15 mm^2 roiSize
        IndexType index;
        SizeType  size;
		double roiSize = 15.0 * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 - roiSize/2.0/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 - roiSize/2.0/m_imageSpacing[1]);
        size[0]  = SizeValueType(roiSize/m_imageSpacing[0]);
        size[1]  = SizeValueType(roiSize/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
    }

    void BuckyPhantom::ComputeMinMaxIntensity()
    {
        if( m_debug ) std::cout << "ComputeMinMaxIntensity" << std::endl;
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
    
    void BuckyPhantom::ComputeRoiStatistics()
    {
        if( m_debug ) std::cout << "ComputeRoiStatistics" << std::endl;
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

    void BuckyPhantom::CopyImage()
	{
		if( m_debug ) std::cout << "CopyImage" << std::endl;
		DuplicatorType::Pointer duplicator = DuplicatorType::New();
		duplicator->SetInputImage(m_image);
		duplicator->Update();
		m_resultImage = duplicator->GetOutput();
	}

    void BuckyPhantom::CreateInverseImage()
    {
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

    void BuckyPhantom::DoAnalyze()
	{
        if( m_debug ) std::cout << "BuckyPhantom::DoAnalyze" << std::endl;

//        this->AnalyzeITK();
//        return;
        
        // Find and show the diaphragm edges
        this->FindDiaphragmEdges();
        this->ShowRoi();
        this->ShowDiaphragmEdges();
        
        // Find the orientation of the phantom
        this->FindOrientation();
        
        // Compute mean and standard deviation at centre and four corners of image
        this->ComputeInhomogeneity();
    }
        
    void BuckyPhantom::AnalyzeITK()
	{
        if( m_debug ) std::cout << "BuckyPhantom::DoAnalyzeITK" << std::endl;
        // 20140209: TESTING ITK APPROACH
        // 1) Edge preserving smoothing: http://www.itk.org/Doxygen/html/classitk_1_1BilateralImageFilter.html
        if( m_feedback ) std::cout << "  1) bilateral smoothing" << std::endl;
        double rangeSigma = 10.0;
        double domainSigma = 0.4;
        BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
        bilateralFilter->SetInput( m_image );
        bilateralFilter->SetDomainSigma(domainSigma);
        bilateralFilter->SetRangeSigma(rangeSigma);
//		try
//		{
//			bilateralFilter->Update();
//		}
//		catch (itk::ExceptionObject & e)
//		{
//			std::cerr << e << std::endl;
//			//return EXIT_FAILURE;
//		}
        
        // 2) Edge detection: http://www.itk.org/Doxygen/html/classitk_1_1CannyEdgeDetectionImageFilter.html
        if( m_feedback ) std::cout << "  2) edge detection" << std::endl;
        double variance = 5.0;
        double upperThreshold = 0.0;
        double lowerThreshold = 0.0;
        CastToRealFilterType::Pointer castToRealFilter = CastToRealFilterType::New();
//        castToRealFilter->SetInput( bilateralFilter->GetOutput() );
        castToRealFilter->SetInput( m_image );
        
        CannyEdgeDetectionFilter::Pointer cannyFilter = CannyEdgeDetectionFilter::New();
        cannyFilter->SetInput( castToRealFilter->GetOutput() );
        cannyFilter->SetVariance( variance );
        cannyFilter->SetUpperThreshold( upperThreshold );
        cannyFilter->SetLowerThreshold( lowerThreshold );
        
        SobelEdgeDetectionFilterType::Pointer sobelFilter = SobelEdgeDetectionFilterType::New();
        sobelFilter->SetInput( castToRealFilter->GetOutput() );

        ZeroCrossingBasedEdgeDetectionFilterType::Pointer zeroEdgeFilter = ZeroCrossingBasedEdgeDetectionFilterType::New();
        zeroEdgeFilter->SetInput( castToRealFilter->GetOutput() );
        zeroEdgeFilter->SetVariance( variance );
        //zeroEdgeFilter->SetMaximumError( 0.5 );
        
        RescaleFilter::Pointer rescaleFilter = RescaleFilter::New();
        rescaleFilter->SetOutputMinimum( 0 );
        rescaleFilter->SetOutputMaximum( 255 );
//        rescaleFilter->SetInput( cannyFilter->GetOutput() );
//        rescaleFilter->SetInput( sobelFilter->GetOutput() );
        rescaleFilter->SetInput( zeroEdgeFilter->GetOutput() );
		try
		{
			rescaleFilter->Update();
		}
		catch (itk::ExceptionObject & e)
		{
			std::cerr << e << std::endl;
			//return EXIT_FAILURE;
		}

        // 3) Connected components
        // 4) labeling statistics (diameter?)

        m_resultImage = rescaleFilter->GetOutput();
	}

    std::string BuckyPhantom::DoGetResult( const std::string quantity, const std::string unit, const std::string description )
	{
        if( m_debug ) std::cout << "BuckyPhantom::DoGetResult" << std::endl;

		// Results are identified according to quantity, unit and description
		std::string result;
        
        if ( quantity == "length" && unit == "mm" && description == "diaphragm center x" ) result = to_string( m_diaphragmCentre[0] );
        if ( quantity == "length" && unit == "mm" && description == "diaphragm center y" ) result = to_string( m_diaphragmCentre[1] );
        if ( quantity == "length" && unit == "mm" && description == "diaphragm size x"   ) result = to_string( m_diaphragmSize[0] );
        if ( quantity == "length" && unit == "mm" && description == "diaphragm size y"   ) result = to_string( m_diaphragmSize[1] );

        if ( quantity == "" && unit == "" && description == "orientation" )
        {
            switch( m_orientation )
            {
                case NORTH:
                    result = "north";
                    break;
                case EAST:
                    result = "east";
                    break;
                case SOUTH:
                    result = "south";
                    break;
                case WEST:
                    result = "west";
                    break;
                default:
                    result = "UNKNOWN";
            }
        }
        
        double aver = 0.0;
        for(int i=0;i<5;++i) aver += m_inhomogeneity.mean[i];
        aver /= 5.0;
        if ( quantity == "" && unit == "%" && description == "inhomogeneity top left"     ) result = to_string( 100.0*(m_inhomogeneity.mean[0]-aver) / aver );
        if ( quantity == "" && unit == "%" && description == "inhomogeneity top right"    ) result = to_string( 100.0*(m_inhomogeneity.mean[1]-aver) / aver );
        if ( quantity == "" && unit == "%" && description == "inhomogeneity center"       ) result = to_string( 100.0*(m_inhomogeneity.mean[2]-aver) / aver );
        if ( quantity == "" && unit == "%" && description == "inhomogeneity bottom left"  ) result = to_string( 100.0*(m_inhomogeneity.mean[3]-aver) / aver );
        if ( quantity == "" && unit == "%" && description == "inhomogeneity bottom right" ) result = to_string( 100.0*(m_inhomogeneity.mean[4]-aver) / aver );
                                                                                                                    
        return result;
	}
	
    void BuckyPhantom::DoInitialize( const ImageType *image )
    {
        if( m_debug ) std::cout << "BuckyPhantom::DoInitialize" << std::endl;
		// Create an inverse of the input image
		this->CreateInverseImage();
		
		// Set orientation to unknown
		m_orientation = UNKNOWN;
	}
    
    void BuckyPhantom::DoWriteResultImage()
    {
		if( m_debug ) std::cout << "BuckyPhantom::DoWriteResultImage" << std::endl;
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
    
    void BuckyPhantom::FindDiaphragmEdges()
    {
        if(m_debug) std::cout << "FindDiaphragmEdges" << std::endl;
        // First compute the mean it the image centre
        this->ComputeMeanAtImageCentre();
        
        // Set the edgeValue to 1.25 times the mean value at the image centre
        PixelType edgeValue = 1.25 * m_roiMean; // NORMI13 TERGOOI
//        PixelType edgeValue = 1.3 * m_roiMean; // DIGI13 MCL
        if(m_debug)
        {
            std::cout << "Finding edges of diaphragm using intensity value of " << edgeValue << std::endl;
        }
        
        // Find indices of diaphragm edges
        IndexType currentIndex;
		// 45 mm + offset from centre to avoid dose detector
		// Note: avoiding dose detector in horizontal as well as vertical direction
		//       since we don't know where it is located...
        double offset = 45.0 * m_scaleFactor; // NORMI13 TERGOOI
//		double offset = 0.0 * m_scaleFactor; // DIGI13 MCL
        // left
        IndexValueType indexLeft = m_imageIndex[0] + IndexValueType(m_imageSize[0]/2);
        currentIndex[0] = indexLeft;
        currentIndex[1] = m_imageIndex[1] + IndexValueType(m_imageSize[1]/2 + offset/m_imageSpacing[1]);
        while(indexLeft>m_imageIndex[0] && m_image->GetPixel(currentIndex) < edgeValue )
        {
            indexLeft--;
            currentIndex[0] = indexLeft;
        }
        // right
        IndexValueType indexRight = m_imageIndex[0] + IndexValueType(m_imageSize[0]/2);
        currentIndex[0] = indexRight;
        currentIndex[1] = m_imageIndex[1] + IndexValueType(m_imageSize[1]/2 + offset/m_imageSpacing[1]);
        while(indexRight<m_imageIndex[0]+m_imageSize[0] && m_image->GetPixel(currentIndex) < edgeValue )
        {
            indexRight++;
            currentIndex[0] = indexRight;
        }
        // bottom
        IndexValueType indexBot = m_imageIndex[1] + IndexValueType(m_imageSize[1]/2);
        currentIndex[0] = m_imageIndex[0] + IndexValueType(m_imageSize[0]/2 + offset/m_imageSpacing[1]);
        currentIndex[1] = indexBot;
        while(indexBot>m_imageIndex[1] && m_image->GetPixel(currentIndex) < edgeValue )
        {
            indexBot--;
            currentIndex[1] = indexBot;
        }
        // top
        IndexValueType indexTop = m_imageIndex[1] + IndexValueType(m_imageSize[1]/2);
        currentIndex[0] = m_imageIndex[0] + IndexValueType(m_imageSize[0]/2 + offset/m_imageSpacing[1]);
        currentIndex[1] = indexTop;
        while(indexTop<m_imageIndex[1]+m_imageSize[1] && m_image->GetPixel(currentIndex) < edgeValue )
        {
            indexTop++;
            currentIndex[1] = indexTop;
        }
        
        // Copy results to diaphragm region
        IndexType index;
        index[0] = indexLeft;
        index[1] = indexBot;
        SizeType size;
        size[0] = indexRight - indexLeft + 1;
        size[1] = indexTop - indexBot + 1;
        m_diaphragmRoi.SetIndex( index );
        m_diaphragmRoi.SetSize( size );
        
		// Compute diaphragm centre and size
		m_diaphragmCentre[0] = ( double(index[0] + size[0]/2) - double(m_imageIndex[0] + m_imageSize[0]/2) )* m_imageSpacing[0] / m_scaleFactor;
		m_diaphragmCentre[1] = ( double(index[1] + size[1]/2) - double(m_imageIndex[1] + m_imageSize[1]/2) )* m_imageSpacing[1] / m_scaleFactor;
		m_diaphragmSize[0] = size[0] * m_imageSpacing[0] / m_scaleFactor;
		m_diaphragmSize[1] = size[1] * m_imageSpacing[1] / m_scaleFactor;
        if( m_debug )
        {
            std::cout << "Diaphragm ROI " << m_diaphragmRoi << std::endl;
			std::cout << "Diaphragm centre ("
					  << m_diaphragmCentre[0] << ","
					  << m_diaphragmCentre[1] << ")" << std::endl;
			std::cout << "Diaphragm size ("
					  << m_diaphragmSize[0] << ","
					  << m_diaphragmSize[1] << ")" << std::endl;
        }
    }
    	
    void BuckyPhantom::FindOrientation()
    {
        if(m_debug) std::cout << "FindOrientation" << std::endl;
		// To find the orientation of the phantom in the image, a ROI
		// is placed at the location of the high contrast pattern of
		// the NORMI13 phantom at all possible locations (north, east,
		// south, west).
		//
		// +--------------------------------+
		// |                                |
		// |         HIGH  CONTRAST         |
		// |         P A T T E R N          |
		// |                                |
		// |                                |
		// |  /\                      dete  |
		// |  \/                      ctor  |
		// |                                |
		// |                                |
		// |         LOW   CONTRAST         |
		// |         P A T T E R N          |
		// |                                |
		// +--------------------------------+
		//
		// The high contrast pattern consists of seven 20x20 mm contrast steps.
		// When at the top of the image, the pattern's top left corner is
		// located (70, 80) mm away from the phantom's centre.
		// To find where the pattern is located in the image a ROI the size
		// of the pattern is placed at north, east, south and west locations
		// and a check is performed.

        IndexType index;
        SizeType  size;

		// Location and size of ROI
		double l[2];
		double s[2];

		// Check north
		// NORMI13 TERGOOI
        s[0] = 140.0 * m_scaleFactor;
        s[1] =  20.0 * m_scaleFactor;
        l[0] = -70.0 * m_scaleFactor;
        l[1] = -80.0 * m_scaleFactor;
		// DIGI13 MCL
//		s[0] = 126.0 * m_scaleFactor;
//		s[1] =  18.0 * m_scaleFactor;
//		l[0] = -63.0 * m_scaleFactor;
//		l[1] = (-63.0-18.0) * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        size[0]  = SizeValueType(s[0]/m_imageSpacing[0]);
        size[1]  = SizeValueType(s[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		// this->ShowRoi();
		if( this->CheckForHighContrastPattern(NORTH) )
		{
			m_orientation = NORTH;
			this->ShowRoi();
		}
			
		// Check east
		// NORMI13 TERGOOI
        s[0] =  20.0 * m_scaleFactor;
        s[1] = 140.0 * m_scaleFactor;
        l[0] = +60.0 * m_scaleFactor;
        l[1] = -70.0 * m_scaleFactor;
		// DIGI13 MCL
//		s[0] =  18.0 * m_scaleFactor;
//		s[1] = 126.0 * m_scaleFactor;
//		l[0] = +63.0 * m_scaleFactor;
//		l[1] = -63.0 * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        size[0]  = SizeValueType(s[0]/m_imageSpacing[0]);
        size[1]  = SizeValueType(s[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		// this->ShowRoi();
		if( this->CheckForHighContrastPattern(EAST) )
		{
			m_orientation = EAST;
			this->ShowRoi();
		}
		
		// Check south
		// NORMI13 TERGOOI
        s[0] = 140.0 * m_scaleFactor;
        s[1] =  20.0 * m_scaleFactor;
        l[0] = -70.0 * m_scaleFactor;
        l[1] = +60.0 * m_scaleFactor;
		// DIGI13 MCL
//		s[0] = 126.0 * m_scaleFactor;
//		s[1] =  18.0 * m_scaleFactor;
//		l[0] = -63.0 * m_scaleFactor;
//		l[1] = +63.0 * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        size[0]  = SizeValueType(s[0]/m_imageSpacing[0]);
        size[1]  = SizeValueType(s[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		// this->ShowRoi();
		if( this->CheckForHighContrastPattern(SOUTH) )
		{
			m_orientation = SOUTH;
			this->ShowRoi();
		}
		
		// Check west
		// NORMI13 TERGOOI
        s[0] =  20.0 * m_scaleFactor;
        s[1] = 140.0 * m_scaleFactor;
        l[0] = -80.0 * m_scaleFactor;
        l[1] = -70.0 * m_scaleFactor;
		// DIGI13 MCL
//		s[0] =  18.0 * m_scaleFactor;
//		s[1] = 126.0 * m_scaleFactor;
//		l[0] = (-63.0-18.0) * m_scaleFactor;
//		l[1] = -63.0 * m_scaleFactor;
        index[0] = IndexValueType(m_imageSize[0]/2 + l[0]/m_imageSpacing[0]);
        index[1] = IndexValueType(m_imageSize[1]/2 + l[1]/m_imageSpacing[1]);
        size[0]  = SizeValueType(s[0]/m_imageSpacing[0]);
        size[1]  = SizeValueType(s[1]/m_imageSpacing[1]);
        this->SetRoi( index, size );
        this->ComputeRoiStatistics();
		// this->ShowRoi();
		if( this->CheckForHighContrastPattern(WEST) )
		{
			m_orientation = WEST;
			this->ShowRoi();
		}
		
		if( m_debug )
		{
			std::cout << "Orientation of phantom is ";
			switch( m_orientation )
			{
			case NORTH: std::cout << "north";
				break;
			case EAST: std::cout << "east";
				break;
			case SOUTH: std::cout << "south";
				break;
			case WEST: std::cout << "west";
				break;
			default: std::cout << "UNKNOWN";
			}
			std::cout << std::endl;
		}
	}
    
	void BuckyPhantom::SetRoi( const IndexType index, const SizeType size )
    {
        if( m_debug ) std::cout << "SetRoi" << std::endl;
        m_roi.SetIndex( index );
        m_roi.SetSize( size );

        if( m_debug )
        {
            std::cout << "Index: " << m_roi.GetIndex() << std::endl;
            std::cout << "Size:  " << m_roi.GetSize() << std::endl;
        }
    }

    void BuckyPhantom::SetScaleFactor( const double thickness, const double SID )
    {
		if( m_debug ) std::cout << "SetScaleFactor" << std::endl;
		this->SetThickness( thickness );
		m_scaleFactor = SID / (SID - m_thickness);		
		
        if( m_debug )
        {
            std::cout << "scale factor: " << m_scaleFactor << std::endl;
        }
	}

    void BuckyPhantom::SetThickness( const double thickness )
    {
		if( m_debug ) std::cout << "SetThickness" << std::endl;
		m_thickness = thickness;
        if( m_debug )
        {
            std::cout << "thickness: " << m_thickness << std::endl;
        }
	}

    void BuckyPhantom::ShowDiaphragmEdges()
    {
		if( m_debug ) std::cout << "ShowDiaphragmEdges" << std::endl;
		m_roi = m_diaphragmRoi;
		this->ShowRoi();
		
		// // Create result image with inverted intensities outside diaphragm edges
		// RegionExclusionConstIteratorType iter1( m_inverseImage, m_inverseImage->GetLargestPossibleRegion() );
		// RegionExclusionIteratorType iter2( m_resultImage, m_resultImage->GetLargestPossibleRegion() );
		// iter1.SetExclusionRegion( m_diaphragmRoi );
		// iter2.SetExclusionRegion( m_diaphragmRoi );
		// while(!iter1.IsAtEnd())
		// {
			// iter2.Set( iter1.Get() );
			// ++iter1;
			// ++iter2;
		// }
	}

    void BuckyPhantom::ShowRoi()
    {
		if( m_debug ) std::cout << "ShowRoi" << std::endl;
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
