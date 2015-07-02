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

#ifndef WADPHANTOM_H
#define WADPHANTOM_H

// ITK includes
#include <itkImageDuplicator.h>

// Custom includes
#include "wad.h"

namespace wad
{
/**
 * \class Phantom
 * \brief Abstract base class for representation of a QC test object (i.e. a "phantom")
 *
 * An itk::Image is used to specify the phantom.
 *
 * \todo Template over Dimension
 */
	class Phantom
	{
	public:
		typedef itk::ImageDuplicator< ImageType > DuplicatorType;

        void Analyze()
        {
            if( m_debug ) std::cout << "Phantom::Analyze" << std::endl;
            if( m_initialized )
            {
                if( m_feedback ) std::cout << "Analyzing phantom" << std::endl;
                DoAnalyze();
                m_analyzed = true;
            }
            else
            {
                std::cerr << "ERROR: Can't analyze Phantom, first call Initialize" << std::endl;
            }
        }

		std::string GetResult( const std::string quantity, const std::string unit, const std::string description )
		{
			std::string result = "";
            if( m_debug ) std::cout << "Phantom::GetResult" << std::endl;
            if( m_analyzed )
            {
                if( m_feedback ) std::cout << "Getting result" << std::endl;
                result = DoGetResult( quantity, unit, description );
            }
            else
            {
                std::cerr << "ERROR: Can't get result, first call call Analyze" << std::endl;
            }
            return result;
		}

		void Initialize( const ImageType *image )
        {
            if( m_debug ) std::cout << "Phantom::Initialize" << std::endl;
            if( m_feedback ) std::cout << "Initializing phantom" << std::endl;
            m_image = image;
            m_imageIndex   = m_image->GetLargestPossibleRegion().GetIndex();
            m_imageSize    = m_image->GetLargestPossibleRegion().GetSize();
            m_imageSpacing = m_image->GetSpacing();

            // Duplicate input image to create result image
            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(image);
    		duplicator->Update();
    		m_resultImage = duplicator->GetOutput();

            if( m_debug )
    		{
    			std::cout << "Image index   " << m_imageIndex   << std::endl;
    			std::cout << "Image size    " << m_imageSize    << std::endl;
    			std::cout << "Image spacing " << m_imageSpacing << std::endl;
    		}

        	DoInitialize( image );
            m_initialized = true;
        }

        void SetDebug( const bool debug )
        {
            m_debug = debug;
        }

        void SetFeedback( const bool feedback )
        {
            m_feedback = feedback;
        }

        void WriteResultImage()
        {
            if( m_debug ) std::cout << "Phantom::WriteResultImage" << std::endl;
            if( m_analyzed )
            {
                if( m_feedback ) std::cout << "Writing result image" << std::endl;
                DoWriteResultImage();
            }
            else
            {
                std::cerr << "ERROR: Can't write result image, first call Analyze" << std::endl;
            }
        }

	protected:
		/// Default constructor
		Phantom()
        {
            m_analyzed = false;
            m_debug = false;
            m_feedback = true;
            m_initialized = false;
        }

        virtual ~Phantom() {}
        
        bool m_analyzed;
        bool m_debug;
        bool m_feedback;
        bool m_initialized;

        /// ITK image to store input DICOM image of phantom
        ImageType::ConstPointer m_image;

        /// ITK image to store analyzed result image
        ImageType::Pointer m_resultImage;
		
        /// Image index (lower-left corner)
        IndexType m_imageIndex;
        
        /// Image size (in pixels)
        SizeType m_imageSize;
        
        /// Image spacing
        SpacingType m_imageSpacing;

    private:
		virtual void DoAnalyze() = 0;

		virtual std::string DoGetResult( const std::string quantity, const std::string unit, const std::string description ) = 0;

		virtual void DoInitialize( const ImageType *image ) = 0;

		virtual void DoWriteResultImage() = 0;

	}; // end of class
    
} // end of namespace

#endif
