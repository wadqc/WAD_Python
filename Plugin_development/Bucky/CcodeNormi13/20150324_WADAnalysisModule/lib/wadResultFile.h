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

#ifndef WADRESULTFILE_H
#define WADRESULTFILE_H

// Custom includes
#include "wad.h"
#include "pugixml.hpp"
#include "wadPhantom.h"

// ITK includes
#include "itkGDCMImageIO.h"

namespace wad
{
/**
 * \class ResultFile
 * \brief Results of the QC analysis
 *
 * The pugixml library is used to write the XML file.
 */
	typedef itk::GDCMImageIO ImageIOType;

	class ResultFile
	{
	public:
		/// Default constructor
		ResultFile()
        {
            m_feedback = false;
			this->CreateWadNode();
        }

		/// Constructor for setting feedback
        ResultFile( const bool feedback )
        {
            m_feedback = feedback;
			this->CreateWadNode();
        }

		/// Append the results read from the DICOM IO object using the DICOM tags
		void AppendDicomTagsResults( const std::vector<LimitsStruct> &limits, const std::vector<DicomTagStruct> &dicomTags, ImageIOType *dicomIO );
		
		/// Append the results from the phantom analysis
		void AppendPhantomResults( const std::vector<LimitsStruct> &limits, const std::vector<ResultsPropertiesStruct> &resultsProperties, Phantom &phantom );
		
        /// Set the file name of the result file
		void SetFileName( const std::string fileName );
        
        /// Write the results to the xml file
        void Write();
        
    protected:
        
    private:
		/// Append a child node to the main WAD node
		void AppendChildNode( const std::string type, const unsigned level, const std::string value, const std::string quantity, const std::string unit, const std::string description, double acceptLow, const double acceptHigh, const double criticalLow, const double criticalHigh );

		/// Creates the main WAD node
		void CreateWadNode();

		/// XML document
	    pugi::xml_document m_doc;

        /// Feedback switch
        bool m_feedback;
        
        /// File name of the result file
        std::string m_fileName;
		
		/// Counter for keeping track of child nodes
		unsigned m_nr;

	}; // end of class
    
} // end of namespace

#endif
