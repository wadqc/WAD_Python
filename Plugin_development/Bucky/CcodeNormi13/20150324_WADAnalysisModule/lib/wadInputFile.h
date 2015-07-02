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

#ifndef WADINPUTFILEREADER_H
#define WADINPUTFILEREADER_H

// Custom includes
#include "wad.h"
#include "pugixml.hpp"

namespace wad
{
/**
 * \class InputFile
 * \brief Input parameters for QC analysis
 *
 * The pugixml library is used to read the XML file.
 * Currently, the input parameters for the analysis of a NORMI13 phantom
 * are stored.
 */
	class InputFile
	{
	public:
		/// Default constructor
		InputFile()
        {
            m_feedback = false;
        }

		/// Constructor for setting feedback
        InputFile( const bool feedback )
        {
            m_feedback = feedback;
        }

		std::string GetAnalysisLevel() const;
		std::string GetConfigFileName() const;
		std::string GetInstanceFileName() const;
		std::string GetInstanceNumber() const;
		std::string GetPatientId() const;
		std::string GetPatientName() const;
		std::string GetResultFileName() const;
		std::string GetSeriesDescription() const;
		std::string GetSeriesNumber() const;
		std::string GetStudyDescription() const;
		std::string GetStudyUid() const;
		std::string GetVersion() const;

        /// Set the file name
		void SetFileName( const std::string fileName );
        
        /// Read the input XML file
        void Read();
        
        /// Show the contents of the input file
		void Show();
		
    protected:
        
    private:
        /// Feedback switch
        bool m_feedback;
        
        /// File name of the result file
        std::string m_fileName;

		/// XML document
	    pugi::xml_document m_doc;

		/// Main wad node
		pugi::xml_node m_wadNode;

	}; // end of class
    
} // end of namespace

#endif
