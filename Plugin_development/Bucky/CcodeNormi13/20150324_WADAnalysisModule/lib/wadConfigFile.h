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

#ifndef WADCONFIGFILE_H
#define WADCONFIGFILE_H

// Custom includes
#include "wad.h"
#include "pugixml.hpp"

namespace wad
{
/**
 * \class ConfigFile
 * \brief Abstract base class for configuration file
 *
 * The pugixml library is used to read the XML config file.
 */
	class ConfigFile
	{
	public:
		std::vector<DicomTagStruct> GetDicomTags()
		{
			return m_dicomTags;
		}
		
		std::vector<LimitsStruct> GetLimits()
		{
			return m_limits;
		}
		
		std::vector<ResultsPropertiesStruct> GetResultsProperties()
		{
			return m_resultsProperties;
		}

        /// Read the config xml file
        void Read()
    	{
	        if( m_debug ) std::cout << "ConfigFile::Read" << std::endl;
    		if ( !m_doc.load_file( m_fileName.c_str() ) )
    		{
    			std::cerr << "Failed to read " << m_fileName << std::endl;
    		}
    		else
    		{
    			this->SetDicomTags();
                this->SetLimits();
    			this->SetResultsProperties();
                this->SetParameters();
    		}
    	}

        void SetDebug( const bool debug )
        {
            m_debug = debug;
        }

        void SetFeedback( const bool feedback )
        {
            m_feedback = feedback;
        }

        /// Set the file name
		void SetFileName( const std::string fileName )
		{
			m_fileName = fileName;
		}

    protected:
		/// Default constructor
		ConfigFile()
        {
            m_debug = false;
            m_feedback = true;
        }

        virtual ~ConfigFile() {}
        
        bool m_debug;
        bool m_feedback;

		/// DICOM tags
		std::vector<DicomTagStruct> m_dicomTags;
		
		/// Limits
		std::vector<LimitsStruct> m_limits;

		/// Results properties
		std::vector<ResultsPropertiesStruct> m_resultsProperties;

        /// File name of the config file
        std::string m_fileName;

		/// XML document
	    pugi::xml_document m_doc;

    private:
		/// Set the DICOM tags
		void SetDicomTags()
		{
	        if( m_debug ) std::cout << "ConfigFile::SetDicomTags" << std::endl;

			DicomTagStruct tagStruct;
			pugi::xml_node tags = m_doc.child("WAD_config").child("dicomtags");
			for (pugi::xml_node tag = tags.child("tag"); tag; tag = tag.next_sibling("tag"))
			{
				tagStruct.type  = tag.attribute("type").value();
				tagStruct.level = tag.attribute("level").as_int();
				tagStruct.tag = tag.attribute("tag").value();
				m_dicomTags.push_back( tagStruct );

				if( m_feedback )
				{
					std::cout << "DICOM tag: "
						<< "type " << tagStruct.type
						<< ", level " << tagStruct.level
						<< ", tag " << tagStruct.tag << std::endl;
				}
			}
		}

		/// Set the Limits
		void SetLimits()
		{
	        if( m_debug ) std::cout << "ConfigFile::SetLimits" << std::endl;
            
			LimitsStruct limitsStruct;
			pugi::xml_node wad = m_doc.child("WAD_config");
			for (pugi::xml_node limit = wad.child("limit"); limit; limit = limit.next_sibling("limit"))
			{
				limitsStruct.quantity = limit.child("quantity").text().get();
				limitsStruct.unit = limit.child("unit").text().get();
				limitsStruct.description = limit.child("description").text().get();
                limitsStruct.acceptableLow = atof( limit.child("acceptable_low").text().get() );
                limitsStruct.acceptableHigh = atof( limit.child("acceptable_high").text().get() );
                limitsStruct.criticalLow = atof( limit.child("critical_low").text().get() );
                limitsStruct.criticalHigh = atof( limit.child("critical_high").text().get() );
				m_limits.push_back( limitsStruct );
                
				if( m_feedback )
				{
					std::cout << "Limits: " << std::endl
                    << "quantity: " << limitsStruct.quantity << std::endl
                    << "unit: " << limitsStruct.unit << std::endl
                    << "description: " << limitsStruct.description << std::endl
                    << "acceptableLow: " << limitsStruct.acceptableLow << std::endl
                    << "acceptableHigh: " << limitsStruct.acceptableHigh << std::endl
                    << "criticalLow: " << limitsStruct.criticalLow << std::endl
                    << "criticalHigh: " << limitsStruct.criticalHigh << std::endl;
				}
			}
		}

		/// Set the parameters
		virtual void SetParameters() = 0;

		/// Set the results properties
		virtual void SetResultsProperties() = 0;

	}; // end of class
    
} // end of namespace

#endif
