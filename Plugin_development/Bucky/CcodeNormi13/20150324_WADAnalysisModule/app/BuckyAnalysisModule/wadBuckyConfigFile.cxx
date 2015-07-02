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

#include "wadBuckyConfigFile.h"

namespace wad
{
	void BuckyConfigFile::SetParameters()
	{
        if( m_debug ) std::cout << "BuckyConfigFile::SetParameters" << std::endl;
		
		// Copy parameters to local variables
		pugi::xml_node parameters = m_doc.child("WAD_config").child("parameters");
		m_phantomThickness = atof( parameters.child("PhantomThickness").text().get() );
        if( m_feedback )
        {
			std::cout << m_phantomThickness << std::endl;
		}
	}

	void BuckyConfigFile::SetResultsProperties()
	{
        /// Defines the quantity, unit and description fields
        /// For the limits to be set correctly, the XML configuration file should use the exact same strings.
        if( m_debug ) std::cout << "BuckyConfigFile::SetResultsProperties" << std::endl;

		ResultsPropertiesStruct properties;

		DicomTagStruct tagStruct;
		pugi::xml_node limits = m_doc.child("WAD_config").child("limits");

		properties.type = "float";
		properties.level = 1;
		properties.quantity = "length";
		properties.unit = "mm";
		properties.description = "diaphragm center x";
		m_resultsProperties.push_back( properties );
		
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "length";
		properties.unit = "mm";
		properties.description = "diaphragm center y";
		m_resultsProperties.push_back( properties );

		properties.type = "float";
		properties.level = 1;
		properties.quantity = "length";
		properties.unit = "mm";
		properties.description = "diaphragm size x";
		m_resultsProperties.push_back( properties );

		properties.type = "float";
		properties.level = 1;
		properties.quantity = "length";
		properties.unit = "mm";
		properties.description = "diaphragm size y";
		m_resultsProperties.push_back( properties );

		properties.type = "char";
		properties.level = 1;
        properties.quantity = "";
		properties.unit = "";
        properties.description = "orientation";
		m_resultsProperties.push_back( properties );
        
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "";
		properties.unit = "%";
		properties.description = "inhomogeneity top left";
		m_resultsProperties.push_back( properties );
        
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "";
		properties.unit = "%";
		properties.description = "inhomogeneity top right";
		m_resultsProperties.push_back( properties );
        
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "";
		properties.unit = "%";
		properties.description = "inhomogeneity center";
		m_resultsProperties.push_back( properties );
        
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "";
		properties.unit = "%";
		properties.description = "inhomogeneity bottom left";
		m_resultsProperties.push_back( properties );
        
		properties.type = "float";
		properties.level = 1;
		properties.quantity = "";
		properties.unit = "%";
		properties.description = "inhomogeneity bottom right";
		m_resultsProperties.push_back( properties );
		
	}
	
} // end namespace
