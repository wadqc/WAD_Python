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

#include "wadResultFile.h"

namespace wad
{
	void ResultFile::AppendChildNode( const std::string type, const unsigned level, const std::string value, const std::string quantity, const std::string unit, const std::string description, double acceptableLow, const double acceptableHigh, const double criticalLow, const double criticalHigh )
	{
		pugi::xml_node node;
		node = m_doc.child("WAD").append_child("results");
		node.append_child("volgnummer").text().set( m_nr );
		node.append_child("type").text().set( type.c_str() );
		node.append_child("niveau").text().set( level );
		node.append_child("waarde").text().set( value.c_str() );
		if( type == "float" )
		{
			node.append_child("grootheid").text().set( quantity.c_str() );
			node.append_child("eenheid").text().set( unit.c_str() );
		}
		node.append_child("omschrijving").text().set( description.c_str() );
		if( type == "float" )
		{
			node.append_child("grens_acceptabel_onder").text().set( acceptableLow );
			node.append_child("grens_acceptabel_boven").text().set( acceptableHigh );
			node.append_child("grens_kritisch_onder").text().set( criticalLow );
			node.append_child("grens_kritisch_boven").text().set( criticalHigh );
		}		
		++m_nr;
	}

	void ResultFile::AppendDicomTagsResults( const std::vector<LimitsStruct> &limits, const std::vector<DicomTagStruct> &dicomTags, ImageIOType *dicomIO )
	{
		 // allow loading of long binary stream tags
		 dicomIO->SetMaxSizeLoadEntry(0xffff);
		
		for(unsigned i=0; i<dicomTags.size();++i)
		{
			std::string type        = dicomTags.at(i).type;
			unsigned    level       = dicomTags.at(i).level;
			std::string value;
            std::string quantity    = "";
            std::string unit        = "";
			std::string description;
            double acceptableLow  = 0.0;
            double acceptableHigh = 0.0;
            double criticalLow    = 0.0;
            double criticalHigh   = 0.0;

			std::string tag = dicomTags.at(i).tag;
			if( itk::GDCMImageIO::GetLabelFromTag( tag, description ) )
			{
				if( dicomIO->GetValueFromTag( tag, value ) )
				{
                    // Find the corresponding limits using only the description string
                    // NOTE: quantity and unit string are set from the corresponding limits fields
                    if( type == "float" )
                    {
                        for(unsigned j=0; j<limits.size();++j)
                        {
                            if( limits.at(j).description == description )
                            {
                                quantity       = limits.at(j).quantity;
                                unit           = limits.at(j).unit;
                                acceptableLow  = limits.at(j).acceptableLow;
                                acceptableHigh = limits.at(j).acceptableHigh;
                                criticalLow    = limits.at(j).criticalLow;
                                criticalHigh   = limits.at(j).criticalHigh;
                            }
                        }
                    }

                    AppendChildNode( type, level, value, quantity, unit, description, acceptableLow, acceptableHigh, criticalLow, criticalHigh );
				}
				else
				{
					std::cout << "(No Value Found in File)" << std::endl;
				}
			}
			else
			{
				std::cerr << "Trying to access non existing DICOM tag: " << tag << std::endl;
			}
		}
	}

	void ResultFile::AppendPhantomResults( const std::vector<LimitsStruct> &limits, const std::vector<ResultsPropertiesStruct> &resultsProperties, Phantom &phantom )
	{
		for(unsigned i=0; i<resultsProperties.size();++i)
		{
			std::string type        = resultsProperties.at(i).type;
			unsigned    level       = resultsProperties.at(i).level;
            std::string quantity    = resultsProperties.at(i).quantity; // only non-empty for "float" type
            std::string unit        = resultsProperties.at(i).unit;     // only non-empty for "float" type
			std::string description = resultsProperties.at(i).description;
            double acceptableLow  = 0.0;
            double acceptableHigh = 0.0;
            double criticalLow    = 0.0;
            double criticalHigh   = 0.0;
            std::string value = phantom.GetResult( quantity, unit, description );

            // Find the corresponding limits using the quantity, unit and description strings
			if( type == "float" )
			{
                for(unsigned j=0; j<limits.size();++j)
                {
                    if( limits.at(j).quantity    == quantity    &&
                        limits.at(j).unit        == unit        &&
                        limits.at(j).description == description )
                    {
                        acceptableLow  = limits.at(j).acceptableLow;
                        acceptableHigh = limits.at(j).acceptableHigh;
                        criticalLow    = limits.at(j).criticalLow;
                        criticalHigh   = limits.at(j).criticalHigh;
                    }
                }
			}
            
            AppendChildNode( type, level, value, quantity, unit, description, acceptableLow, acceptableHigh, criticalLow, criticalHigh );
		}
	}
	
	void ResultFile::CreateWadNode()
	{
		m_doc.append_child("WAD");
		
		// Counter for keeping track of child nodes
		// Should start at "1" according to WAD specs...
		m_nr = 1;
	}

	void ResultFile::SetFileName( const std::string fileName )
	{
		m_fileName = fileName;
	}

	void ResultFile::Write()
	{
		if( m_doc.save_file( m_fileName.c_str() ) )
		{
			std::cout << "Saved results to " << m_fileName << std::endl;
		}
		else
		{
			std::cerr << "Failed to write to " << m_fileName << std::endl;
		}
	}

} // end namespace
