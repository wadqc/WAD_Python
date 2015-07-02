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

#ifndef WADEXAMPLECONFIGFILE_H
#define WADEXAMPLECONFIGFILE_H

// Custom includes
#include "wadConfigFile.h"

namespace wad
{
/**
 * \class ExampleConfigFile
 * \brief Derived class implementing example configuration parameters
 *
 */
	class ExampleConfigFile : public ConfigFile
	{
	public:
        
    protected:
        
    private:
		/// Set the parameters
		void SetParameters();

		/// Set the results properties
		void SetResultsProperties();
		
	}; // end of class
    
} // end of namespace

#endif
