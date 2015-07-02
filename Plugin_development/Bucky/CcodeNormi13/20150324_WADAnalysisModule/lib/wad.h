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

#ifndef WAD_H
#define WAD_H

// Standard includes
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

// ITK includes
#include "itkImage.h"

namespace wad
{
    typedef signed short PixelType;
    const unsigned int   Dimension = 2;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef ImageType::IndexType IndexType;
    typedef IndexType::IndexValueType IndexValueType;
    typedef ImageType::RegionType RegionType;
    typedef ImageType::SizeType SizeType;
    typedef SizeType::SizeValueType SizeValueType;
    typedef ImageType::SpacingType SpacingType;
	
	struct DicomTagStruct
	{
		std::string tag;
		std::string type;
		unsigned level;
	};
	
	struct LimitsStruct
	{
		std::string quantity;
		std::string unit;
		std::string description;
		double      acceptableLow;
		double      acceptableHigh;
		double      criticalLow;
		double      criticalHigh;
	};
	
	struct ResultsPropertiesStruct
	{
		std::string type;
		unsigned    level;
		std::string quantity;
		std::string unit;
		std::string description;
	};
	
	template <typename T>
	std::string to_string(T const& value)
	{
		std::stringstream sstr;
		sstr << value;
		return sstr.str();
	}

} // namespace

#endif
