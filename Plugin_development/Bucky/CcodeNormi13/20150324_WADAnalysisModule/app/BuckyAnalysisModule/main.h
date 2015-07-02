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

// WAD includes
#include "wad.h"
#include "wadBuckyConfigFile.h"
#include "wadBuckyPhantom.h"
#include "wadInputFile.h"
#include "wadResultFile.h"

// ITK includes
#include "itkImageFileReader.h"
#include "itkGDCMImageIO.h"

typedef itk::ImageFileReader< wad::ImageType > ReaderType;
typedef itk::GDCMImageIO ImageIOType;
