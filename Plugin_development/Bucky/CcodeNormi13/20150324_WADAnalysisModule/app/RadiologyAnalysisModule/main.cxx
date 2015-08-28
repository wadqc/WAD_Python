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

#include "main.h"

int main(int argc, char *argv[])
{
    // Check for correct number of command line arguments
	if( argc!=2 )
	{
		std::cerr << "Usage: " << argv[0] << " InputFile.xml" << std::endl;
		return EXIT_FAILURE;
	}

	bool debug = true;
	bool feedback = true;

	// Read the input file
	wad::InputFile inputFile( feedback );
	inputFile.SetFileName( argv[1] );
	inputFile.Read();
	if(debug) inputFile.Show();

	// Read the config file
	wad::RadiologyConfigFile configFile;
	configFile.SetDebug( debug );
	configFile.SetFeedback( feedback );
	configFile.SetFileName( inputFile.GetConfigFileName() );
	configFile.Read();

	// Read the DICOM file into an ITK image
	std::string dicomFileName = inputFile.GetInstanceFileName();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( dicomFileName );	
	GDCMImageIOType::Pointer gdcmImageIO = GDCMImageIOType::New();
	gdcmImageIO->SetMaxSizeLoadEntry(0xffff); // allow loading of long binary stream tags
	reader->SetImageIO( gdcmImageIO );
	std::cout << "Reading DICOM file " << dicomFileName << std::endl;
	try
	{
		reader->Update();
    }
	catch (itk::ExceptionObject & e)
    {
		std::cerr << e << std::endl;
		return EXIT_FAILURE;
    }	

	// Initialize a phantom object with the ITK image
    wad::RadiologyPhantom phantom;
    phantom.SetDebug( debug );
    phantom.SetFeedback( feedback );
    phantom.Initialize( reader->GetOutput() );

	// Analyze phantom image
	phantom.Analyze();

	// Write the result image
	phantom.WriteResultImage();
	
	// Write the results to an XML file
	wad::ResultFile resultFile( feedback );
	resultFile.SetFileName( inputFile.GetResultFileName() );
	resultFile.AppendPhantomResults( configFile.GetLimits(), configFile.GetResultsProperties(), phantom );
	resultFile.AppendDicomTagsResults( configFile.GetLimits(), configFile.GetDicomTags(), gdcmImageIO );
	resultFile.Write();
	
	return EXIT_SUCCESS;
}