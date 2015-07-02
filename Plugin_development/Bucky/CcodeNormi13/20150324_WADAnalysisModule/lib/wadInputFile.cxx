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

#include "wadInputFile.h"

namespace wad
{
	std::string InputFile::GetAnalysisLevel() const
	{
		std::string result = m_wadNode.child("analyselevel").text().get();
		return result;
	}

	std::string InputFile::GetConfigFileName() const
	{
		std::string result = m_wadNode.child("analysemodule_cfg").text().get();
		return result;
	}

	std::string InputFile::GetInstanceFileName() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("series").child("instance").child("filename").text().get();
		return result;
	}

	std::string InputFile::GetInstanceNumber() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("series").child("instance").child("number").text().get();
		return result;
	}

	std::string InputFile::GetPatientId() const
	{
		std::string result = m_wadNode.child("patient").child("id").text().get();
		return result;
	}

	std::string InputFile::GetPatientName() const
	{
		std::string result = m_wadNode.child("patient").child("name").text().get();
		return result;
	}

	std::string InputFile::GetResultFileName() const
	{
		std::string result = m_wadNode.child("analysemodule_output").text().get();
		return result;
	}

	std::string InputFile::GetSeriesDescription() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("series").child("description").text().get();
		return result;
	}

	std::string InputFile::GetSeriesNumber() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("series").child("number").text().get();
		return result;
	}

	std::string InputFile::GetStudyDescription() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("description").text().get();
		return result;
	}

	std::string InputFile::GetStudyUid() const
	{
		std::string result = m_wadNode.child("patient").child("study").child("uid").text().get();
		return result;
	}

	std::string InputFile::GetVersion() const
	{
		std::string result = m_wadNode.child("version").text().get();
		return result;
	}

	void InputFile::Read()
	{
		// Read the XML file
		if ( !m_doc.load_file( m_fileName.c_str() ) )
		{
			std::cerr << "Failed to read " << m_fileName << std::endl;
		}
		else
		{
			// Set the main wad node
			m_wadNode = m_doc.child("WAD");
		}
	}

	void InputFile::SetFileName( const std::string fileName )
	{
		m_fileName = fileName;
	}

	void InputFile::Show()
	{
		std::cout << this->GetAnalysisLevel() << std::endl;
		std::cout << this->GetConfigFileName() << std::endl;
		std::cout << this->GetInstanceFileName() << std::endl;
		std::cout << this->GetInstanceNumber() << std::endl;
		std::cout << this->GetPatientId() << std::endl;
		std::cout << this->GetPatientName() << std::endl;
		std::cout << this->GetResultFileName() << std::endl;
		std::cout << this->GetSeriesDescription() << std::endl;
		std::cout << this->GetSeriesNumber() << std::endl;
		std::cout << this->GetStudyDescription() << std::endl;
		std::cout << this->GetStudyUid() << std::endl;
		std::cout << this->GetVersion() << std::endl;
	}
	
} // end namespace
