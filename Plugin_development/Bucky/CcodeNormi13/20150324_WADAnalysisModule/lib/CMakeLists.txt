SET( WADAnalysisModule_INCLUDES
  wad.h
  wadConfigFile.h
  wadInputFile.h
  wadPhantom.h
  wadResultFile.h
  pugixml.hpp
)

SET( WADAnalysisModule_SOURCES
  wadInputFile.cxx
  wadResultFile.cxx
  pugixml.cpp
)

ADD_LIBRARY(WADAnalysisModule ${WADAnalysisModule_INCLUDES} ${WADAnalysisModule_SOURCES} )

TARGET_LINK_LIBRARIES(WADAnalysisModule ${ITK_LIBRARIES})

INSTALL( TARGETS WADAnalysisModule
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
