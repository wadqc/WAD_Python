SET( WADRadiologyAnalysisModule_INCLUDES
  main.h
  wadRadiologyConfigFile.h
  wadRadiologyPhantom.h
)

SET( WADRadiologyAnalysisModule_SOURCES
  main.cxx
  wadRadiologyConfigFile.cxx
  wadRadiologyPhantom.cxx
)

ADD_EXECUTABLE( WADRadiologyAnalysisModule ${WADRadiologyAnalysisModule_INCLUDES} ${WADRadiologyAnalysisModule_SOURCES} )

TARGET_LINK_LIBRARIES( WADRadiologyAnalysisModule
  WADAnalysisModule
  ${ITK_LIBRARIES}
)

INSTALL_TARGETS(/bin WADRadiologyAnalysisModule)