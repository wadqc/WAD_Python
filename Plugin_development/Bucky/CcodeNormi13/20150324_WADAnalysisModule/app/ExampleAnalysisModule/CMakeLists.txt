SET( WADExampleAnalysisModule_INCLUDES
  main.h
  wadExampleConfigFile.h
  wadExamplePhantom.h
)

SET( WADExampleAnalysisModule_SOURCES
  main.cxx
  wadExampleConfigFile.cxx
  wadExamplePhantom.cxx
)

ADD_EXECUTABLE( WADExampleAnalysisModule ${WADExampleAnalysisModule_INCLUDES} ${WADExampleAnalysisModule_SOURCES} )

TARGET_LINK_LIBRARIES( WADExampleAnalysisModule
  WADAnalysisModule
  ${ITK_LIBRARIES}
)

INSTALL_TARGETS(/bin WADExampleAnalysisModule)