SET( WADBuckyQCAnalysisModule_INCLUDES
  main.h
  wadBuckyConfigFile.h
  wadBuckyPhantom.h
)

SET( WADBuckyAnalysisModule_SOURCES
  main.cxx
  wadBuckyConfigFile.cxx
  wadBuckyPhantom.cxx
)


ADD_EXECUTABLE( WADBuckyAnalysisModule ${WADBuckyAnalysisModule_INCLUDES} ${WADBuckyAnalysisModule_SOURCES} )

TARGET_LINK_LIBRARIES( WADBuckyAnalysisModule
  WADAnalysisModule
  ${ITK_LIBRARIES}
)

INSTALL_TARGETS(/bin WADBuckyAnalysisModule)
