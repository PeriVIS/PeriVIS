CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ExtractBrokenBonds)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

PARAVIEW_ADD_PLUGIN( ExtractBrokenBonds "1.0"
  SERVER_MANAGER_XML ExtractBrokenBonds.xml
  SERVER_MANAGER_SOURCES ExtractBrokenBonds.cxx  
)
  
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo 
                                 RUNTIME DESTINATION bin/RelWithDebInfo 
                                 LIBRARY DESTINATION bin/RelWithDebInfo 
                                 )