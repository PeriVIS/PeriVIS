CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CountBonds)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( CountBonds "1.0"
  SERVER_MANAGER_XML CountBonds.xml
  SERVER_MANAGER_SOURCES CountBonds.cxx  
)
  
INSTALL( TARGETS CountBonds CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS CountBonds CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo 
                                RUNTIME DESTINATION bin/RelWithDebInfo 
                                LIBRARY DESTINATION bin/RelWithDebInfo 
                                )