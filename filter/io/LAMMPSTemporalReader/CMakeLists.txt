cmake_minimum_required(VERSION 2.8)
 
PROJECT(LAMMPSTemporalReader)
 
FIND_PACKAGE(ParaView REQUIRED)

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

paraview_plugin_scan(
  PLUGIN_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  PROVIDES_PLUGINS  plugins
  REQUIRES_MODULES  required_modules)


#PARAVIEW_ADD_PLUGIN(LAMMPSTemporalReader "1.0"
#  SERVER_MANAGER_XML LAMMPSTemporalReader.xml
#  SERVER_MANAGER_SOURCES LAMMPSTemporalReader.cxx
#  REQUIRED_ON_SERVER)

paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  PLUGINS ${plugins})

#TARGET_LINK_LIBRARIES( LAMMPSTemporalReader PRIVATE pqCore )
  
#INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
#                                RUNTIME DESTINATION bin/Release 
#                                LIBRARY DESTINATION bin/Release
#                                )

#INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
#                                 RUNTIME DESTINATION bin/Debug 
#                                 LIBRARY DESTINATION bin/Debug
#    )
                             #INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )
