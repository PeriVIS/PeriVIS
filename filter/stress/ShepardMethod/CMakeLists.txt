cmake_minimum_required(VERSION 2.8)

PROJECT(ShepardMethod)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( ShepardMethod "1.0"
  SERVER_MANAGER_XML ShepardMethod.xml
  SERVER_MANAGER_SOURCES vtkShepardTensor.cxx  
  )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
IF(MSVC)
  INSTALL ( FILES  ${PROJECT_BINARY_DIR}/Debug/${PROJECT_NAME}.pdb  CONFIGURATIONS Debug
                                                                    DESTINATION bin/Debug
                  
  )
ENDIF(MSVC)