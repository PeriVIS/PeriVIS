#MIT License
#
#Copyright (c) 2021 Michael Bussler
#Copyright (c) 2021 Patrick Diehl
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cmake_minimum_required(VERSION 2.8)
 
PROJECT(LAMMPSTemporalReader)
 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

ADD_PARAVIEW_PLUGIN(LAMMPSTemporalReader "1.0"
  SERVER_MANAGER_XML LAMMPSTemporalReader.xml
  SERVER_MANAGER_SOURCES LAMMPSTemporalReader.cxx
  REQUIRED_ON_SERVER)

TARGET_LINK_LIBRARIES( LAMMPSTemporalReader PRIVATE pqCore )
  
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )
