cmake_minimum_required(VERSION 2.8)

project(PV-Plugins)

FIND_PACKAGE(ParaView REQUIRED)


option(BUILD_IO_LAMMPS "Build reader and writer for LAMMPS file format" ON)
if (BUILD_IO_LAMMPS)
	add_subdirectory( io )
endif()

option(BUILD_FRAGMENT "Build the plugins for extraction of fragments" OFF)
if (BUILD_FRAGMENT)
	add_subdirectory( fragment )
endif()

option(BUILD_STRESS "Build the plugins for visualization of the stress" OFF)
if (BUILD_STRESS)
	add_subdirectory( stress )
endif()
