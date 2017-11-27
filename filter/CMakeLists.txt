cmake_minimum_required(VERSION 2.8)

project(PV-Plugins)

FIND_PACKAGE(ParaView REQUIRED)

option(BUILD_FRAGMENT "Build the plugins for extraction of fragments" OFF)
if (BUILD_IO_LAMMPS)
	add_subdirectory( fragment )
endif()

