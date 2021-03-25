#####################################################################
## CMakeLists.txt
## This file is part of the Comfor project
##
## Author: Eduardo Guzman Maldonado
## Copyright (C) 2019 - Innovamics
######################################################################
cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tnt)
target_include_directories(${EXE_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tnt)

target_sources(${EXE_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Main_tnt.cpp)

