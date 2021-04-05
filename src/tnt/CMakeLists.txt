#####################################################################
## CMakeLists.txt
## This file is part of the Comfor project
##
## Author: Eduardo Guzman Maldonado
## Copyright (C) 2019 - Innovamics
######################################################################
cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

target_sources(${EXE_TARGET} PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/tnt.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_version.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_math_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array1d.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array2d.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array1d_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array2d_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array1d_linalg.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array2d_linalg.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_tools.h)