cmake_minimum_required(VERSION 3.12)

#-------------------------------------------------------------------------------
# Project name and version
#-------------------------------------------------------------------------------
set(NUM_VERSION_MAJOR "3")
set(NUM_VERSION_MINOR "0")
set(NUM_VERSION_PATCH "2")
set(NUM_VERSION  "${NUM_VERSION_MAJOR}.${NUM_VERSION_MINOR}.${NUM_VERSION_PATCH}" 
    CACHE INTERNAL "Software version number")

project(TNT VERSION ${NUM_VERSION} 
        DESCRIPTION "Template Numerical Kit"
        LANGUAGES CXX)


#-----------------------------------------------------------------------------
# C++17 Support
# Require C++17 and disable compiler-specific extensions
#-----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#-----------------------------------------------------------------------------
#                               INNER BUILD
#-----------------------------------------------------------------------------
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOC_DIR ${PROJECT_SOURCE_DIR}/doc)

#-----------------------------------------------------------------------------
# Add the cmake folder so the FindSphinx module is found
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#-----------------------------------------------------------------------------
# Define the target : executable
#-----------------------------------------------------------------------------
set(EXE_TARGET ${PROJECT_NAME})
if (MSVC)
    # MSVC creates a folder for each build type
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_LIST_DIR}/Bin")
else()
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_LIST_DIR}/Bin/${CMAKE_BUILD_TYPE}")
endif()

#Creates the executable target
add_executable(${EXE_TARGET} "")

#-----------------------------------------------------------------------------
# Add Source code from subdirectories
#-----------------------------------------------------------------------------
add_subdirectory(${SOURCE_DIR})

#-----------------------------------------------------------------------------
# Add source for doc generation
#-----------------------------------------------------------------------------
if(1)
add_subdirectory(${DOC_DIR})
endif()

