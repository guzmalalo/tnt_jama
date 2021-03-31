cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

#-----------------------------------------------------------------------------
# Googletest framework for unit-tests - C++
#-----------------------------------------------------------------------------
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
  googletest

  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        master)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(
    ${googletest_SOURCE_DIR}
    ${googletest_BINARY_DIR}
    EXCLUDE_FROM_ALL)
endif()

## Include the test cases 
add_subdirectory(TNT_General)
add_subdirectory(Array1D)
add_subdirectory(Array2D)