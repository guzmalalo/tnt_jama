cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
#-----------------------------------------------------------------------------

## Set element test variable 
set(A1D_TARGET TNT_ARRAY1D)

## Set the target to this test 
add_executable(${A1D_TARGET} "")

# Add the sources to this test
target_include_directories(${A1D_TARGET} PRIVATE ${SOURCE_DIR}/tnt)

target_sources(${A1D_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/A1D_Struct.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/A1D_Basic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/A1D_Utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/A1D_Linalg.cpp
)

## Link libraries 
target_link_libraries(${A1D_TARGET} PUBLIC gtest_main)

## Add the test to Ctest
add_test(NAME ${A1D_TARGET}
         COMMAND ${A1D_TARGET}
)
