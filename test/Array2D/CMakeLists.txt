cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
#-----------------------------------------------------------------------------

## Set element test variable 
set(A2D_TARGET TNT_ARRAY2D)

## Set the target to this test 
add_executable(${A2D_TARGET} "")

# Add the sources to this test
target_include_directories(${A2D_TARGET} PRIVATE ${SOURCE_DIR}/tnt)

target_sources(${A2D_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/A2D_Struct.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/A2D_Basic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/A2D_Utils.cpp
)

## Link libraries 
target_link_libraries(${A2D_TARGET} PUBLIC gtest_main)

## Add the test to Ctest
add_test(NAME ${A2D_TARGET}
         COMMAND ${A2D_TARGET}
)
