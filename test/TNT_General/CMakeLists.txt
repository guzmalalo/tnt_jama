cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
#-----------------------------------------------------------------------------

## Set element test variable 
set(TNT_GEN_TARGET TNT_GENERAL)

## Set the target to this test 
add_executable(${TNT_GEN_TARGET} "")

# Add the sources to this test
target_include_directories(${TNT_GEN_TARGET} PRIVATE ${SOURCE_DIR}/tnt)

target_sources(${TNT_GEN_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/TNT_General.cpp
)

## Link libraries 
target_link_libraries(${TNT_GEN_TARGET} PUBLIC gtest_main)

## Add the test to Ctest
add_test(NAME ${TNT_GEN_TARGET}
         COMMAND ${TNT_GEN_TARGET}
)
