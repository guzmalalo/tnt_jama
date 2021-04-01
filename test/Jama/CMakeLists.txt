cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
#-----------------------------------------------------------------------------

## Set element test variable 
set(JAMA_TARGET JAMA_GENERAL)

## Set the target to this test 
add_executable(${JAMA_TARGET} "")

# Add the sources to this test
target_include_directories(${JAMA_TARGET} PRIVATE ${SOURCE_DIR}/tnt)
target_include_directories(${JAMA_TARGET} PRIVATE ${SOURCE_DIR}/jama)

target_sources(${JAMA_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Jama_Struct.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Jama_LU.cpp
)

## Link libraries 
target_link_libraries(${JAMA_TARGET} PUBLIC gtest_main)

## Add the test to Ctest
add_test(NAME ${JAMA_TARGET}
         COMMAND ${JAMA_TARGET}
)
