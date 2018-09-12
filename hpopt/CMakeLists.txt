
#-----------------------------------------------------------#

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

#-----------------------------------------------------------#

project(HPOPT)

#-----------------------------------------------------------#

set(CMAKE_CXX_COMPILER "g++")

#-----------------------------------------------------------#

set(CMAKE_CXX_FLAGS "-O3 -g -Wall -fPIC -std=c++11" )

#-----------------------------------------------------------#

include_directories (${PROJECT_SOURCE_DIR}/include)

#-----------------------------------------------------------#

add_subdirectory(src)
add_subdirectory(tools)

#-----------------------------------------------------------#

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
   )
   install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/html DESTINATION docs)
endif(DOXYGEN_FOUND)

#-----------------------------------------------------------#
