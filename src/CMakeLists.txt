cmake_minimum_required (VERSION 2.6)
project (freestream)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

find_package (GSL REQUIRED)

set (CMAKE_CXX_FLAGS "${GSL_CXX_FLAGS} -g -O3 -Wall")

SET(INCLUDE_DIRS ${GSL_INCLUDE_DIR})
SET(LIBS ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

set (SOURCES
  FreestreamMilne.cpp
  )

add_library(freestream_lib SHARED ${SOURCES})
#add_executable(freestream.e ${SOURCES})
target_link_libraries (freestream_lib ${LIBS})
#target_link_libraries (freestream.e ${LIBS})

#install(TARGETS freestream.e DESTINATION ${CMAKE_HOME_DIRECTORY})
