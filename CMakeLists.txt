project(DHECR CXX)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize" )
set( DHECR_INCLUDES )
set( DHECR_LIBRARIES )

# CRPROPA (provided)
add_subdirectory(libs/crpropa)
list(APPEND DHECR_LIBRARIES crpropa)
list(APPEND DHECR_INCLUDES libs/crpropa/include)

# HepPID (provided)
add_subdirectory(libs/HepPID)
list(APPEND DHECR_LIBRARIES HepPID)
list(APPEND DHECR_INCLUDES libs/HepPID/include)

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND DHECR_LIBRARIES kiss)
list(APPEND DHECR_INCLUDES libs/kiss/include)

# FFTW3F (external)
find_package(FFTW3F)
if(FFTW3F_FOUND)
	list(APPEND DHECR_INCLUDES ${FFTW3F_INCLUDE_DIR})
	list(APPEND DHECR_LIBRARIES ${FFTW3F_LIBRARY})
	add_definitions(-DCRPROPA_HAVE_FFTW3F)
endif(FFTW3F_FOUND)

# add source files
set( DHECR_SRCS
     src/DiffusionCoefficients.cpp
     src/Trajectory.cpp
)

# make library
include_directories (include ${DHECR_INCLUDES})
add_library( DHECR_LIB ${DHECR_SRCS} ${DHECR_INCLUDES} )

# make executables
add_executable(dhecr src/main.cpp)
target_link_libraries (dhecr DHECR_LIB ${DHECR_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "DHECR_INCLUDES ...... = ${DHECR_INCLUDES}")
message (STATUS "DHECR_LIBRARIES ..... = ${DHECR_LIBRARIES}")
