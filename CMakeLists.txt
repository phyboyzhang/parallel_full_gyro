
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
PROJECT(ELECTROSTATIC_EXP)
ENABLE_LANGUAGE(Fortran)


set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_CC_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

OPTION(OPENMP_ENABLED       "Enable use of OpenMP library" on)
OPTION(MPI_ENABLED          "Enable use of MPI library"  ON)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH} )

SET(CMAKE_FORTRAN_FLAGS "-llpack -lblas -ldfftpack")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Define folder to be deleted by make clean
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/include/) 
ADD_CUSTOM_TARGET( clear
COMMAND find . -name *.mod -delete
COMMAND find . -name *.a -delete
COMMAND find . -name *.o -delete
COMMAND find . -name *.cmake -delete
COMMAND find . -name *.xmf -delete
COMMAND find . -name *.h5 -delete
COMMAND find . -name *.gnu -delete
COMMAND find . -name *.dat -delete
COMMAND rm -rf bin
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)


ADD_SUBDIRECTORY( external )

find_package(LAPACK)
find_package(BLAS)
find_library(FFTW3_LIBRARIES fftw3)
find_library(DFFTPACK_LIBRARIES dfftpack)

if(NOT FFTW3_LIBRARIES-NOTFOUND)
  MESSAGE(STATUS "FFTW3 LIBRARIES IS FOUND")
endif()

if(NOT DFFTPACK_LIBRARIES-NOTFOUND)
  MESSAGE(STATUS "DFFTPACK LIBRARIES IS FOUND")
  MESSAGE(STATUS "THE PATH OF DFFTPACK IS ${DFFTPACK_LIBRARIES}")
endif()

FIND_PACKAGE(OpenMP)
#include_directories(OpenMP_Fortran_dir)
if( OpenMP_FOUND )
MESSAGE(STATUS "The path of OpenMP Fortran is ${OpenMP_Fortran_dir}")
#SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
MESSAGE(STATUS "OpenMP fortran flags is ${OpenMP_Fortran_FLAGS}")
endif()
 
INCLUDE(MPIConfig)


FILE(GLOB_RECURSE MACROS ${PROJECT_SOURCE_DIR}/src/*.h)
FILE(COPY ${MACROS} DESTINATION ${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)



include(CTest)
enable_testing()



ADD_SUBDIRECTORY( src )

file(GLOB TARGETS "*.EXE *.dll")
set_target_properties(
  ${GLOB}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  )

ADD_SUBDIRECTORY( simulations )
