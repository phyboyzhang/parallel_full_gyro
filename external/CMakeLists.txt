# Allow compiler-specific intrinsic functions (non-standard)
IF (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  SET( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -fall-intrinsics" )
  SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fall-intrinsics" )
ENDIF()

ADD_SUBDIRECTORY(fftpack)
#ADD_SUBDIRECTORY(fftpack5)
#ADD_SUBDIRECTORY(pppack)
#ADD_SUBDIRECTORY(nufft)
#ADD_SUBDIRECTORY(burkardt)
#ADD_SUBDIRECTORY(clapp)

#IF(BUILD_MUDPACK)
#  ADD_SUBDIRECTORY(mudpack)
#ENDIF(BUILD_MUDPACK)

INCLUDE(ExternalProject)

IF(BUILD_HDF5)

  FIND_PACKAGE(MPI REQUIRED CXX)
  EXTERNALPROJECT_ADD( hdf5
     URL https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.bz2
     SOURCE_DIR ${CMAKE_BINARY_DIR}/hdf5
     BINARY_DIR ${CMAKE_BINARY_DIR}/hdf5
     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/hdf5/configure --enable-fortran --prefix=${CMAKE_BINARY_DIR} --enable-parallel CXX=${MPI_CXX_COMPILER} CC=${MPI_C_COMPILER} FC=${MPI_Fortran_COMPILER}
  )

ENDIF(BUILD_HDF5)


IF(BUILD_PARALUTION)
  INCLUDE(ExternalProject)
  SET(PARALUTION_VERSION 1.1.0 CACHE STRING "Paralution version number")
  EXTERNALPROJECT_ADD( paralution
     URL  http://www.paralution.com/downloads/paralution-${PARALUTION_VERSION}.tar.gz
     CMAKE_ARGS -DSUPPORT_CUDA=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
     INSTALL_COMMAND ""
  )
ENDIF()

IF($ENV{PFUNIT})
  INCLUDE(ExternalProject)
  EXTERNALPROJECT_ADD(
    pFUnit
    PREFIX pFunit
    GIT_REPOSITORY git://git.code.sf.net/p/pfunit/code
    CMAKE_ARGS -DINSTALL_PATH=${CMAKE_BINARY_DIR})
ENDIF()

IF(BUILD_SPRNG)
  INCLUDE(BuildSPRNG)
  ADD_SUBDIRECTORY(sprng)
ENDIF()
