# Set correct compilation flags for FFTPACK, depending on Fortran compiler
IF(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   set(CMAKE_Fortran_FLAGS_DEBUG    "-g -w -cpp -fdefault-real-8 ")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -w -cpp -fdefault-real-8 ")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -r8 -nowarn")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-r8 -nowarn -O3 -xHost -ip -qopenmp")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES IBM)
   SET(CMAKE_Fortran_COMPILER xlf)
   SET(CMAKE_Fortran_FLAGS_DEBUG   "-qfixed -qrealsize=8")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-qfixed -qrealsize=8")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 ")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-fast -r8 ")
ENDIF()

# New version of FFTPACK (this will substitute old 'dfftpack')
ADD_LIBRARY(fftpack5 STATIC fftpack5.f)

# Every external library should be installed (otherwise 'package' does not work)
INSTALL(TARGETS fftpack5 ARCHIVE DESTINATION lib)

# CTest
IF(BUILD_TESTING)
  ADD_EXECUTABLE(test_fftpack5 fftpack5_prb.f)
  TARGET_LINK_LIBRARIES(test_fftpack5 fftpack5)
ENDIF(BUILD_TESTING)
