  IF($ENV{HOSTNAME} MATCHES "hydra*")
      SET(MPI_Fortran_COMPILER "mpiifort")
   ENDIF()
   
   IF(MPI_ENABLED)
      FIND_PACKAGE(MPI REQUIRED Fortran)
   ENDIF(MPI_ENABLED)
   
IF(MPI_Fortran_FOUND)
     MESSAGE(STATUS "MPI FOUND")
     MESSAGE(STATUS "MPI Fortran Compiler is ${MPI_Fortran_COMPILER}")
     MESSAGE(STATUS "MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")
 #   find_path(MPI_Fortran_MOD_DIR NAMES mpi.mod
 #          PATHS $ENV{MPI_FORTRAN_MOD_DIR} ${MPI_Fortran_INCLUDE_PATH})
 #    if(MPI_Fortran_MOD_DIR)
 #      SET(MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_MOD_DIR} ${MPI_Fortran_INCLUDE_PATH})
 #    endif(MPI_Fortran_MOD_DIR)
   ELSE(MPI_FOUND)
     MESSAGE(STATUS "MPI NOT FOUND")
     SET(MPI_ENABLED OFF CACHE BOOL " " FORCE)
   ENDIF(MPI_Fortran_FOUND)
   include_directories(${MPI_Fortran_INCLUDE_PATH})
   #MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY MPI_LIBRARY MPI_Fortran_MOD_DIR)
  
  #MARK_AS_ADVANCED(CLEAR MPI_Fortran_COMPILER)
