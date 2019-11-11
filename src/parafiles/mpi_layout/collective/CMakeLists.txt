ADD_LIBRARY(collective STATIC m_collective.F90)

TARGET_LINK_LIBRARIES(collective
                     utilities
                     ${MPI_Fortran_LIBRARIES})
		     
