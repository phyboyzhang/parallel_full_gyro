Add_Library(parautilities STATIC  m_parautilities.F90 paradata_utilities.F90)

SET(LINK_SOURCE  piclayout collective utilities ${MPI_Fortran_LIBRARIES} mpilayout paradata_type)
Target_link_Libraries(parautilities ${LINK_SOURCE})

#ADD_SUBDIRECTORY(testing)
