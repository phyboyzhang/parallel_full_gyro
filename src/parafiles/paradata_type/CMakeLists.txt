Add_Library(paradata_type STATIC paradata_type.F90)

Target_Link_Libraries(paradata_type 
                      collective
                      utilities
                      piclayout
          ${MPI_Fortran_LIBRARIES}
                      )

#Add_Library(paradata_layout STATIC paradata_layout.F90)
#Target_Link_Libraries(paradata_layout
#                      meshes
#                      utilities
#                      collective
#                      utilities
#                      piclayout
#          ${MPI_Fortran_LIBRARIES}
#                      parautilities)
