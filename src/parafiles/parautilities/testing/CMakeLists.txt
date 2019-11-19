ADD_EXECUTABLE(test_alltoallv_per_per test_alltoallv_per_per.F90)
TARGET_LINK_LIBRARIES(test_alltoallv_per_per 
                      mpilayout utilities 
                      parautilities
                      paradata_type paradata_layout
                      )
ADD_EXECUTABLE(test_alltoallv_subroutine test_alltoallv_subroutine.F90)
TARGET_LINK_LIBRARIES(test_alltoallv_subroutine
                      parautilities)

ADD_EXECUTABLE(test_paradata_utilities test_paradata_utilities.F90)
TARGET_LINK_LIBRARIES(test_paradata_utilities
                      mpilayout utilities
                      parautilities
                      piclayout  
                      meshes  paradata_type paradata_layout                  
                      )

ADD_EXECUTABLE(test_point_process test_point_process.F90)
TARGET_LINK_LIBRARIES(test_point_process
                      mpilayout utilities
                      parautilities
                      piclayout
                      meshes  paradata_type paradata_layout
                      )
