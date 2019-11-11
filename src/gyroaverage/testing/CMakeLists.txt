ADD_EXECUTABLE(test_lapack_blas test_lapack_blas.F90)
ADD_EXECUTABLE(test_gyroaverage test_gyroaverage.F90)
#ADD_EXECUTABLE(test_suspect_two test_suspect_two.F90)
ADD_EXECUTABLE(test_para_gyroaverage test_para_gyroaverage.F90)
ADD_EXECUTABLE(test_first_gyroaverage test_first_gyroaverage.F90)



TARGET_LINK_LIBRARIES(test_lapack_blas gyroaverage )
   #                                  ${LAPACK_LIBRARIES}
#				     ${BLAS_LIBRARIES})

TARGET_LINK_LIBRARIES(test_gyroaverage interpolation
                                       meshes				       
                                       gyroaverage
                                       field_mesh
                                       field_test)

#Target_Link_Libraries(test_suspect_two interpolation)

Target_Link_Libraries(test_para_gyroaverage meshes interpolation
                 parautilities piclayout mpilayout gyroaverage paradata_type paradata_layout)

Target_Link_Libraries(test_first_gyroaverage meshes interpolation
                 parautilities piclayout mpilayout gyroaverage paradata_type paradata_layout)
