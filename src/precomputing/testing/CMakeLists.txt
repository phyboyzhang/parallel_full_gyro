Add_Executable(test_precomputing test_precomputing.F90)

Target_Link_Libraries(test_precomputing meshes interpolation
                 parautilities piclayout mpilayout gyroaverage paradata_type paradata_layout fieldsolver
                 precompute field_test picutilities para_random_sample diagnosis)

ADD_EXECUTABLE(test_dgetri test_dgetri.F90)

Target_Link_Libraries(test_dgetri ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

