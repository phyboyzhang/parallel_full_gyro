Add_Executable(test_precomputing test_precomputing.F90)

Target_Link_Libraries(test_precomputing meshes interpolation
                 parautilities piclayout mpilayout gyroaverage paradata_type paradata_layout fieldsolver
                 precompute field_test picutilities para_random_sample diagnosis)
