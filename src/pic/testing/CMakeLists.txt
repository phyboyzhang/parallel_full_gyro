ADD_EXECUTABLE(test_para_pic test_para_pic.F90)

target_link_libraries(test_para_pic fieldsolver moveparticles para_random_sample picutilities meshes interpolation
                  parautilities piclayout mpilayout gyroaverage precompute
                  paradata_type paradata_layout)

ADD_EXECUTABLE(test_para_partgrid_coupl test_para_partgrid_coupl.F90)
Target_Link_Libraries(test_para_partgrid_coupl fieldsolver moveparticles para_random_sample picutilities 
                      meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata
                      paradata_type paradata_layout)


ADD_EXECUTABLE(test_para_partgrid_coupl_one test_para_partgrid_coupl_one.F90)
Target_Link_Libraries(test_para_partgrid_coupl_one fieldsolver moveparticles para_random_sample picutilities 
                      meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata
                      paradata_type paradata_layout)
