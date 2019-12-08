ADD_EXECUTABLE(test_ful_integrated_simulation test_ful_integrated_simulation.F90)
Target_Link_Libraries(test_ful_integrated_simulation fieldsolver moveparticles para_random_sample picutilities 
                      meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata
                      paradata_type paradata_layout)
ADD_SUBDIRECTORY(single_particle)

ADD_SUBDIRECTORY(integ_sim)
