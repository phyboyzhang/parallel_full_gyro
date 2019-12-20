ADD_EXECUTABLE(integ_per_per integ_per_per.F90)
Target_Link_Libraries(integ_per_per fieldsolver moveparticles
                      orbit para_random_sample picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata constants
                      field_test paradata_layout paradata_type diagnosis tp_preparation)

#ADD_EXECUTABLE(test_random_sample test_random_sample.F90)
#Target_link_libraries(test_random_sample para_random_sample moveparticles
#                      orbit picutilities meshes interpolation
#                      parautilities piclayout mpilayout gyroaverage precompute outputdata constants
#                      field_test paradata_layout paradata_type)
