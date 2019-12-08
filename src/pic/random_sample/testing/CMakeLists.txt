ADD_EXECUTABLE(random_sample random_sample.F90)

TARGET_LINK_LIBRARIES(random_sample utilities meshes)

ADD_EXECUTABLE(test_para_random_sample test_para_random_sample.F90)
Target_link_libraries(test_para_random_sample para_random_sample moveparticles
                      orbit picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata constants
                      field_test fieldsolver paradata_layout paradata_type diagnosis)

ADD_EXECUTABLE(random_sample_2d random_sample_2d.F90)

TARGET_LINK_LIBRARIES(random_sample_2d utilities meshes)

ADD_EXECUTABLE(random_sample_mu random_sample_mu.F90)
target_link_libraries(random_sample_mu)


ADD_EXECUTABLE(test_pic_noise test_pic_noise.F90)
Target_Link_Libraries(test_pic_noise piclayout meshes picutilities)
