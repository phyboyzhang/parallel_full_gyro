ADD_EXECUTABLE(test_rksolver test_rksolver.F90)

TARGET_LINK_LIBRARIES(test_rksolver orbit)

ADD_EXECUTABLE(test_interp_elef test_interp_elef.F90)
Target_Link_Libraries(test_interp_elef orbit interpolation)

ADD_EXECUTABLE(test_para_orbit test_para_orbit.F90)
Target_Link_Libraries(test_para_orbit fieldsolver moveparticles
                      orbit para_random_sample picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata 
                      paradata_type paradata_layout)

ADD_EXECUTABLE(test_para_orbit_one test_para_orbit_one.F90)
Target_Link_Libraries(test_para_orbit_one fieldsolver moveparticles
                      orbit para_random_sample picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata
                      paradata_type paradata_layout)
