ADD_EXECUTABLE(orbit_comparison orbit_comparison.F90)
Target_Link_Libraries(orbit_comparison fieldsolver moveparticles
                      orbit para_random_sample picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata constants
                      field_test paradata_layout paradata_type)
Add_EXECUTABLE(test_para_interpo3 test_para_interpo3.F90)
Target_Link_Libraries(test_para_interpo3 fieldsolver moveparticles
                      orbit para_random_sample picutilities meshes interpolation
                      parautilities piclayout mpilayout gyroaverage precompute outputdata constants
                      field_test paradata_layout paradata_type)
