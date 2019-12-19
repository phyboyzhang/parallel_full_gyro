Add_Library(para_random_sample STATIC para_random_sample.F90 congruence_sampling.F90)

Target_Link_Libraries(para_random_sample utilities piclayout constants
                    picutilities moveparticles mpilayout parautilities paradata_type field_test)

ADD_SUBDIRECTORY(testing)
