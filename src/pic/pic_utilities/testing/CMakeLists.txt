Add_executable(test_pic_utilities test_pic_utilities.F90)

Target_Link_Libraries(test_pic_utilities picutilities meshes interpolation
                 parautilities piclayout mpilayout gyroaverage paradata_type paradata_layout)
