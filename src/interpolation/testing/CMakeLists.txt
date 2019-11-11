ADD_EXECUTABLE(test_interpolation test_interpolation.F90)
ADD_EXECUTABLE(test_interpo_gradient test_interpo_gradient.F90)
TARGET_LINK_LIBRARIES(test_interpolation interpolation)
TARGET_LINK_LIBRARIES(test_interpo_gradient
                      interpolation
		       meshes
                       )
     
        #ADD_EXECUTABLE(test_boundary test_boundary.F90)
        #Target_Link_Libraries(test_boundary interpolation)
ADD_EXECUTABLE(test_com_natper test_com_natper.F90)
Target_Link_Libraries(test_com_natper meshes interpolation)

        #ADD_EXECUTABLE(test_interpotwo test_interpotwo.F90)
        #Target_Link_Libraries(test_interpotwo interpolation)

ADD_EXECUTABLE(test_para_interpolation test_para_interpolation.F90)
Target_Link_Libraries(test_para_interpolation meshes interpolation
                 parautilities piclayout mpilayout paradata_type paradata_layout)
