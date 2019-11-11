Add_Library(interpolation STATIC spline_module.F90 m_para_spline.F90)
Target_Link_Libraries(interpolation meshes parautilities utilities piclayout ${MPI_Fortran_LIBRARIES})

ADD_SUBDIRECTORY(testing)
