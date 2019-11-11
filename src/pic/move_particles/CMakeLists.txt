Add_Library(moveparticles STATIC m_moveparticles.F90)


SET(link_source orbit ${OpenMP_Fortran_Libraries}  orbit mpilayout piclayout parautilities picutilities paradata_type)

Target_Link_Libraries(moveparticles ${link_source})
