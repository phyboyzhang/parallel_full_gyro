Add_Library(precompute STATIC m_precompute.F90)

TARGET_LINK_LIBRARIES(precompute interpolation piclayout parautilities paradata_type gyroaverage)
