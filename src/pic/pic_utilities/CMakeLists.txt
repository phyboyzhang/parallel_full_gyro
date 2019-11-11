Add_Library(picutilities STATIC m_picutilities.F90)

Target_Link_Libraries(picutilities utilities mpilayout piclayout parautilities paradata_type)

Add_Subdirectory(testing)
