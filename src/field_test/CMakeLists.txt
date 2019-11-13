ADD_LIBRARY(field_test STATIC field_initialize.F90)

Target_Link_Libraries(field_test constants meshes field_mesh mpilayout parautilities paradata_type)
