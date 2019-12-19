#ADD_LIBRARY(field_solver STATIC field_2d_mesh.F90)
#TARGET_LINK_LIBRARIES(field_solver interpolation meshes)

#ADD_SUBDIRECTORY(qne_solver )

ADD_LIBRARY(fieldsolver STATIC m_fieldsolver.F90)
Target_Link_Libraries(fieldsolver mpilayout interpolation piclayout
                                  parautilities paradata_type field_test)
