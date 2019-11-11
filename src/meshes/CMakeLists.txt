SET( MESH_SOURCES 
    abstract_mesh_base.F90 
    cartesian_mesh.F90)

ADD_LIBRARY( meshes STATIC  ${MESH_SOURCES} )

IF(BUILD_TESTING)
ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
