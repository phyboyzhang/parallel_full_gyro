Add_Library(field_mesh STATIC field_2d_mesh.F90)

Target_Link_Libraries(field_mesh utilities meshes)
