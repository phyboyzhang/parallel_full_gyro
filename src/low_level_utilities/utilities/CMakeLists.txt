Add_Library(utilities STATIC utilities_module.F90)
Target_Link_Libraries(utilities  constants ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_SUBDIRECTORY(testing)
