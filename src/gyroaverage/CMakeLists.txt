ADD_LIBRARY(gyroaverage STATIC gyroaverage_2d_base.F90
                               gyroaverage_2d.F90
                               gyroaverage_utilities.F90 
                               para_gyroaverage_2d_one.F90)
# Add_Library(gyroaverage STATIC gyroaverage_cartesian_2d.F90)
set(gyro_source meshes 
              interpolation
              ${BLAS_LIBRARIES}
              ${LAPACK_LIBRARIES} 
              ${DFFTW_LIBRARIES}
              dfftpack 
              field_mesh
              piclayout
              mpilayout
              parautilities)
TARGET_LINK_LIBRARIES( gyroaverage ${gyro_source})
#IF(BUILD_TESTING)
ADD_SUBDIRECTORY(testing)
#ENDIF(BUILD_TESTING)
