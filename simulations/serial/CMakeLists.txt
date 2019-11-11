ADD_Library(start_simul STATIC orbit_storing.F90)
SET(link_source interpolation meshes field_mesh orbit)
TARGET_LINK_LIBRARIES(start_simul ${link_source})

ADD_EXECUTABLE(test_large_epotential test_large_epotential.F90)

set(test_source meshes gyroaverage field_mesh orbit start_simul field_test)

TARGET_LINK_LIBRARIES(test_large_epotential ${test_source})
