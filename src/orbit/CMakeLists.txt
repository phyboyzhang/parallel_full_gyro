SET(LINK_SOURCE utilities interpolation gyroaverage field_mesh piclayout
                meshes parautilities picutilities)

SET(ORBIT_SOURCE  boris_rk_orbit.F90 orbit_2d_2v.F90 m_para_orbit.F90
                  m_tp_para_orbit.F90)

ADD_LIBRARY(orbit STATIC ${ORBIT_SOURCE} )
TARGET_LINK_LIBRARIES(orbit ${LINK_SOURCE})

ADD_SUBDIRECTORY(testing)
