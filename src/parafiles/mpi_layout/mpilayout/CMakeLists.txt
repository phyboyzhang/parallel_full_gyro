IF(PGI)

  FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python)
  FIND_PROGRAM(FPP_EXECUTABLE NAMES gfortran)
  SET(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")
  
  SET(_FILE remapper)
  MESSAGE(STATUS "Cut long lines of file ${_FILE}.F90")

  EXECUTE_PROCESS(COMMAND ${FPP_EXECUTABLE} "-Iinclude" "-I${CMAKE_CURRENT_SOURCE_DIR}" "-E" "-w" "-P" "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}.F90"
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}_pgi.F90")

  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT} 
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  ADD_LIBRARY(remapper STATIC remapper_pgi.F90 )

ELSE(PGI)

   Add_Library(mpilayout STATIC m_mpilayout.F90)

ENDIF(PGI)




Target_Link_Libraries(mpilayout
                      collective
		      utilities
                      piclayout
          ${MPI_Fortran_LIBRARIES}
		     )

ADD_SUBDIRECTORY(testing)
