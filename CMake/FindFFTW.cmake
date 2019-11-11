# FFTW_INCLUDE_DIR = fftw3.f03
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found


SET(FFTW3_DIR /usr/local/Cellar/fftw/3.3.7_1/include)
SET(DFFTPACK_DIR /usr/local/lib)
SET(TRIAL_PATHS ${FFTW3_DIR}
                ${DFFTPACK_DIR}
 )


FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3.f 
                          PATHS ${TRIAL_PATHS})
                          #  PATH_SUFFIXES lib lib64)
			  
SET(FFTW_LIBRARIES ${FFTW_LIBRARY})

FIND_LIBRARY(DFFTPACK_LIBRARIES NAMES libdfftpack.a
                                PATHS ${TRIAL_PATHS})


                              #ENDIF(FFTW_FOUND)
