program test_rksolver
#include "work_precision.h"
!  use boris_rk_orbit, only: fulrksolve_noninterpolation

  real8 :: vec(6), elef(3), dtreal, time
  int4  :: flag=0
  int4 ::  i
  vec=(/1.0,0.0,0.0,1.0,0.0,0.0/)
  elef=(/1.0,0.0,0.0/)
  dtreal = 0.1
  time = 0.0


!!$  do i=1,1000
!!$  call fulrksolve_noninterpolation(vec,elef,dtreal,time,flag)
!!$
!!$  print*, vec(:)
!!$
!!$  end do

  

end program test_rksolver
