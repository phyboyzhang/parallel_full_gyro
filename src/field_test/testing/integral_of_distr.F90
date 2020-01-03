program integral_distr
#include "work_precision.h"
use constants, only: pi_ 
use field_initialize, only: test_equdistr
implicit none

  real8 :: ymax,y(2), delta
  real8 :: integ=0.0,mean(2)
  int4 :: i,j,N=1000
  character(100) :: distr="gaussian"
  real8 :: sigma=25.0

  ymax=2.0*pi_  
  delta=ymax/real(1000,8)
  mean(1)=pi_
 
  do j=1,N
    y(2)=delta*real(j-1,8)
    do i=1,N
      y(1)=delta*real(i-1,8)
      integ=integ+test_equdistr(y,sigma,mean,distr)*delta**2
    enddo
  enddo

  print*, integ
end program integral_distr
