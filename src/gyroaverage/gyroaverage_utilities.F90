module gyroaverage_utilities
#include "work_precision.h"
use constants, only: pi_
implicit none

public :: s_compute_shape_circle
contains
  subroutine s_compute_shape_circle(points,N_points)
    int4,intent(in) :: N_points
    real8,dimension(:,:),pointer :: points
    int4 :: i
    real8 :: x
    do i=1,N_points
      x = 2.0d0*pi_*real(i,8)/(real(N_points,8))
      points(1,i) = cos(x)
      points(2,i) = sin(x)
      points(3,i) = 1.0d0/real(N_points,8)
    enddo
   
  end subroutine s_compute_shape_circle



end module gyroaverage_utilities


