program test_lapack_blas
#include "work_precision.h"
use gyroaverage_2d, only: inverse_cubic_splines_matrix, & 
                                    matrix_solve_polar_splines
! inverse_cubic_splines_matrix
!use lib, only: s_inverse_cubic_splines_matrix
implicit none
  
  comp8,dimension(1,3,3) :: mat, mat_inverse
  int4 :: i=3,j
  real8, dimension(3,3) :: phi_in,phi_out  
  mat(1,:,:)= reshape((/(1.0_f64,1.0_f64),(2.0d0,2.0d0),(3.0d0,3.0d0), &
       (4.0d0, 4.0d0),(5.0d0,5.0d0),(6.0d0,6.0d0), &
       (7.0d0, 7.0d0),(8.0d0,8.0d0),(9.0d0,9.0d0)/),shape(mat(1,:,:)))
  phi_in(:,:)=reshape((/1.0d0,2.0d0,3.0d0, &
                        1.0d0,2.0d0,3.0d0, &
                        1.0d0,2.0d0,3.0d0/),shape(phi_in(:,:)))
  mat_inverse=mat

 !call inverse_cubic_splines_matrix(mat_inverse,i-1,1)
  call matrix_solve_polar_splines(phi_in,&
                                 phi_out,&
                                 mat_inverse, &
                                 i-1, &
                                 i)
                                 
!  mat(1,:,:)=matmul(mat(1,:,:),mat_inverse(1,:,:))

  do j=1,3
     print*, phi_out(j,:)
  end do

! print*,phi_out(1,:)
end program test_lapack_blas
