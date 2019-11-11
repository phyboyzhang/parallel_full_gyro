program test_interpolation
#include "work_precision.h"
  use spline_module, only: &
       compute_spl2d_double_per_weight, &
       compute_spl2d_nat_per_weight, &
       compute_spl2d_firstorder_derivative_point_per_per,&
       compute_spl2d_firstorder_derivative_point_nat_per

  implicit none

  int4 :: num_cell=4
  real8,dimension(:,:), pointer :: field_init
  real8,dimension(:,:), pointer :: mesh_field_weight
  real8,dimension(:,:), pointer :: mesh_field_weight_1
  int4 :: i,j
  real8 :: eta_min(2), eta_max(2)
  real8 :: eta_delta(2)
  int4  :: NC(2)
  real8 :: deri_firstorder(2),deri_firstorder_1(2)
  real8 :: x(2)

  x=(/1.0,1.1/)
  eta_min = 0.0
  eta_max = 2.0
  NC=num_cell
  allocate(field_init(num_cell,num_cell))
  allocate(mesh_field_weight_1(num_cell,num_cell), mesh_field_weight(num_cell,num_cell))
!  field_init= reshape((/1.0, 2.0, 3.0, 4.0, &
!                   1.1, 2.1, 3.1, 4.1, &
!                   1.2, 2.2, 3.2, 3.2, &
!                   1.3, 2.3, 3.3, 3.4, &
!                   1.4, 2.4, 3.4, 4.4/),(/4,4/))

!  field_init= reshape((/1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0/),(/4,4/))  

  field_init= reshape((/1.0, 2.0, 3.0, 4.0, &
                   1.0, 2.0, 3.0, 4.0, &
                   1.0, 2.0, 3.0, 4.0, &
                   1.0, 2.0, 3.0, 4.0, &
                   1.0, 2.0, 3.0, 4.0/),(/4,4/))
!

! field_init=1.0_f64
 mesh_field_weight=0._f64 
 call compute_spl2d_double_per_weight(&
      field_init, &
      mesh_field_weight,&
      num_cell, &
      num_cell)

  do i=1, num_cell
    print*,"per_per", mesh_field_weight(i,:)
  end do
 print*, "+++++++++++++++"
 call compute_spl2d_nat_per_weight(&
      field_init, &
      mesh_field_weight_1,&
      num_cell-1, &
      num_cell)

  do i=1, num_cell
    print*, "nat_per",mesh_field_weight_1(i,:)
  end do



! do i=1,num_cell
!    do j=1,num_cell
!      print*, mesh_field_weight_1(i,j)-mesh_field_weight(i,j)
!    end do
! end do
 

!!!!!===>
  do i=1, 2
 ! eta_delta(i)=(eta_max(i)-eta_min(i))/real(num_cell,8)
  end do
  eta_delta(1)=1.0
  eta_delta(2)=0.0
  print*, "eta_max", eta_max(1), eta_max(2)
  print*, "eta_delta", eta_delta(:)
  print*, "x", x(:)
  call compute_spl2d_firstorder_derivative_point_per_per(x,&
       field_init, &
       eta_min, &
       eta_max, &
       eta_delta, &
       NC,      &
       deri_firstorder)

  call compute_spl2d_firstorder_derivative_point_nat_per(x,&
       field_init, &
       eta_min, &
       eta_max, &
       eta_delta, &
       NC,      &
       deri_firstorder_1)

  print*, "deri=",deri_firstorder_1(1), deri_firstorder_1(2)

end program test_interpolation
