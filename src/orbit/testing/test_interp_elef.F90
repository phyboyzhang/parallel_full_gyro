program test_interp_elef
#include "work_precision.h"
use boris_rk_orbit, only: obtain_interpolation_elefield_per_per
use spline_module
use cartesian_mesh, only: cartesian_mesh_1d
implicit none

 int4 :: num_cell=50
  real8,dimension(:,:), pointer :: field_init
  real8,dimension(:,:), pointer :: mesh_field_weight
  real8,dimension(:,:), pointer :: mesh_field_weight_1
  int4 :: i,j
  real8 :: eta_min(2), eta_max(2)
  real8 :: eta_delta(2)
  int4  :: NC(2)
  real8 :: deri_firstorder(2),deri_firstorder_1(2)
  real8 :: x(2),elef(3)


  x=(/2.0,3.1/)
!  eta_min = 0.0
!  eta_max = 2.0

  eta_min=(/0.0,0.0/)
  eta_max=(/2.0,2.0/)

  eta_delta=(/(eta_max(1)-eta_min(1))/real(num_cell,8),(eta_max(2)-eta_min(2))/real(num_cell,8)/)


  NC=num_cell
  allocate(field_init(num_cell,num_cell))
  allocate(mesh_field_weight_1(num_cell,num_cell),mesh_field_weight(num_cell,num_cell))
!  field_init= reshape((/1.0, 2.0, 3.0, 4.0, &
!                   1.1, 2.1, 3.1, 4.1, &
!                   1.2, 2.2, 3.2, 4.2, &
!                   1.3, 2.3, 3.3, 4.3/),(/4,4/))

!  field_init= reshape((/1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0, &
!                   1.0, 1.0, 1.0, 1.0/),(/4,4/))  

!  field_init= reshape((/1.0, 2.0, 3.0, 4.0, &
!                   1.0, 2.0, 3.0, 4.0, &
!                   1.0, 2.0, 3.0, 4.0, &
!                   1.0, 2.0, 3.0, 4.0, &
!                   1.0, 2.0, 3.0, 4.0/),(/4,4/))

       x=0.0_f64
       do i=0,num_cell-1
         x(1)=x(1)+ eta_delta(1)
         x(2)=0.0_f64
          do j=0,num_cell-1
            x(2)=x(2)+eta_delta(2)
           field_init(i+1,j+1) = x(2)
          end do
       end do

!print*, field_init(10,:) 

 mesh_field_weight=0._f64
 call compute_spl2d_double_per_weight(&
      field_init, &
      mesh_field_weight,&
      num_cell, &
      num_cell)

!  do i=1, num_cell
!    print*,"per_per", mesh_field_weight(i,:)
!  end do
! print*, "+++++++++++++++"
 call compute_spl2d_nat_per_weight(&
      field_init, &
      mesh_field_weight_1,&
      num_cell-1, &
      num_cell)

!  do i=1, num_cell
!    print*, "nat_per",mesh_field_weight_1(i,:)
!  end do


call obtain_interpolation_elefield_per_per(x,    &
                                             elef,  &
                                             mesh_field_weight_1,  &
                                             eta_min,    &
                                             eta_max,    &
                                             eta_delta,  &
                                             NC)  

print*, "elef=", elef

  do i=1, num_cell
    print*, "i=",i,mesh_field_weight_1(i,:)
  end do

!print*, "mesh_field_weight_1", mesh_field_weight_1
end program
