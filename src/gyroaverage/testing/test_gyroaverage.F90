program test_gyroaverage
#include "work_precision.h"
  use spline_module, only: &
       compute_spl2d_double_per_weight, &
       compute_spl2d_nat_per_weight, &
       compute_spl2d_field_point_per_per, &
       compute_spl2d_firstorder_derivative_point_per_per
  
  use cartesian_mesh, only: cartesian_mesh_1d, &
                            f_new_cartesian_mesh_1d
  use  field_2d_mesh, only: field_2d_plan, &
       initialize_field_2d_plan
     
!       compute_spl2D_field_weight_2d_plan
  use field_initialize, only: &
      compute_field_2d_mesh

  use gyroaverage_2d, only: &
      gyroaverage_2d_plan, &
      initialize_gyroaverage_2d_plan,   &
      compute_gyroaverage_2d_point_double_per
  
!!$  use field_initialize, only: &
!!$      epotential_2d_analytical, &
!!$       compute_field_2d_mesh
  
  implicit none

  class(cartesian_mesh_1d),pointer :: m_x1,m_x2
  class(field_2d_plan),pointer :: field_2d
  class(gyroaverage_2d_plan), pointer :: gyro
  int4 :: num_cells=50
  character(99) :: boundary = "double_per", geometry = "cartesian"
  int4 :: N_points=100
  real8 :: mu=0.2
  real8 :: contr = 4.0_f64
  real8 :: eta_min(2), eta_max(2),delta_eta(2)
  real8 :: eta_delta(2)
  int4  :: NC(2)
  real8 :: deri_firstorder(2)
  real8 :: x(2)
  int4 :: flag=0
  int4 :: i, j
  real8,dimension(:,:), allocatable :: vec
  
  real8 :: wam=1.0e-3,wn_one=1.0,wn_two=1.0

  real8, dimension(:,:), pointer :: epotential_init, ep_weight_init


  allocate(vec(num_cells,3))
  allocate(epotential_init(num_cells,num_cells))
  allocate(ep_weight_init(num_cells,num_cells))
  
 x=(/3.0,3.1/)
  eta_min(1:2) = 0.0
  eta_max(1:2) = 8.0
  delta_eta(:)=eta_max(:)/real(num_cells,8)
  NC(1:2)=num_cells


  do i=1,num_cells
    vec(i,1)=3.0
    vec(i,2)=eta_min(2)+real(i-1,8)*delta_eta(2)  
  enddo

  m_x1=>f_new_cartesian_mesh_1d(num_cells,eta_min(1),eta_max(1),"periodic")
  m_x2=>f_new_cartesian_mesh_1d(num_cells,eta_min(2),eta_max(2),"periodic")
  field_2d=>initialize_field_2d_plan( &
       m_x1,        &
       m_x2,        &
       trim(geometry))

  call  compute_field_2d_mesh( &
       field_2d, &
       m_x1,  &
       m_x2,  &
       wam,     &
       wn_one,  &
       wn_two,  &
       contr,   &
       geometry)
!!$
!!$  do i=1,num_cells
!!$     do j=1,num_cells
!!$  !      field_2d%epotential_init(i,j)=1.0
!!$  !      field_2d%bf_init_3rd(i,j)=1.0
!!$        epotential_init(i,j)=1.0
!!$        
!!$     end do
!!$  end do
!!$  
!!$  call compute_spl2D_double_per_weight( &
!!$  !     field_2d%epotential_init, &
!!$       !     field_2d%ep_weight_init,  &
!!$       epotential_init, &
!!$       ep_weight_init, &
!!$       num_cells, &
!!$       num_cells)
!!$

  call compute_spl2D_double_per_weight(&
       field_2d%epotential_init, &
       field_2d%ep_weight_init,  &
       m_x1%num_cells,  &
       m_x2%num_cells)

  call compute_spl2D_double_per_weight(&
       field_2d%Bf_init_3rd, &
       field_2d%Bf_weight_3rd,  &
       m_x1%num_cells,  &
       m_x2%num_cells)  

  gyro => initialize_gyroaverage_2d_plan(N_points,mu,m_x1,m_x2,field_2d,geometry,boundary)

!  do i=1, num_cells
!    print*, "i=",i,field_2d%ep_weight_init(:,i)
!  end do

!  do i=1,num_cells
!    print*, "i=",i,gyro%ep_weight_gyro(:,i)
!  end do

! do i=1,num_cell
!    print*, "i=",i,field_2d%epotential_init(i,:)
!  end do

!  do i=1,num_cells
!   vec(i,3)=compute_gyroaverage_2d_point_double_per(vec(i,1:2),m_x1,m_x2,gyro,field_2d,mu)
!   vec(i,4)= compute_spl2d_field_point_per_per(vec(i,1:2), &
!       field_2d%ep_weight_init,  &
!       eta_min,    &
!       eta_max,    &
!       delta_eta,  &
!       NC)
!   print*, vec(i,3), vec(i,4)
!  end do



!!$    call compute_spl2d_firstorder_derivatve_point_per_per(x,&
!!$       field_2d%ep_weight_init, &
!!$       eta_min, &
!!$       eta_max, &
!!$       delta_eta, &
!!$       NC,      &
!!$       deri_firstorder)
!!$
!!$       print*, "no gyroaverage", deri_firstorder
!!$  
  do i=1, num_cells
    call compute_spl2d_firstorder_derivative_point_per_per(vec(i,1:2),&
       gyro%ep_weight_gyro, &
       eta_min, &
       eta_max, &
       delta_eta, &
       NC,      &
       deri_firstorder)
    print*, "i=",i,deri_firstorder
  end do
!!$
!!$    print*, deri_firstorder

end program test_gyroaverage
