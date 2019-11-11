program test_initial_field
#include "work_precision.h"

!  use spline_module
  use cartesian_mesh, only: cartesian_mesh_1d, &
                            f_new_cartesian_mesh_1d
  use  field_2d_mesh, only: field_2d_plan, &
       initialize_field_2d_plan

  
  implicit none

  class(cartesian_mesh_1d),pointer :: m_x1,m_x2
  class(field_2d_plan),pointer :: field_2d
!  class(gyroaverage_2d_plan), pointer :: gyro
  int4 :: num_cells=20
  character(99) :: boundary = "double_per", geometry = "cartesian"
  int4 :: N_points=100
  real8 :: contr = 4.0_f64
  real8 :: eta_min(2), eta_max(2),delta_eta(2)
  real8 :: eta_delta(2)
  int4  :: NC(2)
  real8 :: deri_firstorder(2)
  real8 :: x(2)
  int4 :: flag=0
  int4 :: i, j
  
  real8 :: wam=1.0e-3,wn_one=1.0,wn_two=1.0

  real8, dimension(:,:), pointer :: epotential_init, ep_weight_init



  allocate(epotential_init(num_cells,num_cells))
  allocate(ep_weight_init(num_cells,num_cells))
  
 x=(/3.0,3.1/)
  eta_min(1:2) = 0.0
  eta_max(1:2) = 4.0
  delta_eta(:)=eta_max(:)/real(num_cells,8)
  NC(1:2)=num_cells

  m_x1=>f_new_cartesian_mesh_1d(num_cells,eta_min(1),eta_max(1),"periodic")
  m_x2=>f_new_cartesian_mesh_1d(num_cells,eta_min(2),eta_max(2),"periodic")

  field_2d=>initialize_field_2d_plan( &
       m_x1,        &
       m_x2,        &
       trim(geometry))  


end program test_initial_field
