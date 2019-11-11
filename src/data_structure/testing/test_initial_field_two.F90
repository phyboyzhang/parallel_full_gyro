program test_initial_field_two
#include "work_precision.h"
  use cartesian_mesh, only: cartesian_mesh_1d, &
       f_new_cartesian_mesh_1d
  
  use field_2d_mesh,  only: &
       field_2d_plan, &
       initialize_field_2d_plan
  implicit none

  int4 :: num_cell=4
  real8 :: eta_min=0.0_f64, eta_max=4.0_f64
  real8 :: amp=1.0_f64, wave_one=1.0_f64, wave_two=1.0_f64
  class(cartesian_mesh_1d), pointer :: m_x1
  class(cartesian_mesh_1d), pointer :: m_x2
  class(field_2d_plan),pointer :: field_2d
  int4 :: flag = 0
  character(99) :: geometry="periodic", boundary

  select case(geometry)
  case ("periodic")
     boundary= "periodic"
  case ("polar")
     boundary= "natural"
  case default
     stop
  end select

!!$  allocate(m_x1)
!!$  allocate(m_x2)
!!$  allocate(field_2d)
  
  m_x1=>f_new_cartesian_mesh_1d(num_cell,eta_min,eta_max, boundary)
  m_x2=>f_new_cartesian_mesh_1d(num_cell,eta_min,eta_max, "periodic")

  field_2d => initialize_field_2d_plan(&
       m_x1,  &
       m_x2,  &
       geometry)

 print*, "epotential_init", field_2d%epotential_init
 print*, "ep_weight_init", field_2d%ep_weight_init


end program test_initial_field_two
