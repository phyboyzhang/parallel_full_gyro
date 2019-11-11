program test_large_epotential
#include "work_precision.h"

  use cartesian_mesh, only: cartesian_mesh_1d, &
                               f_new_cartesian_mesh_1d
  
  use gyroaverage_2d, only: gyroaverage_2d_plan, &
          initialize_gyroaverage_2d_plan, &
          compute_spl2d_driftsquare_weight_per_per, &
          compute_spl2d_driftsquare_weight_nat_per
  
  use orbit_storing, only: borisrun, &
                          fulrkorbit, &
                          gyrorkorbit, &
                          borisrun_anal
!!$  
  use field_2d_mesh, only: field_2d_plan, &
       initialize_field_2d_plan
  
  use orbit_2d_2v, only: &
  particle_orbit_cartesian_2d2v, &
  initialize_particle_orbit_cartesian_2d2v, &
  gyro_orbit_cartesian_2d2v, &
  initialize_gyro_orbit_cartesian_2d2v

  use field_initialize, only: &
       compute_field_2d_mesh

  use spline_module, only: &
       compute_spl2D_double_per_weight, &
       compute_spl2D_nat_per_weight

       
  
  implicit none
  
 type simulation_orbit_comparison
  int4  :: nodes_x1, nodes_x2    
  real8 :: x1_min, x1_max, x2_min, x2_max
  real8 :: amplitude, wavenumber_one, wavenumber_two
  int4  :: iter_number
  int4  :: N_par
  int4  :: num_time
  real8 :: dt_gc
  int4  :: gyropoints
  real8 :: mu
  int4  :: flag
  character(99) :: geometry, boundary
  int4   :: orbitorder
  real8 :: contr
  
  class(cartesian_mesh_1d), pointer :: m_x1, m_x2
  class(field_2d_plan), pointer :: field_2d
  class(gyroaverage_2d_plan), pointer :: gyro
  class(particle_orbit_cartesian_2d2v), pointer :: full2d2v
  class(gyro_orbit_cartesian_2d2v), pointer :: gyro2d2v

end type simulation_orbit_comparison

  int4   :: nodes_x1, nodes_x2
  real8 :: x1_min, x1_max, x2_min, x2_max
  real8 :: amplitude, wavenumber_one, wavenumber_two
  int4  :: iter_number
  int4  :: N_par
  int4  :: num_time
  real8 :: dt_gc
  int4  :: gyropoints
  real8 :: mu
  real8 :: x1,x2,x3
  int4 :: flag
  character(99) :: boundary, geometry
  real8 :: contr
  int4  :: orbitorder 

  type(simulation_orbit_comparison) :: sim
  character(len=256) :: filename
  int4, parameter :: input_file = 99
  int4 :: IO_stat
  real8 :: x_pos(3), x4(3),v4(3)
  int4 :: i


  real8 :: eta_min(2),eta_max(2),delta_eta(2),deri_firstorder(2),x(2)
  int4 :: NC(2)

  
!!!=========
!  int4 :: nodes(2), num_cells(2)
!  real8 :: delta_eta(2)
    

  namelist /mesh/    &
       nodes_x1, &
       nodes_x2, &
       x1_min,       &
       x1_max,       &
       x2_min,       &
       x2_max      
  namelist /wave/      &
       amplitude,      &
       wavenumber_one, &
       wavenumber_two 
  namelist /iteration/&
       iter_number,   & 
       dt_gc,         &
       num_time       
  namelist /gyro/     &
       gyropoints,    &
       mu,            &
       orbitorder          
  namelist /init_position_number/  &
       x1,                  &
       x2,                  &
       x3,                  &
       N_par
  namelist /params/  &
       flag,   &
       contr,  &
       boundary,  &
       geometry

  call get_command_argument(1,filename)

  open(unit = input_file, file=trim(filename),IOstat=IO_stat)
  if(IO_stat /= 0) then
     print*, "#it's failed to open the initial data input file"
     stop
  end if

  read(input_file, mesh)
  read(input_file, wave)
  read(input_file, iteration)
  read(input_file, gyro)
  read(input_file, init_position_number)
  read(input_file, params)
  
  close(input_file)

  print*, "#nodes_x1=", nodes_x1
  print*, "#nodes_x2=", nodes_x2
  print*, "#x1_min=",       x1_min      
  print*, "#x1_max=",       x1_max
  print*, "#x2_min=",       x2_min
  print*, "#x2_max=",       x2_max     
  print*, "#amplitude=",    amplitude
  print*, "#wavenumber_one", wavenumber_one
  print*, "#wavenumber_two", wavenumber_two
  print*, "#iter_number",    iter_number
  print*, "#dt_gc=",         dt_gc
  print*, "#num_time=",     num_time
  print*, "#gyropoints=",   gyropoints
  print*, "#mu=",           mu
  print*, "#x1=",           x1
  print*, "#x2=",           x2
  print*, "#x3=",           x3
  print*, "#N_par=",         N_par
  print*, "#flag=",          flag
  print*, "#contr=",         contr
  print*, "#boundary=",      boundary
  print*, "#geometry=",      geometry
  print*, "#orbitorder=",    orbitorder  
  
  sim%nodes_x1=nodes_x1
  sim%nodes_x2=nodes_x2
  sim%x1_min=x1_min
  sim%x1_max=x1_max
  sim%x2_min=x2_min
  sim%x2_max=x2_max
  sim%amplitude=amplitude
  sim%wavenumber_one=wavenumber_one
  sim%wavenumber_two=wavenumber_two
  sim%iter_number=iter_number
  sim%num_time=num_time
  sim%dt_gc=dt_gc
  sim%gyropoints=gyropoints
  sim%flag=flag
  sim%N_par=N_par
  sim%mu=mu
  sim%geometry=geometry
  sim%boundary=boundary
  sim%orbitorder=orbitorder

  sim%m_x1=>f_new_cartesian_mesh_1d(nodes_x1,x1_min,x1_max,"periodic")
  sim%m_x2=>f_new_cartesian_mesh_1d(nodes_x2,x2_min,x2_max,"periodic")

!  print*, "#delta_eta(1),(2)=", sim%m_x1%delta_eta,sim%m_x2%delta_eta
  !!!! allocate the memory of field_2d
  sim%field_2d=>initialize_field_2d_plan( &
       sim%m_x1,        &
       sim%m_x2,        &
       geometry)

  eta_min(1)= sim%m_x1%eta_min
  eta_min(2)= sim%m_x2%eta_min
  eta_max(1)= sim%m_x1%eta_max
  eta_max(2)= sim%m_x2%eta_max
  delta_eta(1)= sim%m_x1%delta_eta
  delta_eta(2)= sim%m_x2%delta_eta
  NC(1)= sim%m_x1%num_cells
  NC(2)= sim%m_x2%num_cells

print*, 1
  !!!! initialize the epotenital mesh and Bfield mesh
  call  compute_field_2d_mesh( &
       sim%field_2d, &
       sim%m_x1,  &
       sim%m_x2,  &
       amplitude,        &
       wavenumber_one,  &
       wavenumber_two,  &
       contr,           &
       geometry)

!!!! compute the steaf matrix of eletrostatic potental and the magnetic field
print*,2

  select case(geometry)
  case("cartesian")
   select case(boundary)
     case ("double_per")
     call compute_spl2D_double_per_weight(&
          sim%field_2d%epotential_init, &
          sim%field_2d%ep_weight_init,  &
          sim%m_x1%num_cells,  &
          sim%m_x2%num_cells)

!print*, sim%field_2d%epotential_init(10,:) 

     call compute_spl2D_double_per_weight(&
          sim%field_2d%Bf_init_3rd, &
          sim%field_2d%Bf_weight_3rd,  &
          sim%m_x1%num_cells,  &
          sim%m_x2%num_cells)
    
     case("nat_per")
     call compute_spl2D_nat_per_weight(&
          sim%field_2d%epotential_init, &
          sim%field_2d%ep_weight_init,  &
          sim%m_x1%num_cells,  &
          sim%m_x2%num_cells)

     call compute_spl2D_nat_per_weight(&
          sim%field_2d%Bf_init_3rd, &
          sim%field_2d%Bf_weight_3rd,  &
          sim%m_x1%num_cells,  &
          sim%m_x2%num_cells)

     case default
       stop
     end select

   case("polar")

  case default
     stop
  end select
   
!       print*, "field_2d%ep_weight_init=", sim%field_2d%ep_weight_init

  eta_min(1)= sim%m_x1%eta_min
  eta_min(2)= sim%m_x2%eta_min
  eta_max(1)= sim%m_x1%eta_max
  eta_max(2)= sim%m_x2%eta_max
  delta_eta(1)= sim%m_x1%delta_eta
  delta_eta(2)= sim%m_x2%delta_eta
  NC(1)= sim%m_x1%num_cells
  NC(2)= sim%m_x2%num_cells
!!$  x(1) = 2.1
!!$  x(2) = 2.3
!!$    
       sim%gyro=>initialize_gyroaverage_2d_plan(gyropoints,mu,sim%m_x1,sim%m_x2,sim%field_2d,geometry,boundary)

! do i=1, sim%m_x1%nodes
!    print*, "i=",i, sim%gyro%ep_weight_gyro(i,:)
! end do

 

  print*,5
  x_pos(1)=x1
  x_pos(2)=x2
  x_pos(3)=x3
  sim%full2d2v=>initialize_particle_orbit_cartesian_2d2v(iter_number, &
       dt_gc/real(num_time,8), &
       N_par, &
       x_pos, &
       mu,geometry )
 
  sim%gyro2d2v => initialize_gyro_orbit_cartesian_2d2v(iter_number, &
       dt_gc, &
!       N_par, &
       x_pos)
 
  do i=1,3
  x4(i)=sim%full2d2v%init_pos_full(1,i)
  v4(i)=sim%full2d2v%init_pos_full(1,i+3)
  end do

!if(boundary=="double_per") then
!  open(12,file="/Users/zhang/program/debug/electrostatic_exp/run/per_per_elef.txt", status="replace")
!else 
!  open(11,file="/Users/zhang/program/debug/electrostatic_exp/run/nat_per_elef.txt", status="replace")
!end if

if(boundary=="double_per") then
 open(12,file="/PARA/blsc950/electrostatic_exp/run/per_per_elef.txt", status="replace")
else 
 open(11,file="/PARA/blsc950/electrostatic_exp/run/nat_per_elef.txt", status="replace")
end if

  call borisrun(x4,v4, &
       sim%field_2d, &
       sim%m_x1, &
       sim%m_x2, &
       sim%iter_number, &
       sim%dt_gc/real(sim%num_time,8), &
       geometry,boundary)
  print*,6

  do i=1,3
  x4(i)=sim%full2d2v%init_pos_full(1,i)
  v4(i)=sim%full2d2v%init_pos_full(1,i+3)
  end do
 
  call borisrun_anal(x4,v4,sim%amplitude,sim%wavenumber_one,&
  sim%wavenumber_two, sim%contr, &
  sim%dt_gc/real(sim%num_time,8), &
  sim%iter_number, sim%geometry)
  
  call fulrkorbit( &
       sim%full2d2v%init_pos_full(1,1:3), &
       sim%full2d2v%init_pos_full(1,4:6), &
       sim%field_2d, &
       sim%m_x1, &
       sim%m_x2, &
       sim%mu,   &
       sim%iter_number, &
       sim%dt_gc/real(sim%num_time,8), &
       geometry,boundary)
 
 if(boundary=="double_per")  then
  close(12)
else
  close(11)
end if
  
   
  print*, 7
 
  if(orbitorder==2) then 
  select case(geometry)
  case ("cartesian")
    select case(boundary)
      case("double_per")
      call compute_spl2d_driftsquare_weight_per_per(sim%field_2d%driftsquare_weight,sim%field_2d%Bf_init_3rd,  &
                                                sim%gyro%ep_weight_gyro,eta_min,eta_max,delta_eta,NC)
      case("nat_per")
      call compute_spl2d_driftsquare_weight_nat_per(sim%field_2d%driftsquare_weight,sim%field_2d%Bf_init_3rd,  &
                                              sim%gyro%ep_weight_gyro,eta_min,eta_max,delta_eta,NC,geometry)
      case default
        stop
    end select
  case ("polar")
      call compute_spl2d_driftsquare_weight_nat_per(sim%field_2d%driftsquare_weight,sim%field_2d%Bf_init_3rd,  &
                                               sim%gyro%ep_weight_gyro,eta_min,eta_max,delta_eta,NC,geometry)
  case default
    stop
  end select
  
  end if
  
print*,9 
  call gyrorkorbit( &
       sim%gyro2d2v%init_pos_gyro(1,1:2), &
       sim%field_2d%Bf_weight_3rd, &
       sim%gyro, &
       sim%field_2d%driftsquare_weight, &
       sim%m_x1, &
       sim%m_x2, &
       sim%mu,   &
       sim%iter_number/sim%num_time, &
       sim%dt_gc, &
       geometry,boundary,orbitorder)
    print*, 8

!  print*,sim%gyro%ep_weight_gyro(:,:)


  
  
end program test_large_epotential


