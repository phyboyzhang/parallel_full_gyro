module orbit_2d_2v
#include "work_precision.h"

  implicit none

  public :: &
  particle_orbit_cartesian_2d2v, &
  initialize_particle_orbit_cartesian_2d2v, &
  gyro_orbit_cartesian_2d2v, &
  initialize_gyro_orbit_cartesian_2d2v

  type particle_2d2v_base   ! initial particle 2d+2v phase space
     real8, dimension(:,:), pointer :: init_pos_full  !!! first the number of particles,
                                                    !!! second x_1,x_2,v_1,v_2
     real8, dimension(:,:,:), pointer :: orbit_full   !!! first the number of particles,
                                                    !!! second the iter step number
                                                   !!! third  x_1,x_2,v_1,v_2
  end type particle_2d2v_base
  
  type gyro_2d2v_base
     real8, dimension(:,:), pointer :: init_pos_gyro  !!! first the number of particles,
                                                         !!! second x_1,x_2,\mu
     real8, dimension(:,:,:), pointer :: orbit_gyro     
  end type gyro_2d2v_base
 

  type, extends(particle_2d2v_base) ::  particle_orbit_cartesian_2d2v
     int4 :: iter_number    !! the iteration number
     real8 :: delta_t       !! time_step
 !  contains
 !     procedure :: initialize_particle_orbit_cartesian_2d2v
 !     procedure, pass(orbit_2d2v) :: compute_particle_orbit    
   end type particle_orbit_cartesian_2d2v

   type, extends(gyro_2d2v_base) :: gyro_orbit_cartesian_2d2v
     int4 :: iter_number
     real8 :: delta_t
 !  contains
 !    procedure :: initialize_gyro_orbit_cartesian_2d2v
 !    procedure, pass(gyro_orbit_2d2v_base) :: compute_gyro_orbit
  end type gyro_orbit_cartesian_2d2v
  
contains
  
  function initialize_particle_orbit_cartesian_2d2v(iter_number, &
       delta_t, &
       N_par, &
       x, &
       mu, geometry ) result(full2d2v)
    type(particle_orbit_cartesian_2d2v), pointer  :: full2d2v
    real8, intent(in) :: delta_t, mu, x(3)
    int4,  intent(in) :: N_par
    int4,  intent(in) :: iter_number
    character(len=*),intent(in) :: geometry
    real8 :: rho

    allocate(full2d2v)
    rho=sqrt(mu)
    full2d2v%iter_number=iter_number
    full2d2v%delta_t=delta_t

!    select case (geometry)
!    case ("cartesian")
      if(N_par==1) then
       allocate(full2d2v%init_pos_full(1,6))
       full2d2v%init_pos_full(1,1)=x(1)
       full2d2v%init_pos_full(1,2)=x(2)+rho
       full2d2v%init_pos_full(1,3)=x(3)
       full2d2v%init_pos_full(1,4)=sqrt(mu)
       full2d2v%init_pos_full(1,5)=0.0_f64
       full2d2v%init_pos_full(1,6)=0.0_f64
      end if
 !    case ("polar")

 !    case default
 !      stop
 !    end select 
  
    return
  end function initialize_particle_orbit_cartesian_2d2v

  
  
   function initialize_gyro_orbit_cartesian_2d2v(iter_number, &
       delta_t, &
       x) result(gyro2d2v)
     type(gyro_orbit_cartesian_2d2v),pointer :: gyro2d2v
     int4,  intent(in) :: iter_number
     real8, intent(in) :: delta_t,x(3)
!    int4,  intent(in) :: N_par
!    real8 :: rho
     allocate(gyro2d2v)
    gyro2d2v%iter_number=iter_number
    gyro2d2v%delta_t=delta_t

!    if(N_par==1) then
       allocate(gyro2d2v%init_pos_gyro(1,3))
       gyro2d2v%init_pos_gyro(1,1:3)=x(1:3)
!    end if
  end function initialize_gyro_orbit_cartesian_2d2v
 




end module
