module gyroaverage_2d_base
#include "work_precision.h"
!use cartesian_mesh, only: cartesian_mesh_1d
!use m_work_precision
implicit none

  public :: gyropoint_node, &
            gyropoint, &
            gyroweight            

!  type, abstract :: s_gyroaverage_2d_base
!     int4  :: N_points    !! number of points on the circle
!     real8, dimension(:,:), pointer :: Bfield !! the magnetic field on the mesh
!     real8, dimension(:,:), pointer :: efield  !! the electrostatic field on the mesh & after the gyroaverage
!     real8, dimension(:), pointer :: weight    !! the weight of the spline function on the mesh     
!   contains
!     procedure(compute_gyroaverage_2d), deferred, pass(gyroaverage) :: gyroaverage_2d_point
!  end type s_gyroaverage_2d_base
!
!  abstract interface
!     subroutine compute_gyroaverage_2d(m_x1,m_x2,gyro,points,mu,x,poten)
!      use cartesian_mesh, only: cartesian_mesh_1d
!       import s_gyroaverage_2d_base
!       implicit none
!       class(s_gyroaverage_2d_base), target :: gyro
!       class(cartesian_mesh_1d), intent(in) :: m_x1, m_x2
!       real8, intent(in) :: mu
!       real8, dimension(:,:), intent(in) :: points
!       real8, dimension(2),intent(in) :: x
!       real8, intent(inout) :: poten
!
!     end subroutine compute_gyroaverage_2d
!  end interface


  type gyroweight
    int4 :: gridrank, ptrank  ! gridrank: the rank of box where the grid point locates
    real8 :: xpoint(2), xgrid(2)  ! the coordinates of the point and the grid
    real8 :: weight(-1:1,-1:1), val(-1:1,-1:1) !weight is the weight for the spline. 
 !   val is the value of the spline at those points
    type(gyropoint), pointer :: next=>null()
  end type gyroweight

  type gyropoint_node
    type(gyropoint), pointer :: ptr
  end type

  type gyropoint
    int4 :: gridrank  ! gridrank: the rank of box where the grid point locates
    real8 :: xpoint(2)  ! the coordinates of the point 
    int4  :: gridind(2) ! the index of the mesh grid
    real8 :: weight(-1:2,-1:2), val(-1:2,-1:2) !weight is the weight for the
 !   spline. val is the value of the spline at those points
    type(gyropoint), pointer :: next=>null()
  end type gyropoint



!  type,extends(gyroaverage_2d_plan):: para_gyroaverage_2d_plan
!    real8 :: mu,rho
!    real8, dimension(:,:), pointer :: boxextend
!    int4  :: nodes(2)   ! the number of node of boxextend
!    int4  :: lvortindex(4,2)  ! the index of the vortices of the current box within the boxextend
!    int4  :: gvortcoord(4,2)  ! the global coords of the 4 vortices of the boxextend
!  end type para_gyroaverage_2d_plan




end module gyroaverage_2d_base
