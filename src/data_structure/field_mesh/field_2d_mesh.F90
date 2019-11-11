module field_2d_mesh
#include "work_precision.h"
  use cartesian_mesh, only: cartesian_mesh_1d
  
!  USE spline_module, only: &
!       compute_spl2D_double_per_weight, &
!       compute_spl2D_nat_per_weight
  
  use utilities_module,only: gp_error

  implicit none

  public :: field_2d_plan,  &
       pic_field_2d_base, &
       initialize_field_2d_plan, &
       initialize_pic_field_2d_base
       
  
 !           compute_spl2D_field_weight_2d_plan
       

  type :: field_2d_plan
    real8, dimension(:,:), pointer :: Bf_init_1st !! the magnetic field on the mesh
    real8, dimension(:,:), pointer :: Bf_init_2nd !! the magnetic field on the mesh     
    real8, dimension(:,:), pointer :: Bf_init_3rd !! the magnetic field on the mesh
    real8, dimension(:,:), pointer :: epotential_init  !! the electrostatic potential on the mesh before the gyroaverage
    real8, dimension(:,:), pointer :: Bf_weight_3rd
    real8, dimension(:,:), pointer :: ep_weight_init
    real8, dimension(:,:), pointer :: driftsquare_weight
!  contains
!    procedure,  pass(field_mesh) :: init_plan => initialize_field_cartesian_2d_plan
 end type field_2d_plan

 type :: pic_field_2d_base
    real8, dimension(:,:), pointer :: Bf_init_x !! the magnetic field on the mesh
    real8, dimension(:,:), pointer :: Bf_init_y !! the magnetic field on the mesh     
    real8, dimension(:,:), pointer :: Bf_init_z !! the magnetic field on the mesh
    real8, dimension(:,:), pointer :: Bf_weight
    real8, dimension(:,:), pointer :: Bf_weight_x    
    real8, dimension(:,:), pointer :: Bf_weight_y      
    real8, dimension(:,:), pointer :: Bf_weight_z
    real8, dimension(:,:), pointer :: ep  !! the electrostatic potential on the mesh before the gyroaverage    
    real8, dimension(:,:), pointer :: ep_weight    
    real8, dimension(:,:), pointer :: density
 end type pic_field_2d_base


contains

    function initialize_field_2d_plan( &
         m_x1,            &
         m_x2,            &
         geometry) result(field_mesh) ! result(field_2d)
       type(field_2d_plan),pointer  :: field_mesh  
       class(cartesian_mesh_1d), pointer :: m_x1, m_x2
       character(len=*), intent(in) :: geometry
      real8 :: x(2)
      int4 :: ierr
      int4 :: i, j
      allocate(field_mesh)
      select case(geometry)
      case ("cartesian")
!         select case (boundary)
!         case ("double_per") 
         allocate(field_mesh%Bf_init_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
         call gp_error(ierr,"Bf_init_3rd")
         allocate(field_mesh%Bf_weight_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
         call gp_error(ierr,"Bf_weight_3rd")         
         allocate(field_mesh%epotential_init(m_x1%nodes, m_x2%nodes),stat=ierr)
         call gp_error(ierr,"epotential_init")
         allocate(field_mesh%ep_weight_init(m_x1%nodes, m_x2%nodes), stat=ierr)
         call gp_error(ierr,"ep_weight_init")
         
!         case ("nat_per")
!         allocate(field_mesh%Bf_init_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
!         call gp_error(ierr,"Bf_init_3rd")
!         allocate(field_mesh%Bf_weight_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
!         call gp_error(ierr,"Bf_weight_3rd")         
!         allocate(field_mesh%epotential_init(m_x1%nodes, m_x2%nodes),stat=ierr)
!         call gp_error(ierr,"epotential_init")
!         allocate(field_mesh%ep_weight_init(m_x1%nodes, m_x2%nodes), stat=ierr)
!         call gp_error(ierr,"ep_weight_init")
! 
!         case default
!           stop
!         end select

         allocate(field_mesh%driftsquare_weight(m_x1%nodes,m_x2%nodes),stat=ierr)
         call gp_error(ierr,"driftsquare_weight")

      case ("polar")
      
         allocate(field_mesh%Bf_init_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
         call gp_error(ierr,"Bf_init_3rd")
         allocate(field_mesh%Bf_weight_3rd(m_x1%nodes,m_x2%nodes), stat=ierr)
         call gp_error(ierr,"Bf_weight_3rd")    
         allocate(field_mesh%epotential_init(m_x1%nodes, m_x2%nodes), stat=ierr)
         call gp_error(ierr,"epotential_init")
         allocate(field_mesh%ep_weight_init(m_x1%nodes, m_x2%nodes), stat=ierr)
         call gp_error(ierr,"ep_weight_init")
         allocate(field_mesh%driftsquare_weight(m_x1%nodes,m_x2%nodes),stat=ierr)
         call gp_error(ierr,"driftsquare_weight")
         
   case default
      print*, "enter the geometry  condition, geometry=", geometry
      stop
   end select


 end function initialize_field_2d_plan



 function initialize_pic_field_2d_base(m_x1, m_x2) result(field2d)
   class(cartesian_mesh_1d), pointer,intent(in) :: m_x1,m_x2
   type(pic_field_2d_base), pointer  :: field2d
   int4 :: ierr
   allocate(field2d)
   
   allocate(field2d%Bf_init_z(m_x1%nodes,m_x2%nodes), stat=ierr)
   call gp_error(ierr,"Bf_init_z")
   allocate(field2d%Bf_weight(m_x1%nodes,m_x2%nodes), stat=ierr)
   call gp_error(ierr,"Bf_weight")         
   allocate(field2d%ep(m_x1%nodes, m_x2%nodes),stat=ierr)
   call gp_error(ierr,"ep")
   allocate(field2d%ep_weight(m_x1%nodes, m_x2%nodes), stat=ierr)
   call gp_error(ierr,"ep_weight")  

 end function initialize_pic_field_2d_base

   
!!$ subroutine compute_spl2D_field_weight_2d_plan(&
!!$      field_mesh, &
!!$      m_x1, &
!!$      m_x2, &
!!$      geometry)
!!$
!!$   class(field_2d_plan), pointer :: field_mesh
!!$   class(cartesian_mesh_1d), pointer :: m_x1,m_x2
!!$   character(len=*), intent(in) :: geometry
!!$
!!$   select case(geometry)
!!$   case ("cartesian")
!!$
!!$         call compute_spl2D_double_per_weight(&
!!$              field_mesh%epotential_init, &
!!$              field_mesh%ep_weight_init,  &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         
!!$         call compute_spl2D_doulbe_per_weight(&
!!$              field_mesh%Bf_init_3rd, &
!!$              field_mesh%Bf_weight_3rd, &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$
!!$
!!$   case ("polar")
!!$          call compute_spl2D_nat_per_weight(&
!!$              field_mesh%epotential_init, &
!!$              field_mesh%ep_weight_init,  &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         call compute_spl2D_nat_per_weight(&
!!$              field_mesh%Bf_init_3rd, &
!!$              field_mesh%Bf_weight_3rd, &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         
!!$   case default
!!$      print*, "The input geometry is not right, geometry=", geometry
!!$      stop
!!$   end select
!!$   
!!$ end subroutine compute_spl2D_field_weight_2d_plan


!!$ subroutine compute_spl2D_field_weight_2d_plan(&
!!$      epotential_init, &
!!$      ep_weight_init,  &
!!$      Bf_init_3rd,     &
!!$      Bf_weight_3rd,   &
!!$      m_x1, &
!!$      m_x2, &
!!$      geometry)
!!$
!!$   real8, dimension(:,:), pointer,intent(in) :: epotential_init, Bf_init_3rd
!!$   real8, dimension(:,:), pointer,intent(inout) ::  ep_weight_init, Bf_weight_3rd
!!$   class(cartesian_mesh_1d), pointer :: m_x1,m_x2
!!$   character(len=*), intent(in) :: geometry
!!$
!!$   select case(geometry)
!!$   case ("cartesian")
!!$
!!$         call compute_spl2D_double_per_weight(&
!!$              epotential_init, &
!!$              ep_weight_init,  &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         
!!$         call compute_spl2D_doulbe_per_weight(&
!!$              Bf_init_3rd, &
!!$              Bf_weight_3rd, &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$
!!$
!!$   case ("polar")
!!$          call compute_spl2D_nat_per_weight(&
!!$              epotential_init, &
!!$              ep_weight_init,  &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         call compute_spl2D_nat_per_weight(&
!!$              Bf_init_3rd, &
!!$              Bf_weight_3rd, &
!!$              m_x1%num_cells, &
!!$              m_x2%num_cells)
!!$         
!!$   case default
!!$      print*, "The input geometry is not right, geometry=", geometry
!!$      stop
!!$   end select
!!$
!!$end subroutine compute_spl2D_field_weight_2d_plan
 

end module field_2d_mesh

