module cartesian_mesh
#include "work_precision.h"
use constants, only: pi_
use abstract_mesh_base, only: mesh_1d_base, &
                              mesh_2d_base
  
implicit none

public :: &
     cartesian_mesh_1d, &
     f_new_cartesian_mesh_1d, &
     f_new_cartesian_mesh_2d, &
     eta1_node_1d, &
     eta1_node_2d, &
     eta2_node_2d, &
     init_para_cartesian_mesh_1d
  
  
  type, extends(mesh_1d_base) :: cartesian_mesh_1d
   contains
     procedure, pass(mesh) :: eta1_node => eta1_node_1d
!     procedure, pass(mesh) :: eta1_cell => eta1_cell_1d
!     procedure, pass(mesh) :: delete => cartesian_mesh_1d_free
  end type cartesian_mesh_1d

  type, extends(mesh_2d_base) :: cartesian_mesh_2d

  end type cartesian_mesh_2d
  
  contains
!======>initializating the type

  function f_new_cartesian_mesh_1d( &
       nodes, &
       eta_min, &
       eta_max, &
       mesh_boundary) result(m)

    type(cartesian_mesh_1d), pointer :: m
    int4, intent(in) :: nodes
    real8,intent(in) :: eta_min
    real8, intent(in) :: eta_max
    character(len=*), intent(in) :: mesh_boundary
    int4 :: ierr
    allocate(m, stat=ierr)
    call cartesian_mesh_1d_init(m,nodes,eta_min,eta_max,mesh_boundary)
  end function f_new_cartesian_mesh_1d

  function init_para_cartesian_mesh_1d(nodes,eta_min,eta_max,delta_eta) result(mesh1d)
     type(cartesian_mesh_1d), pointer :: mesh1d
     int4,intent(in) :: nodes
     real8, intent(in) :: eta_min,eta_max,delta_eta
     int4 :: ierr
      
     allocate(mesh1d)
     mesh1d%nodes=nodes
     mesh1d%eta_min=eta_min
     mesh1d%eta_max=eta_max
     mesh1d%delta_eta=delta_eta
  end function

  subroutine cartesian_mesh_1d_init(m,nodes,eta_min,eta_max,mesh_boundary)
    class(cartesian_mesh_1d), intent(inout) :: m
    int4, intent(in) :: nodes
    real8, optional, intent(in) :: eta_min
    real8, optional, intent(in) :: eta_max
    character(len=*), intent(in) :: mesh_boundary
    
 select case (mesh_boundary)
  case ("periodic")
    m%num_cells = nodes
    m%nodes   = nodes
    m%eta_min = eta_min
    m%eta_max = eta_max 
    m%delta_eta = (m%eta_max-m%eta_min)/real(nodes, 8)

  case ("natural")
    m%num_cells = nodes-1
    m%nodes   = nodes
    m%eta_min = eta_min
    m%eta_max = eta_max 
    m%delta_eta = (m%eta_max-m%eta_min)/real(nodes-1, 8)
    
  case default
    print*, "The input mesh_boundary is not correct, mesh_boudnary=", mesh_boundary
    stop
 end select
    
    if(m%eta_max .le. m%eta_min) then
       print*, "error, cartesian_mesh_1d_init():", &
            "problem to construct the mesh 1d"
       print*, "because eta_max <= eta_min"
       stop
    end if
    
  end subroutine cartesian_mesh_1d_init


  function eta1_node_1d(mesh, i) result(res)
    class(cartesian_mesh_1d), intent(in) :: mesh
    int4, intent(in) :: i
    real8            :: res
    real8            :: eta_min
    real8            :: delta_eta

    eta_min   = mesh%eta_min
    delta_eta = mesh%delta_eta
    res       = eta_min + real(i-1,8)*delta_eta
  end function eta1_node_1d


  subroutine  cartesian_mesh_1d_free(mesh)
    class(cartesian_mesh_1d), intent(inout) :: mesh
    
  end subroutine cartesian_mesh_1d_free

  
  function f_new_cartesian_mesh_2d( &
       num_cells1, &
       num_cells2, &
       eta1_min, &
       eta1_max, &
       eta2_min, &
       eta2_max ) result(m)
    type(cartesian_mesh_2d), pointer :: m
    int4, intent(in) :: num_cells1
    int4, intent(in) :: num_cells2
    real8, optional, intent(in) :: eta1_min
    real8, optional, intent(in) :: eta1_max
    real8, optional, intent(in) :: eta2_min
    real8, optional, intent(in) :: eta2_max
    int4 :: ierr

    allocate(m,stat=ierr)

    call cartesian_mesh_2d_init(&
         m,&
         num_cells1, &
         num_cells2, &
         eta1_max, &
         eta2_min, &
         eta2_max )

  end function f_new_cartesian_mesh_2d

  subroutine cartesian_mesh_2d_init(&
       m, &
       num_cells1, &
       num_cells2, &
       eta1_min, &
       eta1_max, &
       eta2_min, &
       eta2_max)

    class(cartesian_mesh_2d), intent(inout) :: m
    int4, intent(in) :: num_cells1
    int4, intent(in) :: num_cells2
    real8, optional, intent(in) :: eta1_min
    real8, optional, intent(in) :: eta1_max
    real8, optional, intent(in) :: eta2_min
    real8, optional, intent(in) :: eta2_max

    m%num_cells1 = num_cells1
    m%num_cells2 = num_cells2
    m%eta1_min = eta1_min
    m%eta1_max = eta1_max
    m%eta2_min = eta2_min
    m%eta2_max = eta2_max
    m%delta_eta1   = (m%eta1_max - m%eta1_min)/real(num_cells1,8)
    m%delta_eta2   = (m%eta2_max - m%eta2_min)/real(num_cells2,8)

    if ( m%eta1_max <= m%eta1_min) then
       print*,'cartesian_mesh_2d_init():','Problem to construct the mesh 2d because eta1_max <= eta1_min'
    end if
    if ( m%eta2_max <= m%eta2_min) then
       print*,'cartesian_mesh_2d_init():','Problem to construct the mesh 2d because eta2_max <= eta2_min'
    end if

    end subroutine cartesian_mesh_2d_init
    
  function eta1_node_2d(mesh, i, j) result(res)
    class(cartesian_mesh_2d), intent(in) :: mesh
    int4, intent(in) :: i
    int4, intent(in) :: j
    real8            :: res
    real8            :: eta1_min
    real8            :: delta_eta1


    eta1_min   = mesh%eta1_min
    delta_eta1 = mesh%delta_eta1
    res        = eta1_min + real(i-1,8)*delta_eta1
  end function eta1_node_2d

  function eta2_node_2d(mesh, i, j) result(res)
    class(cartesian_mesh_2d), intent(in) :: mesh
    int4, intent(in) :: i
    int4, intent(in) :: j
    real8            :: res
    real8            :: eta2_min
    real8            :: delta_eta2

    eta2_min   = mesh%eta2_min
    delta_eta2 = mesh%delta_eta2
    res        = eta2_min + real(j-1,8)*delta_eta2
  end function eta2_node_2d

    
end module cartesian_mesh
