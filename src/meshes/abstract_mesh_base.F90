module abstract_mesh_base
#include "work_precision.h"
  implicit none
  
  public :: &
       mesh_1d_base, &
       mesh_2d_base
  

  private


  type, abstract :: mesh_1d_base
     int4 :: num_cells
     int4 :: nodes
     real8 :: eta_min
     real8 :: eta_max
     real8 :: delta_eta
     contains
     procedure(get_geometry_1d), deferred, pass(mesh) :: eta1_node
  !   procedure(get_geometry_1d), deferred, pass(mesh) :: eta1_cell
  !   procedure(display_mesh_1d), deferred, pass :: display
  !   procedure(delete_mesh_1d),  deferred, pass :: delete
  end type

  abstract interface
     function get_geometry_1d(mesh,i) result(res)
       import mesh_1d_base
       class(mesh_1d_base), intent(in) :: mesh
       int4, intent(in) :: i
       real8 :: res 
     end function get_geometry_1d    
  end interface

  type, abstract :: mesh_2d_base
     int4 :: num_cells1
     int4 :: num_cells2
     real8 :: eta1_min
     real8 :: eta1_max
     real8 :: eta2_min
     real8 :: eta2_max
     real8 :: delta_eta1
     real8 :: delta_eta2

  end type mesh_2d_base
  
  
  abstract interface
     subroutine delete_mesh_1d(mesh)
       import mesh_1d_base
       class(mesh_1d_base), intent(inout) :: mesh
     end subroutine delete_mesh_1d
  end interface

end module
