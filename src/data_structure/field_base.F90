module field_base
#include "work_precision.h"

  implicit none

  type, abstract :: field_plan_base
    contains
    procedure(initial_field_plan), deferred, pass(field_mesh) :: init_field
 end type field_plan_base
 
 abstract interface
    subroutine initial_field_plan
      import field_plan_base
      class(field_plan_base), target :: field_2d
      int4, dimension(:), intent(in) :: Nc
      int4,intent(in) :: flag
    end subroutine initial_field_plan
 end interface

end module field_base
