module paradata_type
#include "work_precision.h"
use piclayout, only: pic_field_2d_base, &
                     parameters_set_2d, &
                     ful2d_node, &
                     ful2dsend_node, &
                     gy2d_node,  &
                     gy2dsend_node, & 
                     parameters_array_2d
use m_mpilayout, only: t_layout_2d
use orbit_data_base, only: tp_ful2d_node,tp_gy2d_node
use utilities_module, only: muarray_euler_maclaurin_choice
implicit none

public :: pic_para_2d_base, &
          pic_para_total2d_base

  type pic_para_2d_base
     type(t_layout_2d),        pointer  :: layout2d
     type(pic_field_2d_base),  pointer  :: field2d
     type(parameters_set_2d),  pointer  :: para2d
!     type(parameters_array_2d), pointer :: para2darray
  end type pic_para_2d_base


  type,extends(pic_para_2d_base)  :: pic_para_total2d_base
     type(ful2d_node),   pointer  :: ful2d_head ! store all the particle on one box
     type(tp_ful2d_node), pointer :: tp_ful2d_head
     type(ful2dsend_node), dimension(:), pointer :: ful2dsend_head  ! storethe migration partilces
     type(gy2d_node),     pointer :: gy2d_head
     type(tp_gy2d_node),  pointer :: tp_gy2d_head
     type(gy2dsend_node), dimension(:), pointer  :: gy2dsend_head
  end type pic_para_total2d_base


end module paradata_type
