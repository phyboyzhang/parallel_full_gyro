module orbit_data_base
#include "work_precision.h"
use piclayout, only: ful2d_node,gy2d_node
implicit none


!!!!! The orbit data types for the integrated simuations

  type :: rk4ful2dnode
     real8 :: f1(1:6),f2(1:6),f3(1:6),f4(1:6)
     real8 :: vec(1:6)
!     real8 :: vec0(1:6)
     real8 :: coords(1:4)
     int4 :: at    !! at=0, the point locates at the current box; at=1,out the current box 
     type(rk4ful2dnode), pointer :: next=>null()
  end type

  type rk4ful2dout_node   !fuloutnode
     int4  :: numpoint
     real8 :: coords(2)
     type(rk4ful2dout_node), pointer :: next=>null()
  end type

  type fuloutnode
    type(rk4ful2dout_node), pointer :: ptr
  end type

  type pointin_node
    real8 :: x(2)
    real8 :: numpoint
    type(pointin_node), pointer :: next =>null()
  end type

  type pointin_gy_node
    real8 :: x(3)
    real8 :: numpoint
    type(pointin_gy_node), pointer :: next=>null()
  end type

  type :: rk4gy2dnode
     real8 :: f1(1:2),f2(1:2),f3(1:2),f4(1:2)
     real8 :: vec(1:2)
     real8 :: coords(1:3)
     int4 :: at
     type(rk4gy2dnode), pointer :: next=>null()
  end type

  type  rk4gy2dout_node  !  gyoutnode
     int4  :: numpoint
     real8 :: coords(3)
     type(rk4gy2dout_node), pointer :: next=>null()
  end type

  type  gyoutnode
    type(rk4gy2dout_node),pointer :: ptr
  end type

!!!! The orbit data structures for the test particles

  type :: tp_ful2d_node
     real8 :: coords(1:4)
     int4  :: tp
     type(tp_ful2d_node), pointer :: next=>null()
  end type

  type tp_ful2drank_node    ! used to store particles running out of the current box
     real8 :: coords(1:4)
     int4  :: prank
     int4  :: tp
     type(tp_ful2drank_node), pointer :: next=>null()
  end type tp_ful2drank_node

  type tp_ful2dsend_node
     type(tp_ful2drank_node), pointer :: ptr
  end type tp_ful2dsend_node

  type tp_gy2drank_node
     real8 :: coords(1:3)
     int4  :: prank
     int4  :: tp
     type(tp_gy2drank_node), pointer :: next=>null()
  end type tp_gy2drank_node

  type tp_gy2dsend_node
     type(tp_gy2drank_node),pointer :: ptr
  end type tp_gy2dsend_node

  type :: tp_gy2d_node
     real8 :: coords(1:3)
     int4 :: tp
     type(tp_gy2d_node), pointer :: next=>null()
  end type 

  type :: tp_gy2dallmu_node
     type(tp_gy2d_node),pointer :: ptr
  end type 
 
  type :: tp_rk4ful2dnode
     real8 :: f1(6),f2(6),f3(6),f4(6)
     real8 :: vec(6)
     real8 :: coords(4)
     int4 :: at    !! at=0, the point locates at the current box; at=1,out the current box 
     int4 :: tp   !! The index number used to denote the tp-th particle
     type(tp_rk4ful2dnode), pointer :: next=>null()
  end type 

  type :: tp_rk4gy2dnode
     real8 :: f1(2),f2(2),f3(2),f4(2)
     real8 :: vec(2)
     real8 :: coords(3)
     int4 :: at
     int4 :: tp
     type(tp_rk4gy2dnode), pointer :: next=>null()
  end type

  public :: pointin_node, pointin_gy_node,fuloutnode,rk4ful2dnode,gyoutnode,rk4gy2dnode,&
            tp_rk4gy2dnode, tp_ful2dsend_node, tp_gy2dsend_node, &
            tp_gy2dallmu_node
        


end module orbit_data_base
