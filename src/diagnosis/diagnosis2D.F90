module diagnosis2D
#include "work_precision.h"
use m_parautilities, only: gather_field_to_rootprocess_per_per
use paradata_type, only: pic_para_total2d_base
use m_mpilayout,only: get_layout_2d_box_index
implicit none
include "mpif.h"
   
    public :: compare_density_to_initnumber_gy

contains

   subroutine compare_density_to_initnumber_gy(boxfield,pic2d)
     class(pic_para_total2d_base), pointer, intent(in) :: pic2d
     real8, dimension(:,:), pointer,intent(in) :: boxfield
     int4 :: rank,size,boxindex(4),global_sz(2),globind
     real8 :: gxmin(2),delta(2),x(2)
     real8, dimension(:), pointer :: buff
     int4 :: i,j  
     real8 :: sum

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
    global_sz=(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/)
    gxmin=pic2d%para2d%gxmin
    delta(1)=pic2d%para2d%m_x1%delta_eta
    delta(2)=pic2d%para2d%m_x2%delta_eta

    allocate(buff(global_sz(1)*global_sz(2)))
    call gather_field_to_rootprocess_per_per(buff,boxfield,rank,size,boxindex, &
         pic2d%para2d%numproc,pic2d%layout2d)    
    sum=0.0
    if(rank==0) then
     do i=1,global_sz(1)*global_sz(2)
        sum=sum+buff(i)
      end do
    end if

    if(rank==0) then
      print*, "The total number of particles after partition is", sum
    end if
    deallocate(buff)
  end subroutine



end module diagnosis2D
