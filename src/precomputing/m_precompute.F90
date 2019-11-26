module m_precompute
#include "work_precision.h"
use spline_module, only: compute_D_spl2D_per_per_noblock 
use piclayout, only: root_precompute_data
use paradata_type, only: pic_para_total2d_base
use m_parautilities, only: gather_field_to_rootprocess_per_per
implicit none

public :: precompute_ASPL
          
contains

  subroutine precompute_ASPL(rank,global_sz,ASPL)
    int4, dimension(:), intent(in) :: global_sz
    real8, dimension(:,:), intent(inout) :: ASPL
    int4, intent(in) :: rank  

    if(rank==0) then
    call compute_D_spl2D_per_per_noblock( &
         global_sz(1), &
         global_sz(2), &
         ASPL)
    end if

  end subroutine  precompute_ASPL


 


end module m_precompute
