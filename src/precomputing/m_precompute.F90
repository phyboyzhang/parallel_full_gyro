module m_precompute
#include "work_precision.h"

use spline_module, only: compute_D_spl2D_per_per_noblock 
use piclayout, only: root_precompute_data
use paradata_type, only: pic_para_total2d_base
use m_parautilities, only: gather_field_to_rootprocess_per_per
implicit none

public :: precompute_ASPL, &
          precompute_doublegyroaverage_matrix 
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

  end subroutine


 
  subroutine precompute_doublegyroaverage_matrix(rootdata,pic2d)
      class(root_precompute_data), pointer, intent(inout) :: rootdata
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      real8, dimension(:), pointer :: density
      int4 :: rank,size,boxindex(4)
      int4 :: ierr, numdim,i,j
      real8, dimension(:,:), pointer :: buf
  
      numdim=pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2 
      rank=pic2d%layout2d%collective%rank
      size=pic2d%layout2d%collective%size
      boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
      boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
      boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
      boxindex(4)=pic2d%layout2d%boxes(rank)%j_min
      allocate(density(numdim),stat=ierr)
      allocate(buf(boxindex(2)-boxindex(1)+1,boxindex(4)-boxindex(3)+1))

      buf= pic2d%field2d%den-pic2d%field2d%denequ    
      call gather_field_to_rootprocess_per_per(density,buf,rank, &
           size,boxindex,pic2d%para2d%numproc,pic2d%layout2d)
      rootdata%prematrix=matmul(rootdata%ACONTRI,rootdata%ASPL)
      do i=1, numdim
         rootdata%prematrix(i,:)=rootdata%prematrix(i,:)*density(i)
      end do
      rootdata%prematrix=matmul(rootdata%ASPL,rootdata%prematrix)
      rootdata%prematrix=matmul(rootdata%ACONTRI,rootdata%prematrix)
     
      deallocate(buf)
   end subroutine precompute_doublegyroaverage_matrix



end module m_precompute
