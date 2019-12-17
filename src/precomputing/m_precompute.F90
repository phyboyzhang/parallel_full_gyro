module m_precompute
#include "work_precision.h"
use spline_module, only: compute_D_spl2D_per_per_noblock 
use piclayout, only: root_precompute_data, parameters_array_2d
use paradata_type, only: pic_para_total2d_base
use m_parautilities, only: gather_field_to_rootprocess_per_per
use para_gyroaverage_2d_one, only: store_data_on_rootprocess
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

  end subroutine  precompute_ASPL


 

  subroutine precompute_doublegyroaverage_matrix(rootdata,pic2d,pamearray)
      class(root_precompute_data), pointer, intent(inout) :: rootdata
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      class(parameters_array_2d), pointer, intent(in) :: pamearray
      real8, dimension(:), pointer :: density
      int4 :: rank,size,boxindex(4)
      int4 :: ierr, numdim,i,j,mu_num,num1,num2
      real8, dimension(:,:), pointer :: buf,buf1
      real8, dimension(:,:), allocatable :: buf2
      real8 :: mu

      int4 :: LDA, INFO, LWORK
      int4, dimension(:), allocatable :: ipiv
      real8, dimension(:), allocatable :: work

      numdim=pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2
      rank=pic2d%layout2d%collective%rank
      size=pic2d%layout2d%collective%size
      mu_num=pic2d%para2d%mu_num
      boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
      boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
      boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
      boxindex(4)=pic2d%layout2d%boxes(rank)%j_min
      allocate(density(numdim),stat=ierr)
!      allocate(buf(boxindex(2)-boxindex(1)+1,boxindex(4)-boxindex(3)+1))

!      buf= pic2d%field2d%den-pic2d%field2d%denequ
      call gather_field_to_rootprocess_per_per(density,pic2d%field2d%dengeqtot,rank, &
           size,boxindex,pic2d%para2d%numproc,pic2d%layout2d)
     if(rank==0) then
        allocate(buf1(numdim,numdim),buf2(numdim,numdim))
        buf2=0.0
     end if
      do j=1,mu_num
         mu=pamearray%mu_nodes(j)
         if(rank==0) then     
           buf1=0.0
         end if
if(rank==0) then
print*, "double gyroaverge for the jth mu,#j=",j
end if
         call store_data_on_rootprocess(mu,j,rank,rootdata,pic2d)

         if(rank==0) then
           buf1=matmul(rootdata%ACONTRI,rootdata%ASPL)

           do i=1, numdim
             buf1(i,:)=buf1(i,:)*density(i)   ! ????
           end do
           buf1=matmul(rootdata%ASPL,buf1)
           buf2=buf2+matmul(rootdata%ACONTRI,buf1)*pamearray%mu_weights(j)

         endif
      end do

      if(rank==0) then
        do i=1,numdim
          num1 = modulo(i,pic2d%layout2d%global_sz1)
          num2 = (i-num1)/pic2d%layout2d%global_sz1+1 
          buf2(i,i)=buf2(i,i)+1.0_f64/pamearray%temp_i(num1,num2)+1.0_f64/pamearray%temp_e(num1,num2)
        end do


        LDA= numdim
        LWORK = numdim*64
        allocate(IPIV(NUMDIM))
        allocate(work(LWORK))
        call dgetrf(numdim,numdim,buf2,LDA,IPIV,INFO)

        call dgetri(numdim,buf2,LDA,IPIV,WORK,LWORK,INFO)
   
        rootdata%prematrix=buf2
    
      
        deallocate(buf1,buf2)
      endif       

     end subroutine precompute_doublegyroaverage_matrix

end module m_precompute
