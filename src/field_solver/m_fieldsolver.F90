module m_fieldsolver
#include "work_precision.h"

use piclayout, only :   root_precompute_data, parameters_array_2d
use m_para_spline, only: para_compute_spl2D_weight
use paradata_type, only: pic_para_total2d_base
use m_mpilayout, only: get_layout_2d_box_index
use m_parautilities, only: gather_field_to_rootprocess_per_per, &
                           mpi2d_alltoallv_box_per_per, &
                           scatter_field_from_rootprocess_per_per 
use para_gyroaverage_2d_one, only: para_compute_gyroaverage_mesh_field, &
                                   para_compute_gyroaverage_field_on_mesh
implicit none

public :: solve_weight_of_field_among_processes, &
          solve_field_quasi_neutral, &
          solve_gyfieldweight_from_fulfield, &
          compute_gyrodensity_perturbation, &
          compute_equdensity_for_ful, &
          compute_equdensity_for_gy, &
          solve_field_ful
contains

  !!! This subroutine is used to solve the weight of field. First, gather
  !field data to the root process. Then, obtain the weight coefficient. Last,
  !scatter all weight to each process. 
    subroutine solve_weight_of_field_among_processes(boxfield,rootmatrix,rootdata,pic2d,boxweight, &
               rw,re,rn,rs,rsw,rse,rnw,rne)
      real8,dimension(:,:), pointer, intent(in) :: boxfield
      real8,dimension(:,:), pointer, intent(in) :: rootmatrix
      class(root_precompute_data), pointer, intent(inout) :: rootdata
      class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
      real8, dimension(:,:), pointer, intent(inout) :: boxweight, &
                                                       rw,re,rn,rs,rsw,rse,rnw,rne
      int4 :: rank,size, boxindex(4), row, comm, numproc(2)

      rank=pic2d%layout2d%collective%rank
      size=pic2d%layout2d%collective%size
      row =pic2d%para2d%row
      comm=pic2d%layout2d%collective%comm
      numproc=pic2d%para2d%numproc
      call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)

      call  gather_field_to_rootprocess_per_per(rootdata%field,boxfield,rank,size,boxindex,&
            pic2d%para2d%numproc,pic2d%layout2d)

      call  para_compute_spl2D_weight(rootmatrix,rootdata%field,boxweight,pic2d%para2d%numproc, &
            pic2d%layout2d,pic2d%para2d%boundary)

      call  mpi2d_alltoallv_box_per_per(row,comm,rank,numproc,boxweight,rw,re,rn,rs,  &
            rsw,rse,rnw,rne,boxindex) 

     end subroutine  solve_weight_of_field_among_processes


     subroutine solve_field_quasi_neutral(rank,rootdata,pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       int4, intent(in) :: rank
       real8, dimension(:), pointer :: density
       int4 :: dims,global_sz(2),sizeone,boxindex(4),num1,num2
       real8, dimension(:,:), pointer :: buffer
       int4 :: i,j

       num1=size(pic2d%field2d%dengtot,1)
       num2=size(pic2d%field2d%dengtot,2)
       global_sz(1)=pic2d%layout2d%global_sz1
       global_sz(2)=pic2d%layout2d%global_sz2
       sizeone=pic2d%layout2d%collective%size     
       dims=global_sz(1)*global_sz(2)
       allocate(density(dims))
       allocate(buffer(num1,num2))
       
       do i=1,num1
         do j=1,num2     
            buffer(i,j)=pic2d%field2d%dengtot(i,j)/pic2d%field2d%dengeqtot(i,j)
         end do
       enddo
!!!!!! Here, it may be required to carry out the gyroaverage of
!!!pic2d%field2d%deng
       call get_layout_2d_box_index( pic2d%layout2d, rank, boxindex )       
       call gather_field_to_rootprocess_per_per(density,buffer,rank, &
           sizeone,boxindex,pic2d%para2d%numproc,pic2d%layout2d)
       if(rank==0) then
          rootdata%field=matmul(rootdata%prematrix,density)
       end if 

       call scatter_field_from_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,sizeone, &
            pic2d%para2d%numproc,global_sz,pic2d%layout2d)

       deallocate(density,buffer)
     end subroutine solve_field_quasi_neutral

     subroutine solve_field_ful(rootdata,pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata 
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       int4 :: i,j,num1,num2
 
       num1=size(pic2d%field2d%denf,1)
       num2=size(pic2d%field2d%denf,2)
       do i=1,num1
         do j=1,num2
            pic2d%field2d%denf(i,j)=pic2d%field2d%denf(i,j)*pamearray%temp_e(i,j)/ &
                                    pic2d%field2d%denfeq(i,j) 
         end do
       enddo
       
         
     end subroutine

     subroutine solve_gyfieldweight_from_fulfield(rootdata, pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       real8,dimension(:,:),pointer :: buf,box,re,rs,rw,rn,rne,rse,rsw,rnw 
       int4 :: i,ierr,num1,num2,row 

       num1=size(pic2d%field2d%gep,1)
       num2=size(pic2d%field2d%gep,2)
       row=pic2d%para2d%row
       allocate(buf(num1,num2))       
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2), &
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)
       
       call solve_weight_of_field_among_processes(pic2d%field2d%gep,rootdata%ASPL,rootdata,pic2d, &
       pic2d%field2d%gep_weight, pic2d%field2d%gepwg_w,pic2d%field2d%gepwg_e,pic2d%field2d%gepwg_n, &
       pic2d%field2d%gepwg_s, pic2d%field2d%gepwg_sw,pic2d%field2d%gepwg_se, &
       pic2d%field2d%gepwg_nw,pic2d%field2d%gepwg_ne)

       do i=1,pic2d%para2d%mu_num
         buf=0.0
         box=0.0
         rw=0.0
         re=0.0
         rn=0.0
         rs=0.0
         rsw=0.0
         rse=0.0
         rnw=0.0
         rne=0.0 
         call para_compute_gyroaverage_field_on_mesh(pamearray%mu_nodes(i),i,pic2d,buf, &
              pic2d%field2d%gep_weight,pic2d%field2d%gepwg_w,pic2d%field2d%gepwg_e,pic2d%field2d%gepwg_n, &
              pic2d%field2d%gepwg_s, pic2d%field2d%gepwg_sw,pic2d%field2d%gepwg_se, &
              pic2d%field2d%gepwg_nw,pic2d%field2d%gepwg_ne)

         call solve_weight_of_field_among_processes(buf,rootdata%ASPL,rootdata,pic2d, &
              box,rw,re,rn,rs,rsw,rse,rnw,rne)

         pic2d%field2d%epgy_weight(i,:,:)=box
         pic2d%field2d%epgywg_w(i,:,:)=rw
         pic2d%field2d%epgywg_e(i,:,:)=re
         pic2d%field2d%epgywg_n(i,:,:)=rn 
         pic2d%field2d%epgywg_s(i,:,:)=rs
         pic2d%field2d%epgywg_sw(i,:,:)=rsw 
         pic2d%field2d%epgywg_se(i,:,:)=rse
         pic2d%field2d%epgywg_nw(i,:,:)=rnw
         pic2d%field2d%epgywg_ne(i,:,:)=rne
       
       end do    

     end subroutine solve_gyfieldweight_from_fulfield


     subroutine compute_equdensity_for_gy(pic2d)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       int4 :: i      
!       pic2d%field2d%deng=pic2d%field2d%deng-pic2d%field2d%dengeq
       pic2d%field2d%dengeqtot(:,:)=0.0
       pic2d%field2d%dengeqtot_e(:,:)=0.0
       pic2d%field2d%dengeqtot_s(:,:)=0.0
       pic2d%field2d%dengeqtot_w(:,:)=0.0
       pic2d%field2d%dengeqtot_n(:,:)=0.0
       pic2d%field2d%dengeqtot_ne(:,:)=0.0
       pic2d%field2d%dengeqtot_se(:,:)=0.0
       pic2d%field2d%dengeqtot_sw(:,:)=0.0
       pic2d%field2d%dengeqtot_nw(:,:)=0.0       
       do i=1,pic2d%para2d%mu_num
         pic2d%field2d%dengeqtot=pic2d%field2d%dengeqtot+pic2d%field2d%deng(i,:,:)
         pic2d%field2d%dengeqtot_e=pic2d%field2d%dengeqtot_e+pic2d%field2d%deng_e(i,:,:)
         pic2d%field2d%dengeqtot_s=pic2d%field2d%dengeqtot_s+pic2d%field2d%deng_s(i,:,:)
         pic2d%field2d%dengeqtot_w=pic2d%field2d%dengeqtot_w+pic2d%field2d%deng_w(i,:,:)
         pic2d%field2d%dengeqtot_n=pic2d%field2d%dengeqtot_n+pic2d%field2d%deng_n(i,:,:)
         pic2d%field2d%dengeqtot_ne=pic2d%field2d%dengeqtot_ne+pic2d%field2d%deng_ne(i,:,:)
         pic2d%field2d%dengeqtot_se=pic2d%field2d%dengeqtot_se+pic2d%field2d%deng_se(i,:,:)
         pic2d%field2d%dengeqtot_sw=pic2d%field2d%dengeqtot_sw+pic2d%field2d%deng_sw(i,:,:)
         pic2d%field2d%dengeqtot_nw=pic2d%field2d%dengeqtot_nw+pic2d%field2d%deng_nw(i,:,:)
       end do
 

     end subroutine
 

     subroutine compute_equdensity_for_ful(pic2d)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       
       pic2d%field2d%denfeq=pic2d%field2d%denf
       pic2d%field2d%denfeq_e=pic2d%field2d%denf_e
       pic2d%field2d%denfeq_s=pic2d%field2d%denf_s
       pic2d%field2d%denfeq_w=pic2d%field2d%denf_w
       pic2d%field2d%denfeq_n=pic2d%field2d%denf_n
       pic2d%field2d%denfeq_ne=pic2d%field2d%denf_ne
       pic2d%field2d%denfeq_se=pic2d%field2d%denf_se
       pic2d%field2d%denfeq_nw=pic2d%field2d%denf_nw
       pic2d%field2d%denfeq_sw=pic2d%field2d%denf_sw 

     end subroutine

     subroutine compute_gyrodensity_perturbation(rootdata,pic2d,pamearray)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       real8,dimension(:,:), pointer :: buf
       real8,dimension(:,:),pointer :: box,re,rs,rw,rn,rne,rse,rsw,rnw
       int4 :: i,num1,num2,row,comm,rank,numproc(2),boxindex(4),ierr     
    
       row=pic2d%para2d%row
       comm=pic2d%layout2d%collective%comm
       rank=pic2d%layout2d%collective%rank
       numproc=pic2d%para2d%numproc
       num1=size(pic2d%field2d%deng_weight,1)
       num2=size(pic2d%field2d%deng_weight,2)
       call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
       allocate(buf(num1,num2))
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2), &
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)

       do i=1,pic2d%para2d%mu_num
         buf=0.0
         box=0.0
         rw=0.0
         re=0.0
         rn=0.0
         rs=0.0
         rsw=0.0
         rse=0.0
         rnw=0.0
         rne=0.0
         buf=pic2d%field2d%deng(i,:,:)
         call solve_weight_of_field_among_processes(buf,rootdata%ASPL,rootdata, &         
              pic2d,box,rw,re,rn,rs,rsw,rse,rnw,rne)
         pic2d%field2d%deng_weight(i,:,:)=box
         pic2d%field2d%dengwg_w(i,:,:)=rw
         pic2d%field2d%dengwg_e(i,:,:)=re
         pic2d%field2d%dengwg_n(i,:,:)=rn
         pic2d%field2d%dengwg_s(i,:,:)=rs
         pic2d%field2d%dengwg_sw(i,:,:)=rsw
         pic2d%field2d%dengwg_se(i,:,:)=rse
         pic2d%field2d%dengwg_nw(i,:,:)=rnw 
         pic2d%field2d%dengwg_ne(i,:,:)=rne
         
         buf=0.0            
         call para_compute_gyroaverage_field_on_mesh(pamearray%mu_nodes(i),i,pic2d,buf, &
              box,rw,re,rn,rs,rsw,rse,rnw,rne)
   
         pic2d%field2d%dengtot=pic2d%field2d%dengtot+buf
       end do
         
         pic2d%field2d%dengtot=pic2d%field2d%dengtot-pic2d%field2d%dengeqtot 
         call mpi2d_alltoallv_box_per_per(row,comm,rank,numproc, &
              pic2d%field2d%dengtot,pic2d%field2d%dengtot_w,pic2d%field2d%dengtot_e,pic2d%field2d%dengtot_n, &
              pic2d%field2d%dengtot_s,pic2d%field2d%dengtot_sw,pic2d%field2d%dengtot_se,pic2d%field2d%dengtot_nw, &
              pic2d%field2d%dengtot_ne,boxindex)

       deallocate(buf,box,re,rn,rs,rw,rse,rsw,rne,rnw)
     end subroutine

!     subroutine compute_fuldensity_perturbation(pic2d)

end module m_fieldsolver
