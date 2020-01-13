module m_fieldsolver
#include "work_precision.h"

use piclayout, only :   root_precompute_data, parameters_array_2d
use m_para_spline, only: para_compute_spl2D_weight
use paradata_type, only: pic_para_total2d_base
use m_mpilayout, only: get_layout_2d_box_index
use m_parautilities, only: gather_field_to_rootprocess_per_per, &
                           mpi2d_alltoallv_box_per_per, &
                           scatter_field_from_rootprocess_per_per 
use paradata_utilities, only: coords_from_localind
use para_gyroaverage_2d_one, only: para_compute_gyroaverage_mesh_field, &
                                   para_compute_gyroaverage_field_on_mesh
use field_initialize, only: test_trigonfun, test_equdistr
implicit none

public :: solve_weight_of_field_among_processes, &
          solve_field_quasi_neutral, &
          solve_gyfieldweight_from_field, &
          compute_gyrodensity_perturbation, &
          compute_equdensity, &
          compute_equdensity_for_gy, &
          solve_field_ful, &
          compute_equdensity_initdistr_gy
contains

  !!! This subroutine is used to solve the weight of field. First, gather
  !field data to the root process. Then, obtain the weight coefficient. Last,
  !scatter all weight to each process. 
    subroutine solve_weight_of_field_among_processes(boxfield,rootdata,pic2d,boxweight, &
               rw,re,rn,rs,rsw,rse,rnw,rne)
      real8,dimension(:,:), pointer, intent(in) :: boxfield
  !    real8,dimension(:,:), pointer, intent(in) :: rootmatrix
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

      call  para_compute_spl2D_weight(rootdata%ASPL,rootdata%field,boxweight,pic2d%para2d%numproc, &
            pic2d%layout2d,pic2d%para2d%boundary)

      call  mpi2d_alltoallv_box_per_per(row,comm,rank,numproc,boxweight,rw,re,rn,rs,  &
            rsw,rse,rnw,rne,boxindex) 

     end subroutine  solve_weight_of_field_among_processes



!!!! Here, field2d%gep stores the electrostatic potential on the ful-orbit
!spatial space.  
     subroutine solve_field_quasi_neutral(rank,rootdata,pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       int4, intent(in) :: rank
       real8, dimension(:), pointer :: density
       int4 :: dims,global_sz(2),sizeone,boxindex(4),num1,num2
!       real8, dimension(:,:), pointer :: buffer
       int4 :: i,j

       num1=size(pic2d%field2d%dengtot,1)
       num2=size(pic2d%field2d%dengtot,2)
       global_sz(1)=pic2d%layout2d%global_sz1
       global_sz(2)=pic2d%layout2d%global_sz2
       sizeone=pic2d%layout2d%collective%size     
       dims=global_sz(1)*global_sz(2)
       allocate(density(dims))

!!!!!! Here, it may be required to carry out the gyroaverage of
!!!pic2d%field2d%deng
       call get_layout_2d_box_index( pic2d%layout2d, rank, boxindex )       
       call gather_field_to_rootprocess_per_per(density,pic2d%field2d%dengtot,rank, &
           sizeone,boxindex,pic2d%para2d%numproc,pic2d%layout2d)

       if(rank==0) then
          rootdata%field=matmul(rootdata%prematrix,density) 
       end if 


       call scatter_field_from_rootprocess_per_per(rootdata%field,pic2d%field2d%gep,sizeone, &
            pic2d%para2d%numproc,global_sz,pic2d%layout2d)

       deallocate(density)
     end subroutine solve_field_quasi_neutral

     subroutine solve_field_ful(rootdata,pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata 
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       int4 :: i,j,num1,num2,rank
 
       rank=pic2d%layout2d%collective%rank
       num1=size(pic2d%field2d%denf,1)
       num2=size(pic2d%field2d%denf,2)
 
       pic2d%field2d%ep=pic2d%field2d%denf-pic2d%field2d%denfeq
!       if(rank==0) then
!         print*, "denf=", pic2d%field2d%denf
!         print*, 
!         print*, "denfeq=", pic2d%field2d%denfeq 
!       endif
       do i= 1,num1
         do j=1,num2
            pic2d%field2d%ep(i,j)=pic2d%field2d%ep(i,j)*pic2d%para2d%temp_e(i,j)/ &
                                    pic2d%field2d%denfeq(i,j) 
         end do
       enddo
       
         
     end subroutine

     subroutine solve_gyfieldweight_from_field(rootdata, pic2d,pamearray)
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       real8,dimension(:,:),pointer :: buf,box,re,rs,rw,rn,rne,rse,rsw,rnw 
       int4 :: i,ierr,num1,num2,row, rank
 
       rank=pic2d%layout2d%collective%rank
       num1=size(pic2d%field2d%ep,1)
       num2=size(pic2d%field2d%ep,2)
       row=pic2d%para2d%row
       allocate(buf(num1,num2))       
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2), &
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)
 
       pic2d%field2d%epgyro=0.0       
       do i=1,pic2d%para2d%mu_num

         call para_compute_gyroaverage_mesh_field(pamearray%mu_nodes(i),i,pic2d)
 
!if(rank==0) then
!if(i==1.or.i==5) then
!print*, pic2d%field2d%epgyro(i,:,:)
!endif
!endif
         buf=pic2d%field2d%epgyro(i,:,:)

         box=0.0
         rw=0.0
         re=0.0
         rn=0.0
         rs=0.0
         rsw=0.0
         rse=0.0
         rnw=0.0
         rne=0.0

         call solve_weight_of_field_among_processes(buf,rootdata,pic2d, &
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
!if(rank==0) then
!print*, "i=1",pic2d%field2d%epgyro(1,:,:)
!print*, 
!print*, "i=12",pic2d%field2d%epgyro(12,:,:)
!endif   
       deallocate(buf,box,rw,re,rn,rs,rsw,rse,rnw,rne)
     end subroutine solve_gyfieldweight_from_field


     subroutine compute_equdensity_for_gy(pic2d)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       int4 :: i,rank      
!       pic2d%field2d%deng=pic2d%field2d%deng-pic2d%field2d%dengeq
        rank=pic2d%layout2d%collective%rank
        pic2d%field2d%dengeqtot(:,:)=0.0
!       pic2d%field2d%dengeqtot_e(:,:)=0.0
!       pic2d%field2d%dengeqtot_s(:,:)=0.0
!       pic2d%field2d%dengeqtot_w(:,:)=0.0
!       pic2d%field2d%dengeqtot_n(:,:)=0.0
!       pic2d%field2d%dengeqtot_ne(:,:)=0.0
!       pic2d%field2d%dengeqtot_se(:,:)=0.0
!       pic2d%field2d%dengeqtot_sw(:,:)=0.0
!       pic2d%field2d%dengeqtot_nw(:,:)=0.0       
       do i=1,pic2d%para2d%mu_num
!         if(rank==0) then
!          print*, "i=",i,pic2d%field2d%dengeq(i,1,1)
!         endif
         pic2d%field2d%dengeqtot=pic2d%field2d%dengeqtot+pic2d%field2d%dengeq(i,:,:)
!         pic2d%field2d%dengeqtot_e=pic2d%field2d%dengeqtot_e+pic2d%field2d%deng_e(i,:,:)
!         pic2d%field2d%dengeqtot_s=pic2d%field2d%dengeqtot_s+pic2d%field2d%deng_s(i,:,:)
!         pic2d%field2d%dengeqtot_w=pic2d%field2d%dengeqtot_w+pic2d%field2d%deng_w(i,:,:)
!         pic2d%field2d%dengeqtot_n=pic2d%field2d%dengeqtot_n+pic2d%field2d%deng_n(i,:,:)
!         pic2d%field2d%dengeqtot_ne=pic2d%field2d%dengeqtot_ne+pic2d%field2d%deng_ne(i,:,:)
!         pic2d%field2d%dengeqtot_se=pic2d%field2d%dengeqtot_se+pic2d%field2d%deng_se(i,:,:)
!         pic2d%field2d%dengeqtot_sw=pic2d%field2d%dengeqtot_sw+pic2d%field2d%deng_sw(i,:,:)
!         pic2d%field2d%dengeqtot_nw=pic2d%field2d%dengeqtot_nw+pic2d%field2d%deng_nw(i,:,:)
       end do
 

     end subroutine compute_equdensity_for_gy
 

     subroutine compute_equdensity(numparticles,box,re,rs,rw,rn,rne,rse,rsw,rnw,pic2d,eqdistr,statkind)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       real(8), dimension(:,:),pointer :: box,re,rs,rw,rn,rne,rse,rsw,rnw
       int4, intent(in) :: numparticles
       real8 :: gxmin(2), delta(2),coords(2),mean(2)
       int4 :: ND(2),rank,comm,boxindex(4),row,numproc(2),globind(2)
       int4 :: i,j     
       real8 :: summ

       character(len=*), optional,intent(in) :: eqdistr
       int4, intent(in) :: statkind
 
       comm=pic2d%layout2d%collective%comm
       rank=pic2d%layout2d%collective%rank
       delta(1)=pic2d%para2d%m_x1%delta_eta
       delta(2)=pic2d%para2d%m_x2%delta_eta
       ND(1)=pic2d%para2d%m_x1%nodes
       ND(2)=pic2d%para2d%m_x2%nodes
       globind(1)=pic2d%layout2d%global_sz1
       globind(2)=pic2d%layout2d%global_sz2
       gxmin=pic2d%para2d%gxmin
       row=pic2d%para2d%row
       numproc=NINT(sqrt(real(pic2d%layout2d%collective%size,8)))
       mean=(pic2d%para2d%gxmin+pic2d%para2d%gxmax)/2.0
       summ = 0.0  

!!!!! Here, MPI_ALLREDUCE is better.      
       select case(statkind)
         case(1)
         call get_layout_2d_box_index(pic2d%layout2d, rank,boxindex)
         do j=1,globind(2)
           coords(2)=gxmin(2)+delta(2)*real(j-1,8)
           do i=1,globind(1)
             coords(1)=gxmin(1)+delta(1)*real(i-1,8)
!             coords=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)
             summ=summ+test_equdistr(coords,pic2d%para2d%sigma,mean,eqdistr)*delta(1)*delta(2)
           enddo
         enddo

         do j=1,ND(2)
           do i=1,ND(1)
             coords=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)           
             box(i,j)=test_equdistr(coords,pic2d%para2d%sigma,mean,eqdistr)*delta(1) &
                    *delta(2)/summ*real(numparticles,8)
           end do
         enddo

         
       case(2)
         box = pic2d%field2d%denf         
       case default
         stop
       end select
 
!print*,1

       call mpi2d_alltoallv_box_per_per(row,comm,rank,numproc, &
                                        box,rw,re,rn,rs,rsw,rse,rnw,rne,boxindex)
       
     end subroutine compute_equdensity

     subroutine compute_equdensity_initdistr_gy(mu_num,munum_partition,eqdistr,statkind,pic2d)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       int4, dimension(:), intent(in) :: munum_partition
       int4, intent(in) :: mu_num
       character(len=*), optional,intent(in) :: eqdistr
       int4, intent(in) :: statkind
       int4 :: rank,num1,num2,row, sizeone
       int4 :: i,j,ierr

       real8,dimension(:,:),pointer :: box,re,rs,rw,rn,rne,rse,rsw,rnw 

       rank=pic2d%layout2d%collective%rank
       num1=size(pic2d%field2d%ep,1)
       num2=size(pic2d%field2d%ep,2)
       row=pic2d%para2d%row
       sizeone=pic2d%layout2d%collective%size
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2), &
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)       
       
       do i=1, mu_num

         box=0.0
         rw=0.0
         re=0.0
         rn=0.0
         rs=0.0
         rsw=0.0
         rse=0.0
         rnw=0.0
         rne=0.0
         
         call compute_equdensity(sizeone*munum_partition(i),box,re,rs,rw,rn,rne,rse,rsw,rnw, &
                pic2d,eqdistr,statkind)

!   if(rank==0) then
!     print*, "i=",i,rne
!   endif
           pic2d%field2d%dengeq(i,:,:)=box
           pic2d%field2d%dengeq_e(i,:,:)=re 
           pic2d%field2d%dengeq_s(i,:,:)=rs
           pic2d%field2d%dengeq_w(i,:,:)=rw
           pic2d%field2d%dengeq_n(i,:,:)=rn
           pic2d%field2d%dengeq_ne(i,:,:)=rne
           pic2d%field2d%dengeq_se(i,:,:)=rse
           pic2d%field2d%dengeq_sw(i,:,:)=rsw
           pic2d%field2d%dengeq_nw(i,:,:)=rnw
 
       end do

       deallocate(box,re,rs,rw,rn,rne,rse,rsw,rnw)
     end subroutine  compute_equdensity_initdistr_gy

     subroutine compute_gyrodensity_perturbation(rootdata,pic2d,pamearray)
       class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
       class(root_precompute_data), pointer, intent(inout) :: rootdata
       class(parameters_array_2d), pointer, intent(in) :: pamearray
       real8,dimension(:,:), pointer :: buf
       real8,dimension(:,:),pointer :: box,re,rs,rw,rn,rne,rse,rsw,rnw
       int4 :: num1,num2,row,comm,rank,numproc(2),boxindex(4),ierr     
       int4 :: ND(2),i,j   
 
       row=pic2d%para2d%row
       comm=pic2d%layout2d%collective%comm
       rank=pic2d%layout2d%collective%rank
       numproc=pic2d%para2d%numproc
       num1=size(pic2d%field2d%deng_weight,2)
       num2=size(pic2d%field2d%deng_weight,3)
       ND(1)=pic2d%para2d%m_x1%nodes
       ND(2)=pic2d%para2d%m_x2%nodes
       call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
       allocate(buf(num1,num2))
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2), &
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)

       pic2d%field2d%dengtot = 0.0
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

         call solve_weight_of_field_among_processes(buf,rootdata, &         
              pic2d,box,rw,re,rn,rs,rsw,rse,rnw,rne)

!         pic2d%field2d%deng_weight(i,:,:)=box
!         pic2d%field2d%dengwg_w(i,:,:)=rw
!         pic2d%field2d%dengwg_e(i,:,:)=re
!         pic2d%field2d%dengwg_n(i,:,:)=rn
!         pic2d%field2d%dengwg_s(i,:,:)=rs
!         pic2d%field2d%dengwg_sw(i,:,:)=rsw
!         pic2d%field2d%dengwg_se(i,:,:)=rse
!         pic2d%field2d%dengwg_nw(i,:,:)=rnw 
!         pic2d%field2d%dengwg_ne(i,:,:)=rne
         
         buf=0.0            
         call para_compute_gyroaverage_field_on_mesh(pamearray%mu_nodes(i),i,pic2d,buf, &
              box,rw,re,rn,rs,rsw,rse,rnw,rne)
   
         pic2d%field2d%dengtot=pic2d%field2d%dengtot+buf
       end do
       
       do j=1,ND(2)
         do i=1,ND(1)
          pic2d%field2d%dengtot(i,j)=(pic2d%field2d%dengtot(i,j)-pic2d%field2d%dengeqtot(i,j)) &
                                      /pic2d%field2d%dengeqtot(i,j) 
         end do 
       end do 

        call mpi2d_alltoallv_box_per_per(row,comm,rank,numproc, &
              pic2d%field2d%dengtot,pic2d%field2d%dengtot_w,pic2d%field2d%dengtot_e,pic2d%field2d%dengtot_n, &
              pic2d%field2d%dengtot_s,pic2d%field2d%dengtot_sw,pic2d%field2d%dengtot_se,pic2d%field2d%dengtot_nw, &
              pic2d%field2d%dengtot_ne,boxindex)

       deallocate(buf,box,re,rn,rs,rw,rse,rsw,rne,rnw)
     end subroutine compute_gyrodensity_perturbation

!     subroutine compute_fuldensity_perturbation(pic2d)

end module m_fieldsolver
