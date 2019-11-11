module m_fieldsolver
#include "work_precision.h"

use piclayout, only :   root_precompute_data
use m_para_spline, only: para_compute_spl2D_weight
use paradata_type, only: pic_para_total2d_base
use m_mpilayout, only: get_layout_2d_box_index
use m_parautilities, only: gather_field_to_rootprocess_per_per, &
                           mpi2d_alltoallv_box_per_per 
implicit none

public :: solve_weight_of_field_among_processes

contains

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

!    subroutine solve_weight_of_field_among_processes(boxfield,rootmatrix,rootdata,pic2d,boxweight, &
!               bw_w,bw_e,bw_n,bw_s,bw_sw,bw_se,bw_nw,bw_ne)
!      real8,dimension(:,:), pointer, intent(in) :: boxfield
!      real8,dimension(:,:), pointer, intent(in) :: rootmatrix
!      class(root_precompute_data), pointer, intent(inout) :: rootdata
!      class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
!      real8, dimension(:,:), pointer, intent(inout) :: boxweight,bw_w,bw_e,bw_n,bw_s,bw_sw,bw_se,bw_nw,bw_ne
!      int4 :: rank,size, boxindex(4), row, comm
!
!      rank=pic2d%layout2d%collective%rank
!      size=pic2d%layout2d%collective%size
!      row =pic2d%para2d%row
!      comm=pic2d%layout2d%collective%comm
!      call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
!
!      call  gather_field_to_rootprocess_per_per(rootdata%field,boxfield,rank,size,boxindex,&
!            pic2d%para2d%numproc,pic2d%layout2d)
!
!      call  para_compute_spl2D_weight(rootmatrix,rootdata%field,boxweight,pic2d%para2d%numproc, &
!            pic2d%layout2d,pic2d%para2d%boundary)
!
!      call mpi2d_alltoallv_box_per_per(row,comm,rank,pic2d%para2d%numproc, &
!                                        boxweight,bw_w,bw_e,bw_n,bw_s,bw_sw,bw_se,bw_nw,bw_ne,boxindex)
!
!    end subroutine  solve_weight_of_field_among_processes



end module
