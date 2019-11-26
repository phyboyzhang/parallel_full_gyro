module tp_preparation
#include "work_precision.h"
use constants, only: pi_
use paradata_layout, only:    initialize_pic_para_2d_base
use paradata_type, only: pic_para_total2d_base
use orbit_data_base, only: rk4ful2dnode, pointin_node,tp_ful2d_node, &
                           tp_ful2dsend_node, &
                           rk4gy2dnode, tp_gy2dallmu_node,tp_gy2dsend_node
use m_picutilities, only: mpi2d_alltoallv_send_particle_2d,&
                          tp_mpi2d_alltoallv_send_particle_2d, &
                          tp_mpi2d_alltoallv_send_particle_2d_gy
use paradata_utilities, only: compute_process_of_point_per_per
implicit none
include "mpif.h"

contains

 subroutine initialize_test_particles(tpful2d_head,tprk4ful2d_head,tpgy2dmu_head, & 
            numleft_boris,numleft_rk4,numleft_gy,numgr,numgr_gy,pic2d)
   class(pic_para_total2d_base),pointer,intent(inout) :: pic2d
   class(tp_ful2d_node), pointer, intent(inout) :: tpful2d_head,tprk4ful2d_head
   class(tp_ful2d_node), pointer :: tpful2dtmp,tprk4ful2dtmp
   class(tp_gy2dallmu_node),dimension(:), pointer, intent(inout) :: tpgy2dmu_head
   int4, intent(inout) :: numleft_boris,numleft_rk4
   int4, dimension(:), pointer, intent(inout) :: numleft_gy
   class(tp_gy2dallmu_node),dimension(:), pointer :: tpgy2dmutmp 
   class(tp_ful2dsend_node), dimension(:), pointer :: tpful2dsend_head,tpful2dsendtmp, &
                             tprk4ful2dsend_head, tprk4ful2dsendtmp
   class(tp_gy2dsend_node), dimension(:), pointer :: tpgy2dsend_head,tpgy2dsendtmp
   int4,intent(in) :: numgr,numgr_gy
   int4,dimension(:), pointer :: num_boris,num_rk4, num_gy   
   int4 :: size,rank,i,j,numcircle,rank1,mu_num
   real8 :: rho,x1(2),theta,coords(4)  
   
   mu_num=pic2d%para2d%mu_num
   rank=pic2d%layout2d%collective%rank
   size=pic2d%layout2d%collective%size
 
   allocate(tpful2dsend_head(0:size-1),tpful2dsendtmp(0:size-1),tprk4ful2dsend_head(0:size-1), &
           tprk4ful2dsendtmp(0:size-1),tpgy2dsendtmp(0:size-1))
   allocate(num_boris(0:size-1),num_rk4(0:size-1),num_gy(0:size-1))
   allocate(tpgy2dmutmp(mu_num))
   do i=1,mu_num
     tpgy2dmutmp(i)%ptr=>tpgy2dmu_head(i)%ptr
   end do

   do i=0, size-1
     allocate(tpful2dsend_head(i)%ptr)
     allocate(tprk4ful2dsend_head(i)%ptr)
     tpful2dsendtmp(i)%ptr=>tpful2dsend_head(i)%ptr
     tprk4ful2dsendtmp(i)%ptr=>tprk4ful2dsend_head(i)%ptr
   end do

   numcircle=pic2d%para2d%numcircle
   tpful2dtmp=>tpful2d_head
   tprk4ful2dtmp=>tprk4ful2d_head
!   numgr=5
!   numgr_gy=4
   num_boris=0
   num_rk4=0

  do i=1,mu_num
    allocate(tpgy2dsend_head(0:size-1))  
    do j=0,size-1
      allocate(tpgy2dsend_head(i)%ptr)
      tpgy2dsendtmp(i)%ptr=>tpgy2dsend_head(i)%ptr 
    end do 
  num_gy=0
 
  if(rank==0) then
   x1(1)=(pic2d%para2d%gxmin(1)+pic2d%para2d%gxmax(1))/2.0+0.5
   x1(2)=(pic2d%para2d%gxmin(2)+pic2d%para2d%gxmax(2))/2.0
   rho=sqrt(2._f64*pic2d%para2d%mu)

   rank1=compute_process_of_point_per_per(x1(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
         pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
       if(rank1==rank) then
         tpgy2dmutmp(i)%ptr%coords(1:2)=x1(1:2)
         tpgy2dmutmp(i)%ptr%coords(3)  =pic2d%para2d%mu
         tpgy2dmutmp(i)%ptr%tp=1
         num_gy(rank1)=num_gy(rank1)+1
         allocate(tpgy2dmutmp(i)%ptr%next)
         tpgy2dmutmp(i)%ptr=>tpgy2dmutmp(i)%ptr%next
       else
         tpgy2dsendtmp(rank1)%ptr%coords(1:2)=x1(1:2)
         tpgy2dsendtmp(rank1)%ptr%coords(3)  =pic2d%para2d%mu
         tpgy2dsendtmp(rank1)%ptr%tp=1
         num_gy(rank1)=num_gy(rank1)+1
         allocate(tpgy2dsendtmp(rank1)%ptr%next)
         tpgy2dsendtmp(rank1)%ptr=>tpgy2dsendtmp(rank1)%ptr%next
       endif

   end if
       call tp_mpi2d_alltoallv_send_particle_2d_gy(tpgy2dmu_head,numleft_gy(i),tpgy2dsend_head, &
                                                num_gy,numgr_gy,1,pic2d)
       deallocate(tpgy2dsend_head)
  if(rank==0) then
    do j=1,numcircle-1
       theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)
       coords(1)=x1(1)+rho*cos(theta)
       coords(2)=x1(2)+rho*sin(theta)
       coords(3)=-rho*cos(theta+pi_/2.0_f64)
       coords(4)=-rho*sin(theta+pi_/2.0_f64)
!print*, "coords=",coords(1:4)
!       call coordinates_pointoutbound_per_per(coords(1:2),pic2d)
       rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
             pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
       if(rank1==rank) then
         tpful2dtmp%coords(1:4)=coords(1:4)
         tpful2dtmp%tp=j+1
         num_boris(rank1)=num_boris(rank1)+1
         allocate(tpful2dtmp%next)
         tpful2dtmp=>tpful2dtmp%next

         tprk4ful2dtmp%coords(1:4)=coords(1:4)
         tprk4ful2dtmp%tp=j+1
         num_rk4(rank1)=num_rk4(rank1)+1
         allocate(tprk4ful2dtmp%next)
         tprk4ful2dtmp=>tprk4ful2dtmp%next

       else
         tpful2dsendtmp(rank1)%ptr%coords(1:4)=coords(1:4)
         tpful2dsendtmp(rank1)%ptr%tp=j+1
         num_boris(rank1)=num_boris(rank1)+1
         allocate(tpful2dsendtmp(rank1)%ptr%next)
         tpful2dsendtmp(rank1)%ptr=>tpful2dsendtmp(rank1)%ptr%next

         tprk4ful2dsendtmp(rank1)%ptr%coords(1:4)=coords(1:4)
         tprk4ful2dsendtmp(rank1)%ptr%tp=j+1
         num_rk4(rank1)=num_rk4(rank1)+1
         allocate(tprk4ful2dsendtmp(rank1)%ptr%next)
         tprk4ful2dsendtmp(rank1)%ptr=>tprk4ful2dsendtmp(rank1)%ptr%next

       endif
     end do
   end if

  end do   !!!! end mu_num

    call tp_mpi2d_alltoallv_send_particle_2d(tpful2d_head,numleft_boris,tpful2dsend_head,num_boris,numgr,pic2d)
    call tp_mpi2d_alltoallv_send_particle_2d(tprk4ful2d_head,numleft_rk4,tprk4ful2dsend_head,num_rk4,numgr,pic2d)


 end subroutine

 
end module tp_preparation 
