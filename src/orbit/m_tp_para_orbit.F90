module m_tp_para_orbit
#include "work_precision.h"
use orbit_data_base, only: tp_ful2d_node, tp_gy2d_node, &
                           tp_ful2dsend_node,tp_gy2dsend_node    
use piclayout, only: ful2d_node, gy2d_node
use paradata_utilities,only: pic_para_total2d_base
use m_para_orbit, only: borissolve,fulrk4solve,gyrork4solve
use m_picutilities, only: tp_sort_particles_among_ranks, &
                          tp_mpi2d_alltoallv_send_particle_2d, &
                          tp_sort_particles_among_ranks_gy, &
                          tp_mpi2d_alltoallv_send_particle_2d_gy                 
 
implicit none
include "mpif.h"

public :: tp_push_ful_orbit, &
          tp_push_gy_orbit

contains

    subroutine tp_push_ful_orbit(tp_ful2d_head,numleft,numgr,pic2d,pushkind,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
       int4, intent(inout) :: numleft
       int4, intent(in) :: numgr,iter_num,rk4order
       character(len=*),intent(in) :: pushkind
       int4, dimension(:), pointer :: num
       int4 :: size,rank,comm
       class(tp_ful2dsend_node),dimension(:), pointer :: tp_ful2dsend_head 
       int4 :: i,nump
 
       class(tp_ful2d_node), pointer :: tpful2dtmp 
!       class(tp_ful2dsend_node), dimension(:), pointer :: tpful2dsendtmp 
    
       rank=pic2d%layout2d%collective%rank      
       size=pic2d%layout2d%collective%size 
       comm=pic2d%layout2d%collective%comm
       allocate(num(0:size-1),tp_ful2dsend_head(0:size-1))
       num=0 
       do i=0, size-1
         allocate(tp_ful2dsend_head(i)%ptr)
       end do

!if(iter_num==4) then
!tpful2dtmp=>tp_ful2d_head
!do while(associated(tpful2dtmp))
!   if(.not.associated(tpful2dtmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",tpful2dtmp%coords(1:4)
!      tpful2dtmp=>tpful2dtmp%next
!   end if
!end do
!end if
 
        call tp_ful_solve(tp_ful2d_head,pic2d,numleft,pushkind,rk4order,iter_num)

        call tp_sort_particles_among_ranks(tp_ful2d_head,tp_ful2dsend_head,pic2d,num)

        call tp_mpi2d_alltoallv_send_particle_2d(tp_ful2d_head,numleft,tp_ful2dsend_head,num,numgr,pic2d)

       deallocate(tp_ful2dsend_head,num)
!       nullify(tpful2dsendtmp)
    end subroutine tp_push_ful_orbit

    subroutine tp_push_gy_orbit(tp_gy2d_head,numleft,numgr,pic2d,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_gy2d_node), pointer, intent(inout) :: tp_gy2d_head
       int4, intent(inout) :: numleft
       int4, intent(in) :: numgr,iter_num,rk4order
!       character(len=*),intent(in) :: pushkind
       int4, dimension(:), pointer :: num
       int4 :: size,rank,comm
       class(tp_gy2dsend_node),dimension(:), pointer :: tp_gy2dsend_head 
       int4 :: i,nump
 
       class(tp_gy2d_node), pointer :: tpgy2dtmp 
   
       rank=pic2d%layout2d%collective%rank      
       size=pic2d%layout2d%collective%size 
       comm=pic2d%layout2d%collective%comm
       allocate(num(0:size-1),tp_gy2dsend_head(0:size-1))
       num=0 
       do i=0, size-1
         allocate(tp_gy2dsend_head(i)%ptr)
       end do


       call tp_gy_solve(tp_gy2d_head,pic2d,numleft,rk4order,iter_num)

!if(iter_num==4) then
!tpgy2dtmp=>tp_gy2d_head
!do while(associated(tpgy2dtmp)) 
!   if(.not.associated(tpgy2dtmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",tpgy2dtmp%coords(1:3)
!      tpgy2dtmp=>tpgy2dtmp%next
!   end if
!end do
!end if
!print*, rank
!call mpi_barrier(comm)


       call tp_sort_particles_among_ranks_gy(tp_gy2d_head,tp_gy2dsend_head,pic2d,num)
!if(iter_num==6) then
!   print*, "rank=",rank,"num=",num    
!tpful2dtmp=>tp_ful2d_head
!do while(associated(tpful2dtmp)) 
!   if(.not.associated(tpful2dtmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",tpful2dtmp%coords(1:4)
!      tpful2dtmp=>tpful2dtmp%next
!   end if
!end do
!
!allocate(tpful2dsendtmp(0:size-1))
!do i=0,size-1
!   tpful2dsendtmp(i)%ptr=>tp_ful2dsend_head(i)%ptr
!do while(associated(tpful2dsendtmp(i)%ptr))
!   if(.not.associated(tpful2dsendtmp(i)%ptr%next)) then
!      exit
!   else
!      print*, "rank2=",rank,"sendcoords1=",tpful2dsendtmp(i)%ptr%coords(1:4)
!      tpful2dsendtmp(i)%ptr=>tpful2dsendtmp(i)%ptr%next
!   end if
!end do
!end do
!end if

       call tp_mpi2d_alltoallv_send_particle_2d_gy(tp_gy2d_head,numleft,tp_gy2dsend_head,num,numgr,pic2d)
!if(iter_num==3) then
!tpful2dtmp=>tp_ful2d_head
!do while(associated(tpful2dtmp))
!   if(.not.associated(tpful2dtmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",tpful2dtmp%coords(1:4)
!      tpful2dtmp=>tpful2dtmp%next
!   end if
!end do
!end if

       deallocate(tp_gy2dsend_head,num)
    end subroutine tp_push_gy_orbit



    subroutine tp_ful_solve(tp_ful2d_head,pic2d,numleft,pushkind,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
       int4, intent(in) :: numleft,iter_num,rk4order
       character(len=*),intent(in) :: pushkind
       class(ful2d_node), pointer :: ful2d_head, ful2dtmp      
       int4 :: rank,comm  

       class(tp_ful2d_node), pointer :: tpful2dtmp 
       rank=pic2d%layout2d%collective%rank
       comm=pic2d%layout2d%collective%comm
       allocate(ful2d_head)
       if(.not.associated(tp_ful2d_head)) then
          print*, "rank=",rank,"#ERROR: tp_ful2d_head is not allocated."
          stop           
       end if 

       call tp_ful2dlist_to_ful2dlist(ful2d_head,tp_ful2d_head)
!if(iter_num==2) then
!tpful2dtmp=>tp_ful2d_head
!do while(associated(tpful2dtmp))
!   if(.not.associated(tpful2dtmp%next)) then
!      exit
!   else
!      print*, "rank2=",rank,"coords1=",tpful2dtmp%coords(1:4)
!      tpful2dtmp=>tpful2dtmp%next
!   end if
!end do
!end if

       if(pushkind=="boris") then
         call borissolve(ful2d_head,pic2d,numleft)
       else if(pushkind=="rk4") then
         call fulrk4solve(ful2d_head,pic2d,rk4order,iter_num)
       end if

       call ful2dlist_to_tp_ful2dlist(tp_ful2d_head,ful2d_head,rank)
 
       deallocate(ful2d_head)
       nullify(ful2d_head)
    end subroutine tp_ful_solve

    subroutine tp_gy_solve(tp_gy2d_head,pic2d,numleft,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_gy2d_node), pointer, intent(inout) :: tp_gy2d_head
       int4, intent(in) :: numleft,iter_num,rk4order
!       character(len=*),intent(in) :: pushkind
       class(gy2d_node), pointer :: gy2d_head, gy2dtmp      
       int4 :: rank,comm  

       class(tp_gy2d_node), pointer :: tpgy2dtmp 
       rank=pic2d%layout2d%collective%rank
       comm=pic2d%layout2d%collective%comm
       allocate(gy2d_head)
       if(.not.associated(tp_gy2d_head)) then
          print*, "rank=",rank,"#ERROR: tp_ful2d_head is not allocated."
          stop           
       end if 

       call tp_gy2dlist_to_gy2dlist(gy2d_head,tp_gy2d_head)
       call gyrork4solve(gy2d_head,pic2d,iter_num)
       call gy2dlist_to_tp_gy2dlist(tp_gy2d_head,gy2d_head,rank)
       deallocate(gy2d_head)
       nullify(gy2d_head)
    end subroutine tp_gy_solve



    subroutine tp_ful2dlist_to_ful2dlist(ful2d_head,tp_ful2d_head)
       class(tp_ful2d_node), pointer, intent(in) :: tp_ful2d_head
       class(tp_ful2d_node), pointer :: tpful2dtmp
       class(ful2d_node), pointer,intent(inout) :: ful2d_head
       class(ful2d_node), pointer :: ful2dtmp

       if(.NOT.associated(ful2d_head)) then
          print*, "#ERROR: ful2d_head is not allocated."
          stop
       end if
       if(.NOT.associated(tp_ful2d_head)) then
          print*, "#ERROR: tp_ful2d_head is not allocated."
          stop
       end if
       ful2dtmp=>ful2d_head
       tpful2dtmp=>tp_ful2d_head

       do while(associated(tpful2dtmp))
          if(.not.associated(tpful2dtmp%next)) then
            exit
          else
            ful2dtmp%coords(1:4)=tpful2dtmp%coords(1:4)
            allocate(ful2dtmp%next)
            ful2dtmp=>ful2dtmp%next
            tpful2dtmp=>tpful2dtmp%next
          end if
       end do
      
       nullify(ful2dtmp)
       nullify(tpful2dtmp)
    end subroutine tp_ful2dlist_to_ful2dlist

    subroutine tp_gy2dlist_to_gy2dlist(gy2d_head,tp_gy2d_head)
       class(tp_gy2d_node), pointer, intent(in) :: tp_gy2d_head
       class(tp_gy2d_node), pointer :: tpgy2dtmp
       class(gy2d_node), pointer,intent(inout) :: gy2d_head
       class(gy2d_node), pointer :: gy2dtmp

       if(.NOT.associated(gy2d_head)) then
          print*, "#ERROR: ful2d_head is not allocated."
          stop
       end if
       if(.NOT.associated(tp_gy2d_head)) then
          print*, "#ERROR: tp_ful2d_head is not allocated."
          stop
       end if
       gy2dtmp=>gy2d_head
       tpgy2dtmp=>tp_gy2d_head

       do while(associated(tpgy2dtmp))
          if(.not.associated(tpgy2dtmp%next)) then
            exit
          else
            gy2dtmp%coords(1:3)=tpgy2dtmp%coords(1:3)
            allocate(gy2dtmp%next)
            gy2dtmp=>gy2dtmp%next
            tpgy2dtmp=>tpgy2dtmp%next
          end if
       end do
      
       nullify(gy2dtmp)
       nullify(tpgy2dtmp)
    end subroutine



    subroutine ful2dlist_to_tp_ful2dlist(tp_ful2d_head,ful2d_head,rank)
       class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
       class(tp_ful2d_node), pointer :: tpful2dtmp
       class(ful2d_node), pointer, intent(in) :: ful2d_head
       class(ful2d_node), pointer :: ful2dtmp
       int4, intent(in) :: rank

       if(.NOT.associated(ful2d_head)) then
          print*, "#ERROR: ful2d_head is not allocated."
          stop
       end if
       if(.NOT.associated(tp_ful2d_head)) then
          print*,  "rank=",rank,"#ERROR: tp_ful2d_head is not allocated."
          stop
       end if
       ful2dtmp=>ful2d_head
       tpful2dtmp=>tp_ful2d_head

       do while(associated(ful2dtmp)) 
          if(.not.associated(ful2dtmp%next)) then
            exit
          else
            tpful2dtmp%coords(1:4)=ful2dtmp%coords(1:4)
   !         allocate(tpful2dtmp%next)
            tpful2dtmp=>tpful2dtmp%next
            ful2dtmp=>ful2dtmp%next
          end if
       end do
       nullify(ful2dtmp)
       nullify(tpful2dtmp)
    end subroutine

  subroutine gy2dlist_to_tp_gy2dlist(tp_gy2d_head,gy2d_head,rank)
       class(tp_gy2d_node), pointer, intent(inout) :: tp_gy2d_head
       class(tp_gy2d_node), pointer :: tpgy2dtmp
       class(gy2d_node), pointer, intent(in) :: gy2d_head
       class(gy2d_node), pointer :: gy2dtmp
       int4, intent(in) :: rank

       if(.NOT.associated(gy2d_head)) then
          print*, "#ERROR: ful2d_head is not allocated."
          stop
       end if
       if(.NOT.associated(tp_gy2d_head)) then
          print*,  "rank=",rank,"#ERROR: tp_ful2d_head is not allocated."
          stop
       end if
       gy2dtmp=>gy2d_head
       tpgy2dtmp=>tp_gy2d_head

       do while(associated(gy2dtmp))
          if(.not.associated(gy2dtmp%next)) then
            exit
          else
            tpgy2dtmp%coords(1:3)=gy2dtmp%coords(1:3)
   !         allocate(tpful2dtmp%next)
            tpgy2dtmp=>tpgy2dtmp%next
            gy2dtmp=>gy2dtmp%next
          end if
       end do
       nullify(gy2dtmp)
       nullify(tpgy2dtmp)   
   end subroutine gy2dlist_to_tp_gy2dlist   

 end module m_tp_para_orbit
