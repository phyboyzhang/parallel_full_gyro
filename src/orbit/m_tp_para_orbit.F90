module m_tp_para_orbit
#include "work_precision.h"
use orbit_data_base, only: tp_ful2d_node, tp_gy2d_node, &
                           tp_ful2dsend_node,tp_gy2dsend_node, &
                           tp_gy2dallmu_node    
use piclayout, only: ful2d_node, gy2d_node, gy2dmu_node
use paradata_utilities,only: pic_para_total2d_base
use m_para_orbit, only: borissolve,fulrk4solve,gyrork4solveallmu
use m_picutilities, only: tp_sort_particles_among_ranks, &
                          tp_mpi2d_alltoallv_send_particle_2d, &
                          tp_sort_particles_among_ranks_gy, &
                          tp_mpi2d_alltoallv_send_particle_2d_gy                 
 
implicit none
include "mpif.h"

public :: tp_push_ful_orbit, &
          tp_push_gy_orbit_allmu

contains

    subroutine tp_push_ful_orbit(tp_ful2d_head,numleft,numgr,pic2d,pushkind,iter_num, rk4order)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
       int4, intent(inout) :: numleft
       int4, intent(in) :: numgr,iter_num
       int4, optional,intent(in) :: rk4order
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

       call tp_ful_solve(tp_ful2d_head,pic2d,numleft,pushkind,rk4order,iter_num)
       call tp_sort_particles_among_ranks(tp_ful2d_head,tp_ful2dsend_head,pic2d,num)
        call tp_mpi2d_alltoallv_send_particle_2d(tp_ful2d_head,numleft,tp_ful2dsend_head,num,numgr,pic2d)
 
       deallocate(tp_ful2dsend_head,num)
!       nullify(tpful2dsendtmp)
    end subroutine tp_push_ful_orbit

    subroutine tp_push_gy_orbit_allmu(tpgy2dmu_head,numleft,numgr,pic2d,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_gy2dallmu_node),dimension(:), pointer, intent(inout) :: tpgy2dmu_head
       int4,dimension(:), intent(inout) :: numleft
       int4, intent(in) :: numgr,iter_num,rk4order
!       character(len=*),intent(in) :: pushkind
       class(tp_gy2dsend_node),dimension(:), pointer :: tp_gy2dsend_head 
       class(tp_gy2d_node), pointer :: tpgy2d_head, tpgy2dtmp
       class(tp_gy2dallmu_node), dimension(:), pointer :: tpgy2dmutmp
       int4, dimension(:), pointer :: num
       int4 :: size,rank,comm 
       int4 :: i,j,nump,mu_num
   
       rank=pic2d%layout2d%collective%rank      
       size=pic2d%layout2d%collective%size 
       comm=pic2d%layout2d%collective%comm
       mu_num=pic2d%para2d%mu_num
       allocate(num(0:size-1))
       allocate(tpgy2dmutmp(1:mu_num))
   
       do i=1,mu_num
         tpgy2dmutmp(i)%ptr=>tpgy2dmu_head(i)%ptr
       end do

       call tp_gy_solve_allmu(tpgy2dmu_head,pic2d,rk4order,iter_num)
         call mpi_barrier(comm)

      do i=1,mu_num
         allocate(tp_gy2dsend_head(0:size-1))
         num=0
         do j=0, size-1
           allocate(tp_gy2dsend_head(j)%ptr)
         end do 

         call tp_sort_particles_among_ranks_gy(tpgy2dmu_head,tp_gy2dsend_head,i,pic2d,num)
 
         call tp_mpi2d_alltoallv_send_particle_2d_gy(tpgy2dmu_head,numleft(i),tp_gy2dsend_head,num,numgr,mu_num,i,pic2d)

         deallocate(tp_gy2dsend_head)
       end do

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

       deallocate(num,tpgy2dmutmp)
    end subroutine tp_push_gy_orbit_allmu



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

!if(iter_num==1) then
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

    subroutine tp_gy_solve_allmu(tp_gy2dmu_head,pic2d,rk4order,iter_num)
       class(pic_para_total2d_base), pointer :: pic2d
       class(tp_gy2dallmu_node),dimension(:), pointer, intent(inout) :: tp_gy2dmu_head
       int4, intent(in) :: iter_num,rk4order
!       character(len=*),intent(in) :: pushkind
       class(gy2dmu_node), dimension(:),pointer :: gy2dmu_head     
!       class(tp_gy2dallmu_node),dimension(:), pointer :: tpgy2dtmp 
       int4 :: rank,comm,mu_num,i  
 
       rank=pic2d%layout2d%collective%rank
       comm=pic2d%layout2d%collective%comm
       mu_num=pic2d%para2d%mu_num
       allocate(gy2dmu_head(1:mu_num))  !,gy2dmutmpt(1:mu_num), tpgy2dtmp(1:mu_num))
       do i=1, mu_num
         allocate(gy2dmu_head(i)%ptr)
       end do
       do i=1,mu_num
         if(.not.associated(tp_gy2dmu_head(i)%ptr)) then
           print*, "rank=",rank,"#ERROR: tp_gy2dallmu_head is not allocated."
           stop           
         end if 
       end do 

       call tp_gy2dlist_to_gy2dlist_allmu(gy2dmu_head,tp_gy2dmu_head,mu_num)
       call gyrork4solveallmu(gy2dmu_head,pic2d,iter_num)
       call gy2dlist_to_tp_gy2dlist_allmu(tp_gy2dmu_head,gy2dmu_head,rank,mu_num)
       deallocate(gy2dmu_head)

    end subroutine tp_gy_solve_allmu



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

    subroutine tp_gy2dlist_to_gy2dlist_allmu(gy2dmu_head,tp_gy2dallmu_head,mu_num)
       class(tp_gy2dallmu_node),dimension(:), pointer, intent(in) :: tp_gy2dallmu_head
       class(tp_gy2dallmu_node),dimension(:), pointer :: tpgy2dmutmp
       class(gy2dmu_node), dimension(:), pointer,intent(inout) :: gy2dmu_head
       class(gy2dmu_node), dimension(:), pointer :: gy2dmutmp
       int4,intent(in) :: mu_num
       int4 :: i
         
       do i=1,mu_num
         if(.NOT.associated(gy2dmu_head(i)%ptr)) then
           print*, "#ERROR: gy2dmu_head is not allocated."
           stop
         end if
         if(.NOT.associated(tp_gy2dallmu_head(i)%ptr)) then
           print*, "#ERROR: tp_gy2dallmu_head is not allocated."
           stop
         end if
       end do
       
       allocate(tpgy2dmutmp(1:mu_num),gy2dmutmp(1:mu_num))
       do i=1,mu_num
 !        allocate(tpgy2dmutmp(i)%ptr,gy2dmutmp(i)%ptr)
         tpgy2dmutmp(i)%ptr=>tp_gy2dallmu_head(i)%ptr
         gy2dmutmp(i)%ptr=>gy2dmu_head(i)%ptr
       end do 

     do i=1,mu_num
       do while(associated(tpgy2dmutmp(i)%ptr))
          if(.not.associated(tpgy2dmutmp(i)%ptr%next)) then
            exit
          else
            gy2dmutmp(i)%ptr%coords(1:3)=tpgy2dmutmp(i)%ptr%coords(1:3)
            allocate(gy2dmutmp(i)%ptr%next)
            gy2dmutmp(i)%ptr=>gy2dmutmp(i)%ptr%next
            tpgy2dmutmp(i)%ptr=>tpgy2dmutmp(i)%ptr%next
          end if
       end do
     enddo     
       deallocate(gy2dmutmp)
       deallocate(tpgy2dmutmp)
    end subroutine tp_gy2dlist_to_gy2dlist_allmu



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
    end subroutine ful2dlist_to_tp_ful2dlist

  subroutine gy2dlist_to_tp_gy2dlist_allmu(tp_gy2dallmu_head,gy2dmu_head,rank,mu_num)
       class(tp_gy2dallmu_node),dimension(:), pointer, intent(inout) :: tp_gy2dallmu_head
       class(tp_gy2dallmu_node),dimension(:), pointer :: tpgy2dmutmp
       class(gy2dmu_node),dimension(:), pointer, intent(in) :: gy2dmu_head
       class(gy2dmu_node),dimension(:), pointer :: gy2dmutmp
       int4, intent(in) :: rank,mu_num
       int4 :: i

       do i=1,mu_num
       if(.NOT.associated(gy2dmu_head(i)%ptr)) then
          print*, "#ERROR: gy2dmu_head is not allocated."
          stop
       end if
       if(.NOT.associated(tp_gy2dallmu_head(i)%ptr)) then
          print*,  "rank=",rank,"#ERROR: tp_gy2dallmu_head is not allocated."
          stop
       end if
       end do

       allocate(tpgy2dmutmp(1:mu_num),gy2dmutmp(1:mu_num))
       do i=1,mu_num
!         allocate(tpgy2dmutmp(i)%ptr,gy2dmutmp(i)%ptr)
         tpgy2dmutmp(i)%ptr=>tp_gy2dallmu_head(i)%ptr
         gy2dmutmp(i)%ptr=>gy2dmu_head(i)%ptr
       end do

     do i=1,mu_num 
       do while(associated(gy2dmutmp(i)%ptr))
          if(.not.associated(gy2dmutmp(i)%ptr%next)) then
            exit
          else
            tpgy2dmutmp(i)%ptr%coords(1:3)=gy2dmutmp(i)%ptr%coords(1:3)
   !         allocate(tpful2dtmp%next)
            tpgy2dmutmp(i)%ptr=>tpgy2dmutmp(i)%ptr%next
            gy2dmutmp(i)%ptr=>gy2dmutmp(i)%ptr%next
          end if
       end do
     end do
       deallocate(gy2dmutmp)
       deallocate(tpgy2dmutmp)   
   end subroutine gy2dlist_to_tp_gy2dlist_allmu   

 end module m_tp_para_orbit
