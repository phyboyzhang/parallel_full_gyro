module m_moveparticles
#include "work_precision.h"
  use omp_lib
  use m_mpilayout, only: &
       t_layout_2d
  use piclayout,only: &
       ful2drank_node, &
       gy2drank_node, &
!       pic_para_ful2d_base, &
       ful2d_node, &
       gy2d_node, &
       ful2dsend_node, &
       gy2dsend_node

  use paradata_type, only: pic_para_total2d_base
  use m_para_spline, only: &
       para_spl2d_firstorder_derivatve_point_per_per, &
       para_compute_spl2d_field_point_per_per

  use m_parautilities,only: mpi2d_alltoallv_box_per_per
  use m_picutilities, only: mpi2d_alltoallv_send_particle_2d
  use paradata_utilities, only: compute_process_of_point_per_per

  use m_para_orbit, only: borissolve, fulrk4solve,boris_single_per_per
  implicit none

  public :: new_list_ful2d_per_per, &
            push_particle_among_box_ful2d_per_per, &
            new_position_per_per_ful 
  
contains



  subroutine new_list_ful2d_per_per(ful2d_head,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), pointer, intent(inout) :: ful2d_head 
    int4,dimension(:), pointer :: num
    class(ful2dsend_node), dimension(:), pointer :: ful2dsend_head
    int4 :: size
    int4 ::  i

    size=pic2d%layout2d%collective%size
!    rank=pic2d%layout2d%collective%rank
!    comm=pic2d%layout2d%collective%comm
    allocate(num(0:pic2d%layout2d%collective%size-1))
    allocate(ful2dsend_head(0:size-1))
    do i=0,size-1
      allocate(ful2dsend_head(i)%ptr)
    end do
    call push_particle_among_box_ful2d_per_per(ful2dsend_head,pic2d,num)
    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head, &
         num,pic2d)

    deallocate(ful2dsend_head,num)
    
  end subroutine new_list_ful2d_per_per

 

 !! ????? will do: allocate the ful2dsr_head list, which will also be nullified for each step.
    !!===> get the new position of ful2d_list among the current box, and keep the particle still locating in current box in the list, 
!!    while move the particles whose new position locates out of the current box to the full2drank_list.
  subroutine push_particle_among_box_ful2d_per_per(ful2dsend_head,pic2d,num)
    class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
    class(ful2dsend_node), dimension(:), intent(inout) :: ful2dsend_head
  !  class(ful2drank_node),dimension(0:pic2d%para2d%numproc1*pic2d%para2d%numproc2-1), pointer :: outcur
    class(ful2dsend_node),dimension(:), allocatable :: outcur
    class(ful2d_node), pointer :: coordcur, coordcurhead
    int4, dimension(:),pointer, intent(inout) :: num     ! num(0:pic2d%layout2d%collective%size-1)
    int4 :: prank  !! prank is the rank of the processor where particle locates
    int4 :: comm,row, numproc(2),boxindex(4),rank,size
    int4 :: numout=0
    real8 :: x(0:1),v(0:1)
    int4 :: i,h
    
    size=pic2d%layout2d%collective%size
    allocate(outcur(0:pic2d%para2d%numproc(1)*pic2d%para2d%numproc(2)-1))
    coordcur=>pic2d%ful2d_head
    comm=pic2d%layout2d%collective%comm
    rank=pic2d%layout2d%collective%rank
    row=pic2d%para2d%row
    num=0
    numproc(1)=pic2d%para2d%numproc(1)
    numproc(2)=pic2d%para2d%numproc(2)
  
!    if(.not.associated(ful2dsend_head(0)%ptr)) then
!       allocate(ful2dsend_head(0:size-1))
!       do i=0,size-1
!          allocate(pic2d%ful2dsend_head(i)%ptr)
!       end do
!    end if
    do i=0,size-1
       outcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do

    !!The first "do while" is to sort out of particles which are at the beginning of the list but out of the current box, until the appearance of the particle located in the current box.

!$omp parallel num_threads(pic2d%layout2d%collective%thread_level_provided) firstprivate(coordcur) 
    !$omp single 
    do while(associated(coordcur).and.associated(coordcur%next))
             ! new_position_per_per is to obtain the new position of the particle
 !   !$omp task firstprivate(coordcur) private(prank)
       call new_position_per_per_ful(coordcur%coords(1:2),coordcur%coords(3:4),pic2d)       
       prank=compute_process_of_point_per_per(coordcur%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
             pic2d%para2d%gxmax,pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
!print*, "coordcur%coords=",coordcur%coords
 !   !$omp end task

       if(prank.ne.rank) then
          num(prank)=num(prank)+1
         !  add the deleted particle to the particle_rankth node of full2drank_head list
          outcur(prank)%ptr%coords=coordcur%coords
          outcur(prank)%ptr%prank=prank
          allocate(outcur(prank)%ptr%next)
          outcur(prank)%ptr=>outcur(prank)%ptr%next
!if(rank==0) then
!print*, "outcur(prank)=",outcur(prank)%ptr%coords
!end if

       ! delete this particle from the ful2d_head list
          coordcur=>coordcur%next
          pic2d%ful2d_head=>coordcur
          
       else
          !!!!coordcurhead points to pic2d%ful2d_head for the first time
          coordcurhead=>pic2d%ful2d_head          
          coordcur=>coordcur%next   ! Then, coordcurh points before coordcur
          exit
       end if
    
    end do

     do while(associated(coordcur).and.associated(coordcur%next))
       !$omp task firstprivate(coordcur)   
       call new_position_per_per_ful(coordcur%coords(1:2),coordcur%coords(3:4),pic2d)
       prank=compute_process_of_point_per_per(coordcur%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
             pic2d%para2d%gxmax,pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
      !$omp end task
       if(prank.ne.rank) then
          num(prank)=num(prank)+1         
          ! delete this particle from the ful2d_head list
          !  add the deleted particle to full2drank_head list
          outcur(prank)%ptr%coords=coordcur%coords
          outcur(prank)%ptr%prank=prank

!if(rank==1) then
!print*, "outcur(prank)=",outcur(prank)%ptr%coords
!end if
          allocate(outcur(prank)%ptr%next)
          outcur(prank)%ptr=>outcur(prank)%ptr%next
          coordcur=>coordcur%next
       else
          coordcurhead%next=>coordcur
          coordcurhead=>coordcurhead%next
          coordcur=>coordcur%next      
       end if

    end do
  !$omp taskwait
  !$omp end single
!$omp end parallel

 
    deallocate(coordcur,coordcurhead)   
!    nullify(coordcur)
!    nullify(coordcurhead)
    do i=0, size-1
        deallocate(outcur(i)%ptr)
        nullify(outcur(i)%ptr)
    end do

  end subroutine push_particle_among_box_ful2d_per_per

  
  subroutine new_position_per_per_ful(x,v,pic2d)
    real8, dimension(1:2), intent(inout) :: x, v
    class(pic_para_total2d_base), pointer, intent(in) :: pic2d
    real8 ::  magf(1:3), ep_deri(1:3)
    real8 ::  x1(1:3), v1(1:3)
    real8 ::  eta_min(2),eta_max(2),delta_eta(2),dtful
    int4  ::  row,rank,NC(2),numproc(2), comm

    row=pic2d%para2d%row
    rank=pic2d%layout2d%collective%rank
    numproc(1)=pic2d%para2d%numproc(1)
    numproc(2)=pic2d%para2d%numproc(2)
    comm=pic2d%layout2d%collective%comm
    dtful=pic2d%para2d%dtful
    
    x1(1:2)=x(1:2)
    v1(1:2)=v(1:2)
    magf(1:2)=0.0
!    magf(3)=para_compute_spl2d_field_point_per_per(x,pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
!            pic2d%field2d%bf_weight_z(1,:,:),pic2d%field2d%bf_w(1,:,:),pic2d%field2d%bf_e(1,:,:), &
!            pic2d%field2d%bf_n(1,:,:), pic2d%field2d%bf_s(1,:,:), pic2d%field2d%bf_sw(1,:,:),  &
!            pic2d%field2d%bf_se(1,:,:),pic2d%field2d%bf_ne(1,:,:),pic2d%field2d%bf_nw(1,:,:))
    
!    magf(2)=para_compute_spl2d_field_point_per_per(eta_min,eta_max,row,rank,NC,numproc,x, &
!            pic2d%field2d%bf_weight_z(2,:,:),pic2d%field2d%bf_w(2,:,:),pic2d%field2d%bf_e(2,:,:), &
!            pic2d%field2d%bf_n(2,:,:), pic2d%field2d%bf_s(2,:,:),pic2d%field2d%bf_sw(2,:,:), &
!            pic2d%field2d%bf_se(2,:,:),pic2d%field2d%bf_ne(2,:,:),pic2d%field2d%bf_nw(2,:,:),comm)
    
    call para_spl2d_firstorder_derivatve_point_per_per(ep_deri(1:2),x, &
         pic2d%para2d%m_x1,pic2d%para2d%m_x2,row,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w, &
         pic2d%field2d%epwg_e, pic2d%field2d%epwg_n, pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw, &
         pic2d%field2d%epwg_se,pic2d%field2d%epwg_ne, pic2d%field2d%epwg_nw)

!    dtful=pic2d%para2d%dt_gc/real(pic2d%para2d%num_time,8)
    call boris_single(x1,v1,dtful,magf,-ep_deri)

    x(1:2)=x1(1:2)
    v(1:2)=v1(1:2)
    
    return
  end subroutine new_position_per_per_ful

  subroutine advance_particles_among_ranks_ful(ful2d_head,pic2d,rk4order,pushalgorithm)
    class(pic_para_total2d_base), pointer :: pic2d 
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
 !   class(gy2d_node),  pointer, intent(inout) :: gy2d_head
    int4,intent(in) :: rk4order
    character(len=*),intent(in) :: pushalgorithm
        select case(pushalgorithm)
          case ("boris")
             call borissolve(ful2d_head,pic2d)
          case ("rk4")
             call fulrk4solve(ful2d_head,pic2d,rk4order)
          case default
            stop
        end select

  end subroutine

!  subroutine advance_particles_among_ranks_gy(gy2d_head,pic2d,iter_num)
!    class(pic_para_total2d_base), pointer :: pic2d 
!    class(gy2d_node), pointer, intent(inout) :: gy2d_head
!    int4,intent(in) :: iter_num
!
!    call gyrork4solve(gy2d_head,pic2d,iter_num)
!
!  end subroutine


  
end module m_moveparticles
