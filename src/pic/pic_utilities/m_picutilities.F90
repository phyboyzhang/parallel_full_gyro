module m_picutilities
#include "work_precision.h"
use utilities_module, only: gp_error
use piclayout, only:& 
     parameters_set_2d, &
     ful2d_node, &
     ful2drank_node, &
     gy2drank_node, &
     ful2dsend_node, &
     gy2d_node,gy2dsend_node
    
use paradata_type, only: pic_para_total2d_base
use m_mpilayout, only: t_layout_2d, &
                 get_layout_2d_box_index
use paradata_utilities, only: compute_process_of_point_per_per, &
                              coordinates_pointoutbound_per_per, & 
                              prep_for_mpi_alltoallv_with_zeroinput
use m_parautilities, only: combine_boundfield_between_neighbor_alltoall
use orbit_data_base, only: tp_ful2d_node,tp_ful2dsend_node, &
                           tp_gy2d_node,tp_gy2dsend_node
use piclayout, only: ful2d_node,ful2dsend_node 

implicit none
include "mpif.h"

     public ::  mpi2d_alltoallv_send_particle_2d, &
                mpi2d_alltoallv_send_particle_2d_gy, &
                partition_density_to_grid_ful, &
                singlepart_to_mesh, &
                grid_global_2dcoord, &
                particle_to_2dgrid,  &
                low_localindex_of_particle, &
                tp_mpi2d_alltoallv_send_particle_2d, &
                sort_particles_among_ranks, &
                tp_sort_particles_among_ranks, &
                tp_mpi2d_alltoallv_send_particle_2d_gy, &
                tp_sort_particles_among_ranks_gy, &
                sort_particles_among_ranks_gy

    
contains

!!! allocate all the particles located in one box to its mesh nodes  
  subroutine partition_density_to_grid_ful(pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
  !  class(pic_field_2d_base), pointer :: field2d
    class(ful2d_node), pointer :: temp
    real8 :: x(2), delta(2),xmin(2)
    int4  :: boxindex(4),index(2), rank
    int4 :: row=1, n1,n2, numproc1,numproc2
    int4,dimension(:,:),pointer :: rw,re,rn,rs,rsw,rse,rne,rnw
!    real8,dimension(:,:),pointer :: density
    int4 :: ierr, i
           
    delta(1)=pic2d%para2d%m_x1%delta_eta
    delta(2)=pic2d%para2d%m_x2%delta_eta
    xmin(1)=pic2d%para2d%m_x1%eta_min
    xmin(2)=pic2d%para2d%m_x2%eta_min
    
    if(.not.associated(pic2d%ful2d_head)) then
       print*, "error: the particles are not sorted by the rank of processes"
       stop
    end if
    
    temp=>pic2d%ful2d_head
    do while(associated(temp))
      if(associated(temp%next)) then
         x(1:2)=temp%coords(1:2)
         call singlepart_to_mesh(pic2d%field2d%den,x,delta,xmin)      
         temp=>temp%next
      else
         exit
      end if
    end do

    call combine_boundfield_between_neighbor_alltoall(pic2d%field2d%den,pic2d%para2d%numproc,pic2d%layout2d)    
!    deallocate(pic2d%ful2d_head)
    nullify(pic2d%ful2d_head)
  
    return
  end subroutine partition_density_to_grid_ful
 

  !!!!!!!!! the arrangement of the trapezoid
  !!          3 ************* 4
  !!           *  c  *   d   *
  !!          *     f*        *
  !          *       * x       * 
  !         *********************     
  !        *     g   *    h      *
  !       *         e*            *  
  !      *           *             *
  !   1 ***************************** 2
  !            a              b

!!! allocate the particle to the cloest nodes of the mesh  
  subroutine singlepart_to_mesh(density,x,delta,x_min)
    real8, dimension(:,:), pointer, intent(inout) :: density
    real8, dimension(:),  intent(inout) :: x, delta, x_min
!    int4,dimension(:),intent(in) :: index
    real8 :: trap_coord(4,2)
    real8 :: ptc_ratio(4)  !!!ptc->particle
    int4 :: low(2)
    int4 :: i,j, n
    n=size(x_min)

!!$    allocate(trap_coord(2*n,n))
!!$    allocate(ptc_ratio(2*n))
    
    call grid_global_2dcoord(trap_coord,x,delta,x_min)
    call particle_to_2dgrid(ptc_ratio, trap_coord, x)
    call low_localindex_of_particle(low,x,x_min,delta)
    
    !!! here only consider two-dimensional case, n=2
    i=low(1)
    j=low(2)
    density(i+1,j)=density(i+1,j)+ ptc_ratio(3)
    density(i+1,j+1)=density(i+1,j+1)+ptc_ratio(4)
    density(i,j)=density(i,j)+ptc_ratio(1)
    density(i,j+1)=density(i,j+1)+ptc_ratio(2)

!    print*, i,j,ptc_ratio
  end subroutine singlepart_to_mesh

!!!! The lower index of the closet grid  which surrounds the particle
  subroutine low_localindex_of_particle(low,x,x_min,delta)
    real8, dimension(:), intent(in) :: x, delta,x_min !!! x_min is the starting
     !! position of the box 
    int4,  dimension(:), intent(inout) :: low  ! the index of the local low index of the grid where
                                 ! particle locates
    int4 :: n,i,j

    n=size(x)
    do i=1,n
       j=floor((x(i)-x_min(i))/delta(i))
 !      local=(x(i)-x_min(i))/delta(i)-real(i,8)
       low(i)=j+1
    end do
  end subroutine low_localindex_of_particle

  
!!!! The global coordinate of the vertices of the trapezoid which surrounds the particle  
  subroutine grid_global_2dcoord(trape_coord,x,delta,x_min)
    real8, intent(inout) :: trape_coord(4,2)
    real8, intent(in)  ::  x(2), delta(2), x_min(2)
!    int4, intent(in) :: index(2)
    int4 :: i,j
    real8 :: x1(2),x2(2)

    i=floor((x(1)-x_min(1))/delta(1))
    j=floor((x(2)-x_min(2))/delta(2))

    x1(1)=real(i,8)*delta(1)
    x1(2)=real(j,8)*delta(2)

    trape_coord(1,1)=x_min(1)+x1(1)             !real(i+1,8)*delta(1)
    trape_coord(1,2)=x_min(2)+x1(2)                !real(j,8)*delta(2)
    trape_coord(2,1)=trape_coord(1,1)             !real(i+1,8)*delta(1)
    trape_coord(2,2)=trape_coord(1,2)+delta(2)
    trape_coord(3,1)=trape_coord(1,1)+delta(1)
    trape_coord(3,2)=trape_coord(1,2)
    trape_coord(4,1)=trape_coord(1,1)+delta(1)
    trape_coord(4,2)=trape_coord(1,2)+delta(2)
   
  end subroutine grid_global_2dcoord
  
!!!! The partition of the particle to the vertices of the trapezoid which surrounds it  

  subroutine particle_to_2dgrid(ptc_ratio,trap_coord,x_gl)  !!!！ trapezoid
    real8,dimension(:), intent(inout) :: ptc_ratio
    real8, intent(in)    :: trap_coord(4,2)  ! the global coord of the four votex
    real8, intent(in)    :: x_gl(2)     ! the real global coordinate of particle
    real8 :: vol1,vol2,vol3,vol4, vol
    real8 ::  a,b,c,d,e,f,g,h  ! a>c,b>d
    

    a=abs(trap_coord(1,2)-x_gl(2))
    b=abs(trap_coord(2,2)-x_gl(2))
    c=abs(trap_coord(3,2)-x_gl(2))
    d=abs(trap_coord(4,2)-x_gl(2))
    e=abs(trap_coord(1,1)-x_gl(1))
    f=abs(trap_coord(3,1)-x_gl(1))
    g=c+(a-c)*f/(e+f)
    h=d+(b-d)*f/(e+f)

    vol1=(h+d)*f/2.0_f64
    vol2=(c+g)*f/2.0_f64
    vol3=(h+b)*e/2.0_f64
    vol4=(g+a)*e/2.0_f64
    vol=vol1+vol2+vol3+vol4
    
    ptc_ratio(1)=vol1/vol
    ptc_ratio(2)=vol2/vol
    ptc_ratio(3)=vol3/vol
    ptc_ratio(4)=vol4/vol

  end subroutine particle_to_2dgrid


   subroutine sort_particles_among_ranks(ful2d_head,ful2dsend_head, pic2d, num)
     class(pic_para_total2d_base), pointer :: pic2d
     class(ful2d_node), pointer, intent(inout) :: ful2d_head
     class(ful2dsend_node),dimension(:), pointer, intent(inout) :: ful2dsend_head
     int4, dimension(:), pointer, intent(inout) :: num
     class(ful2d_node), pointer :: tmp,tmphead
     class(ful2dsend_node),dimension(:), pointer :: outcur
     int4 :: size,prank,rank
     real8 :: coords(4)
     int4 :: i,j
 
     size=pic2d%layout2d%collective%size
     rank=pic2d%layout2d%collective%rank
     allocate(outcur(0:size-1))
     tmp=>ful2d_head
     do i=0,size-1
        outcur(i)%ptr=>ful2dsend_head(i)%ptr
     end do

     do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(prank.ne.rank) then
            num(prank)=num(prank)+1
            outcur(prank)%ptr%coords=tmp%coords
            allocate(outcur(prank)%ptr%next)
            outcur(prank)%ptr=>outcur(prank)%ptr%next

! delete this particle from the ful2d_head list
            tmp=>tmp%next
            ful2d_head=>tmp

         else 
!!!!tmphead points to pic2d%ful2d_head for the first time
            num(rank)=num(rank)+1
            tmphead=>ful2d_head           
            tmp=>tmp%next
            exit
         end if
       end if
    end do

    do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)         
         if(prank.ne.rank) then
           num(prank)=num(prank)+1
          ! delete this particle from the ful2d_head list
          !  add the deleted particle to full2drank_head list
           outcur(prank)%ptr%coords=tmp%coords
          
           allocate(outcur(prank)%ptr%next)
           outcur(prank)%ptr=>outcur(prank)%ptr%next
           if(.not.associated(tmp%next%next)) then
              tmphead%next=>tmp%next
              tmphead=>tmphead%next
              exit
           else 
              tmp=>tmp%next
           end if
         else
           num(rank)=num(rank)+1
           tmphead%next=>tmp
           tmphead=>tmphead%next
           tmp=>tmp%next
         end if

       endif
    end do

     nullify(tmp)
     nullify(tmphead)
     nullify(outcur)
   end subroutine sort_particles_among_ranks

   subroutine sort_particles_among_ranks_gy(gy2d_head,gy2dsend_head, pic2d, num)
     class(pic_para_total2d_base), pointer :: pic2d
     class(gy2d_node), pointer, intent(inout) :: gy2d_head
     class(gy2dsend_node),dimension(:), pointer, intent(inout) :: gy2dsend_head
     int4, dimension(:), pointer, intent(inout) :: num
     class(gy2d_node), pointer :: tmp,tmphead
     class(gy2dsend_node),dimension(:), pointer :: outcur
     int4 :: size,prank,rank
     real8 :: coords(4)
     int4 :: i,j
 
     size=pic2d%layout2d%collective%size
     rank=pic2d%layout2d%collective%rank
     allocate(outcur(0:size-1))
     tmp=>gy2d_head
     do i=0,size-1
        outcur(i)%ptr=>gy2dsend_head(i)%ptr
     end do

     do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(prank.ne.rank) then
            num(prank)=num(prank)+1
            outcur(prank)%ptr%coords=tmp%coords
            allocate(outcur(prank)%ptr%next)
            outcur(prank)%ptr=>outcur(prank)%ptr%next

! delete this particle from the gy2d_head list
            tmp=>tmp%next
            gy2d_head=>tmp

         else 
!!!!tmphead points to pic2d%gy2d_head for the first time
            num(rank)=num(rank)+1
            tmphead=>gy2d_head           
            tmp=>tmp%next
            exit
         end if
       end if
    end do

    do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)         
         if(prank.ne.rank) then
           num(prank)=num(prank)+1
          ! delete this particle from the gy2d_head list
          !  add the deleted particle to gy2drank_head list
           outcur(prank)%ptr%coords=tmp%coords
          
           allocate(outcur(prank)%ptr%next)
           outcur(prank)%ptr=>outcur(prank)%ptr%next
           if(.not.associated(tmp%next%next)) then
              tmphead%next=>tmp%next
              tmphead=>tmphead%next
              exit
           else 
              tmp=>tmp%next
           end if
         else
           num(rank)=num(rank)+1
           tmphead%next=>tmp
           tmphead=>tmphead%next
           tmp=>tmp%next
         end if

       endif
    end do

     nullify(tmp)
     nullify(tmphead)
     nullify(outcur)
   end subroutine sort_particles_among_ranks_gy




   subroutine tp_sort_particles_among_ranks(tp_ful2d_head,tp_ful2dsend_head, pic2d, num)
     class(pic_para_total2d_base), pointer :: pic2d
     class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
     class(tp_ful2dsend_node),dimension(:), pointer, intent(inout) :: tp_ful2dsend_head
     int4, dimension(:), pointer, intent(inout) :: num
     class(tp_ful2d_node), pointer :: tmp,tmphead
     class(tp_ful2dsend_node),dimension(:), pointer :: outcur
     int4 :: size,rank,prank
     real8 :: coords(4)
     int4 :: i,j
 
     rank=pic2d%layout2d%collective%rank
     size=pic2d%layout2d%collective%size
     allocate(outcur(0:size-1))
     tmp=>tp_ful2d_head
     do i=0,size-1
        outcur(i)%ptr=>tp_ful2dsend_head(i)%ptr
     end do

     do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(prank.ne.rank) then
            num(prank)=num(prank)+1
            outcur(prank)%ptr%coords=tmp%coords
            outcur(prank)%ptr%tp=tmp%tp
            allocate(outcur(prank)%ptr%next)
            outcur(prank)%ptr=>outcur(prank)%ptr%next

! delete this particle from the ful2d_head list
            tmp=>tmp%next
            tp_ful2d_head=>tmp

         else 
!!!!tmphead points to pic2d%ful2d_head for the first time
            num(rank)=num(rank)+1
            tmphead=>tp_ful2d_head           
            tmp=>tmp%next
            exit
         end if
       end if
    end do

    do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)         
         if(prank.ne.rank) then
           num(prank)=num(prank)+1
          ! delete this particle from the ful2d_head list
          !  add the deleted particle to full2drank_head list
           outcur(prank)%ptr%coords=tmp%coords
           outcur(prank)%ptr%tp=tmp%tp 
           allocate(outcur(prank)%ptr%next)
           outcur(prank)%ptr=>outcur(prank)%ptr%next
!           tmp=>tmp%next
!          !!! The following "if" is critical
           if(.not.associated(tmp%next%next)) then
              tmphead%next=>tmp%next
              tmphead=>tmphead%next
              exit
           else 
              tmp=>tmp%next
           end if
         else
           num(rank)=num(rank)+1
           tmphead%next=>tmp
           tmphead=>tmphead%next
           tmp=>tmp%next
         end if

       endif
    end do
    
     nullify(outcur)
     nullify(tmp)
     nullify(tmphead)
   end subroutine tp_sort_particles_among_ranks

   subroutine tp_sort_particles_among_ranks_gy(tp_gy2d_head,tp_gy2dsend_head, pic2d, num)
     class(pic_para_total2d_base), pointer :: pic2d
     class(tp_gy2d_node), pointer, intent(inout) :: tp_gy2d_head
     class(tp_gy2dsend_node),dimension(:), pointer, intent(inout) :: tp_gy2dsend_head
     int4, dimension(:), pointer, intent(inout) :: num
     class(tp_gy2d_node), pointer :: tmp,tmphead
     class(tp_gy2dsend_node),dimension(:), pointer :: outcur
     int4 :: size,rank,prank
     real8 :: coords(4)
     int4 :: i,j
 
     rank=pic2d%layout2d%collective%rank
     size=pic2d%layout2d%collective%size
     allocate(outcur(0:size-1))
     tmp=>tp_gy2d_head
     do i=0,size-1
        outcur(i)%ptr=>tp_gy2dsend_head(i)%ptr
     end do

     do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(prank.ne.rank) then
            num(prank)=num(prank)+1
            outcur(prank)%ptr%coords=tmp%coords
            outcur(prank)%ptr%tp=tmp%tp
            allocate(outcur(prank)%ptr%next)
            outcur(prank)%ptr=>outcur(prank)%ptr%next

! delete this particle from the ful2d_head list
            tmp=>tmp%next
            tp_gy2d_head=>tmp

         else 
!!!!tmphead points to pic2d%ful2d_head for the first time
            num(rank)=num(rank)+1
            tmphead=>tp_gy2d_head           
            tmp=>tmp%next
            exit
         end if
       end if
    end do

    do while(associated(tmp))
       if(.not.associated(tmp%next)) then
         exit
       else
         call coordinates_pointoutbound_per_per(tmp%coords(1:2),pic2d%para2d%gxmin,pic2d%para2d%gxmax)
         prank=compute_process_of_point_per_per(tmp%coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)         
         if(prank.ne.rank) then
           num(prank)=num(prank)+1
          ! delete this particle from the ful2d_head list
          !  add the deleted particle to full2drank_head list
           outcur(prank)%ptr%coords=tmp%coords
           outcur(prank)%ptr%tp=tmp%tp 
           allocate(outcur(prank)%ptr%next)
           outcur(prank)%ptr=>outcur(prank)%ptr%next
!           tmp=>tmp%next
!          !!! The following "if" is critical
           if(.not.associated(tmp%next%next)) then
              tmphead%next=>tmp%next
              tmphead=>tmphead%next
              exit
           else 
              tmp=>tmp%next
           end if
         else
           num(rank)=num(rank)+1
           tmphead%next=>tmp
           tmphead=>tmphead%next
           tmp=>tmp%next
         end if

       endif
    end do
    
     nullify(outcur)
     nullify(tmp)
     nullify(tmphead)
   end subroutine tp_sort_particles_among_ranks_gy



 !!!! send particles which locate out of current box to those boxes where they locate
 !!!  and store the particles on the ful2dsend_head list 
 subroutine mpi2d_alltoallv_send_particle_2d_one(ful2d_head,ful2dsend_head,num,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2dsend_node), dimension(:), pointer, intent(in) :: ful2dsend_head
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    int4, dimension(:), pointer, intent(in) :: num
!    int4, intent(in) :: rank,size,comm
    int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls
    real8, dimension(:), pointer :: sbuf,rbuf
    int4 :: i,j,h,numout,ierr, rank,size,comm,numsend
    class(ful2dsend_node),dimension(:), pointer :: rankcur
    class(ful2d_node),     pointer :: coordcur

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    allocate(rankcur(0:size-1))
    do i=0,size-1
       rankcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do    
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))
    call mpi_alltoall(num,1,mpi_integer,rcounts,1,mpi_integer,comm,ierr)
   numsend=0 
   do i=0,size-1
      numsend=numsend+num(i)
   end do
   if(numsend==0) then
     allocate(sbuf(0:0))
   else
     allocate(sbuf(0:numsend*4-1))
   endif

   numout=0
   do i=0, size-1
      numout=numout+rcounts(i)
   end do
   allocate(rbuf(0:numout*4-1))

   do i=0,size-1
      scounts(i)=num(i)*4
      rcounts(i)=rcounts(i)*4
       if(i==0) then
          sdispls(i)=0
          rdispls(i)=0
       else
          sdispls(i)=sdispls(i-1)+scounts(i-1)
          rdispls(i)=rdispls(i-1)+rcounts(i-1)
       end if
    end do    
!print*, "scounts=",scounts      
    h=0
    do i=0, size-1
       do while(associated(rankcur(i)%ptr))  !  .or.associated(rankcur(i)%ptr%next))
          if(associated(rankcur(i)%ptr%next)) then
            sbuf(4*h:4*h+3)=rankcur(i)%ptr%coords(1:4)
            rankcur(i)%ptr=>rankcur(i)%ptr%next
            h=h+1
          else 
            exit
          end if
       end do
!       nullify(pic2d%ful2dsend_head(i)%ptr)
    end do

     call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)

!!! insert the particles into the tail of the head linked list 
   coordcur=>ful2d_head
   do while(associated(coordcur))
      if(associated(coordcur%next)) then
         coordcur=>coordcur%next
      else 
         exit    
      end if
   end do         

   do i=0,numout-1
          coordcur%coords(1:4)=rbuf(4*i:4*i+3)
          allocate(coordcur%next)
          coordcur=>coordcur%next
   end do

!print*, "rank=",pic2d%layout2d%collective%rank       
  end subroutine mpi2d_alltoallv_send_particle_2d_one


  subroutine mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2dsend_node), dimension(:), pointer, intent(in) :: ful2dsend_head
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    int4, dimension(:), pointer, intent(in) :: num
!    int4, intent(inout) :: numleft
!    int4, intent(in) :: rank,size,comm
    int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls
    real8, dimension(:), pointer :: sbuf,rbuf
    int4 :: i,j,h,numout,ierr, rank,size,comm,numsend
    class(ful2dsend_node),dimension(:), pointer :: rankcur
    class(ful2d_node),     pointer :: coordcur

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    allocate(rankcur(0:size-1))
    do i=0,size-1
       rankcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do    
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))
    call mpi_alltoall(num,1,mpi_integer,rcounts,1,mpi_integer,comm,ierr)
   numsend=0 
   do i=0,size-1
      numsend=numsend+num(i)
   end do
   numsend=numsend-num(rank)
   numout=0
   do i=0, size-1
      numout=numout+rcounts(i)
   end do
!   numleft=numout
   numout=numout-rcounts(rank)
 
   if(numsend==0) then
     allocate(sbuf(0:0))
     do i=0,size-1
        scounts(i)=0
        sdispls(i)=0
     end do
   else
     allocate(sbuf(0:numsend*4-1))
     do i=0,size-1
        if(i==rank) then
          scounts(i)=0
        else
          scounts(i)=num(i)*4
        end if
        if(i==0) then
          sdispls(i)=0
        else
          sdispls(i)=sdispls(i-1)+scounts(i-1)
        end if
     end do
   endif

   if(numout==0) then
     allocate(rbuf(0:0))
     do i=0,size-1
        rcounts(i)=0
        rdispls(i)=0
     end do
   else
     allocate(rbuf(0:numout*4-1))
     do i=0,size-1
       if(i==rank) then
         rcounts(i)=0
       else
         rcounts(i)=rcounts(i)*4
       end if

       if(i==0) then
         rdispls(i)=0
       else
         rdispls(i)=rdispls(i-1)+rcounts(i-1)
       end if
     end do
  end if

  if(numsend==0) then
       call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)    
  else
    h=0
    do i=0, size-1
       if(size==rank) then
          goto 60
       else
         do while(associated(rankcur(i)%ptr))  !  .or.associated(rankcur(i)%ptr%next))
            if(associated(rankcur(i)%ptr%next)) then
               sbuf(4*h:4*h+3)=rankcur(i)%ptr%coords(1:4)
               rankcur(i)%ptr=>rankcur(i)%ptr%next
               h=h+1
            else 
              exit
            end if
         end do
       end if
60  end do

     call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
   end if

!!! insert the particles into the tail of the head linked list 
   if(.not.associated(ful2d_head)) then
       print*, "ful2d_head is not allocated"
       stop
   end if
   coordcur=>ful2d_head
   do while(associated(coordcur))
      if(associated(coordcur%next)) then
         coordcur=>coordcur%next
      else 
         exit    
      end if
   end do         

   if(numout==0) then
     goto 20
   else
     do i=0,numout-1
        coordcur%coords(1:4)=rbuf(4*i:4*i+3)
        allocate(coordcur%next)
        coordcur=>coordcur%next
     end do
20   end if

   deallocate(rankcur)
   nullify(coordcur)

 end subroutine mpi2d_alltoallv_send_particle_2d

 subroutine mpi2d_alltoallv_send_particle_2d_gy(gy2d_head,gy2dsend_head,num,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(gy2dsend_node), dimension(:), pointer, intent(in) :: gy2dsend_head
    class(gy2d_node), pointer, intent(inout) :: gy2d_head
    int4, dimension(:), pointer, intent(in) :: num
    int4 :: numgr
!    int4, intent(in) :: rank,size,comm
    int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls,rbuf0
    real8, dimension(:), pointer :: sbuf,rbuf
    int4 :: i,j,h,numout,ierr, rank,size,comm,numsend
    class(gy2dsend_node),dimension(:), pointer :: rankcur
    class(gy2d_node),     pointer :: coordcur

    numgr=3
    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    allocate(rankcur(0:size-1))
    do i=0,size-1
       rankcur(i)%ptr=>gy2dsend_head(i)%ptr
    end do    
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1),rbuf0(0:size-1))
    call mpi_alltoall(num,1,mpi_integer,rbuf0,1,mpi_integer,comm,ierr)
    call prep_for_mpi_alltoallv_with_zeroinput(rank,size,numgr,num,rbuf0,scounts, &
             rcounts,sdispls,rdispls,numsend,numout)
    if(numsend==0) then
       allocate(sbuf(0:0))
    else
       allocate(sbuf(0:numgr*numsend-1))
    end if
    if(numout==0) then
       allocate(rbuf(0:0))
    else
       allocate(rbuf(0:numgr*numout-1))
    end if
    

  if(numsend==0) then
       call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)    
  else
    h=0
    do i=0, size-1
       if(size==rank) then
          goto 60
       else
         do while(associated(rankcur(i)%ptr))  !  .or.associated(rankcur(i)%ptr%next))
            if(associated(rankcur(i)%ptr%next)) then
               sbuf(numgr*h:numgr*h+2)=rankcur(i)%ptr%coords(1:3)
               rankcur(i)%ptr=>rankcur(i)%ptr%next
               h=h+1
            else 
              exit
            end if
         end do
       end if
60  end do

     call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
   end if

!!! insert the particles into the tail of the head linked list 
   if(.not.associated(gy2d_head)) then
       print*, "gy2d_head is not allocated"
       stop
   end if
   coordcur=>gy2d_head
   do while(associated(coordcur))
      if(associated(coordcur%next)) then
         coordcur=>coordcur%next
      else 
         exit    
      end if
   end do         

   if(numout==0) then
     goto 20
   else
     do i=0,numout-1
        coordcur%coords(1:3)=rbuf(4*i:4*i+2)
        allocate(coordcur%next)
        coordcur=>coordcur%next
     end do
20   end if

   deallocate(rankcur)
   nullify(coordcur)

 end subroutine mpi2d_alltoallv_send_particle_2d_gy


 

  subroutine tp_mpi2d_alltoallv_send_particle_2d(tp_ful2d_head,numleft, &
             tp_ful2dsend_head,num,numgr,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(tp_ful2dsend_node), dimension(:),pointer, intent(in) :: tp_ful2dsend_head
    class(tp_ful2d_node), pointer, intent(inout) :: tp_ful2d_head
    int4, dimension(:), pointer, intent(in) :: num
    int4, intent(inout) :: numleft
    int4, intent(in) :: numgr
    int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls
    real8, dimension(:), pointer :: sbuf,rbuf
    class(tp_ful2dsend_node),dimension(:), pointer :: rankcur
    class(tp_ful2d_node),     pointer :: coordcur
    int4 :: i,j,h,ierr,rank,size,comm, nums,numre
 
    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    allocate(rankcur(0:size-1))
    do i=0,size-1
       rankcur(i)%ptr=>tp_ful2dsend_head(i)%ptr
    end do    
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))
   call mpi_alltoall(num,1,mpi_integer,rcounts,1,mpi_integer,comm,ierr) 
   nums=0 
   do i=0,size-1
      nums=nums+num(i)
   end do
   nums=nums-num(rank)

!   allocate(sbuf(0:numout*numgr-1))
   numre=0
   do i=0, size-1
      numre=numre+rcounts(i)
   end do

   numleft=numre  !!!! numleft is use to show whether there exists particles in
                  !!!! the this process.
   numre=numre-rcounts(rank)
!   allocate(rbuf(0:numre*numgr-1))
  
   if(nums==0) then
     allocate(sbuf(0:0))
     do i=0,size-1
        scounts(i)=0
        sdispls(i)=0
     end do
   else
     allocate(sbuf(0:nums*numgr-1))
     do i=0,size-1
        if(i==rank) then
           scounts(i)=0
        else
           scounts(i)=num(i)*numgr
        end if
        if(i==0) then
           sdispls(i)=0
        else
           sdispls(i)=sdispls(i-1)+scounts(i-1) 
        endif
     enddo
    endif 

   if(numre==0) then
     allocate(rbuf(0:0))
     do i=0,size-1
        rcounts(i)=0
        rdispls(i)=0
     end do
   else
     allocate(rbuf(0:numre*numgr-1))
     do i=0,size-1
        if(i==rank) then
           rcounts(i)=0
        else
           rcounts(i)=rcounts(i)*numgr
        end if
        if(i==0) then
           rdispls(i)=0
        else
           rdispls(i)=rdispls(i-1)+rcounts(i-1) 
        endif
     enddo
    endif 
!print*, rank

    if(nums==0) then
       call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)

!print*, rank
    else
      h=0
      do i=0, size-1
        if(size==rank) then
           goto 40
        else 
          do while(associated(rankcur(i)%ptr))  !  .or.associated(rankcur(i)%ptr%next))
            if(associated(rankcur(i)%ptr%next)) then
              sbuf(numgr*h:numgr*h+numgr-2)=rankcur(i)%ptr%coords(1:numgr-1)
              sbuf(numgr*h+numgr-1)=real(rankcur(i)%ptr%tp,8)
              rankcur(i)%ptr=>rankcur(i)%ptr%next
              h=h+1
            else 
              exit
            end if
          end do
40        end if
    end do
      call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
    end if

!!! insert the particle into the head linked list 
   if(.not.associated(tp_ful2d_head)) then
      print*, "pic2d%ful2d_head is not allocated"
      stop
   end if
   coordcur=>tp_ful2d_head
   do while(associated(coordcur))
      if(associated(coordcur%next)) then
         coordcur=>coordcur%next
      else 
         exit 
      end if
   end do

   if(numre==0) then
     goto 10
   else
     do i=0,numre-1
          coordcur%coords(1:numgr-1)=rbuf(numgr*i:numgr*i+numgr-2)
          coordcur%tp = NINT(rbuf(numgr*i+numgr-1))
          allocate(coordcur%next)
          coordcur=>coordcur%next
     end do
10   end if

   deallocate(rankcur)
   nullify(coordcur)
   end subroutine tp_mpi2d_alltoallv_send_particle_2d

 subroutine tp_mpi2d_alltoallv_send_particle_2d_gy(tp_gy2d_head,numleft, &
             tp_gy2dsend_head,num,numgr,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(tp_gy2dsend_node), dimension(:),pointer, intent(in) :: tp_gy2dsend_head
    class(tp_gy2d_node), pointer, intent(inout) :: tp_gy2d_head
    int4, dimension(:), pointer, intent(in) :: num
    int4, intent(inout) :: numleft
    int4, intent(in) :: numgr
    int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls
    real8, dimension(:), pointer :: sbuf,rbuf
    class(tp_gy2dsend_node),dimension(:), pointer :: rankcur
    class(tp_gy2d_node),     pointer :: coordcur
    int4 :: i,j,h,ierr,rank,size,comm, nums,numre
 
    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    allocate(rankcur(0:size-1))
    do i=0,size-1
       rankcur(i)%ptr=>tp_gy2dsend_head(i)%ptr
    end do    
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))
   call mpi_alltoall(num,1,mpi_integer,rcounts,1,mpi_integer,comm,ierr) 
   nums=0 
   do i=0,size-1
      nums=nums+num(i)
   end do
   nums=nums-num(rank)

!   allocate(sbuf(0:numout*numgr-1))
   numre=0
   do i=0, size-1
      numre=numre+rcounts(i)
   end do

   numleft=numre  !!!! numleft is use to show whether there exists particles in
                  !!!! the this process.
   numre=numre-rcounts(rank)
!   allocate(rbuf(0:numre*numgr-1))
  
   if(nums==0) then
     allocate(sbuf(0:0))
     do i=0,size-1
        scounts(i)=0
        sdispls(i)=0
     end do
   else
     allocate(sbuf(0:nums*numgr-1))
     do i=0,size-1
        if(i==rank) then
           scounts(i)=0
        else
           scounts(i)=num(i)*numgr
        end if
        if(i==0) then
           sdispls(i)=0
        else
           sdispls(i)=sdispls(i-1)+scounts(i-1) 
        endif
     enddo
    endif 

   if(numre==0) then
     allocate(rbuf(0:0))
     do i=0,size-1
        rcounts(i)=0
        rdispls(i)=0
     end do
   else
     allocate(rbuf(0:numre*numgr-1))
     do i=0,size-1
        if(i==rank) then
           rcounts(i)=0
        else
           rcounts(i)=rcounts(i)*numgr
        end if
        if(i==0) then
           rdispls(i)=0
        else
           rdispls(i)=rdispls(i-1)+rcounts(i-1) 
        endif
     enddo
    endif 
!print*, rank

    if(nums==0) then
       call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)

!print*, rank
    else
      h=0
      do i=0, size-1
        if(size==rank) then
           goto 40
        else 
          do while(associated(rankcur(i)%ptr))  !  .or.associated(rankcur(i)%ptr%next))
            if(associated(rankcur(i)%ptr%next)) then
              sbuf(numgr*h:numgr*h+numgr-2)=rankcur(i)%ptr%coords(1:numgr-1)
              sbuf(numgr*h+numgr-1)=real(rankcur(i)%ptr%tp,8)
              rankcur(i)%ptr=>rankcur(i)%ptr%next
              h=h+1
            else 
              exit
            end if
          end do
40        end if
    end do
      call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
    end if

!!! insert the particle into the head linked list 
   if(.not.associated(tp_gy2d_head)) then
      print*, "tp_gy2d_head is not allocated"
      stop
   end if
   coordcur=>tp_gy2d_head
   do while(associated(coordcur))
      if(associated(coordcur%next)) then
         coordcur=>coordcur%next
      else 
         exit 
      end if
   end do

   if(numre==0) then
     goto 10
   else
     do i=0,numre-1
          coordcur%coords(1:numgr-1)=rbuf(numgr*i:numgr*i+numgr-2)
          coordcur%tp = NINT(rbuf(numgr*i+numgr-1))
          allocate(coordcur%next)
          coordcur=>coordcur%next
     end do
10   end if

   deallocate(rankcur)
   nullify(coordcur)
   end subroutine tp_mpi2d_alltoallv_send_particle_2d_gy



end module m_picutilities
