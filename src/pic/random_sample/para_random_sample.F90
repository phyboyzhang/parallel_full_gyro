module para_random_sample
#include "work_precision.h"
  use omp_lib
  use utilities_module, only: gp_error
  use piclayout, only: &
       ful2dsend_node, &
!       pic_para_total2d_base, &
       ful2d_node
  use paradata_type, only: pic_para_total2d_base

  use constants, only: pi_
  use m_moveparticles, only: mpi2d_alltoallv_box_per_per 
  use paradata_utilities, only: compute_process_of_point_per_per
  use m_picutilities, only: mpi2d_alltoallv_send_particle_2d
  use paradata_utilities, only: coordinates_pointoutbound_per_per
!  use orbit_data_base, only: ful2d_node,ful2dsend_node
  implicit none

  public:: para_accept_reject_gaussian1d_ful2d_per_per, &
           para_accprej_gaus2d2v_fulgyro_unifield_per_per, &
           para_accprej_gaus1d2v_fulgyro_unifield_per_per
  
contains

  subroutine para_accept_reject_gaussian1d_ful2d_per_per(ful2d_head,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer ::currk, ful2dsend_head
    real8 :: mean,y1,y2,py,pymax,x,x1,py1,gmin(2),gmax(2)
    int4,dimension(:),pointer :: num
    class(ful2d_node), pointer :: ful2dtmp
    real8 :: theta, sigma, coords(2),rho,tempt
    int4 :: rank1,numcircle,size,rank,comm
    int4 :: ierr,i,j

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    rho=sqrt(2._f64*pic2d%para2d%mu)
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    allocate(num(0:size-1))
    allocate(currk(0:size-1),ful2dsend_head(0:size-1))
    do i=0,size-1
       currk(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax

    if(.not.associated(pic2d%ful2d_head)) then
      print*, "pic2d%ful2d_head is not allocated"
      stop
    end if

!    call random_seed()       
    do i=1,pic2d%para2d%numparticle
       pymax=1.0/sqrt(2.0*pi_)/sigma
       x=1.0
       py=0.0
       mean=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64 
       do while(x.gt.py) 
         y1=gmin(1)+(gmax(1)-gmin(1))*rand()    !generate_random_number()
         py=exp(-(y1-mean)*(y1-mean)/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x=pymax*rand()
       end do

       x1=1.0
       py1=0.0
       mean=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      
       do while(x1.gt.py1)
         y2=gmin(2)+(gmax(2)-gmin(2))*rand()        
         py1=exp(-(y2-mean)*(y2-mean)/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x1=pymax*rand()
       end do

     do j=0,numcircle-1
           theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=y1+rho*cos(theta)
           coords(2)=y2+rho*sin(theta)       
!           call coordinates_pointoutbound_per_per(coords,pic2d) 
           rank1=compute_process_of_point_per_per(coords,pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(rank1==rank) then
           ful2dtmp%coords(1:2)=coords(1:2)
           allocate(ful2dtmp%next)
           ful2dtmp=>ful2dtmp%next
           num(rank)=num(rank)+1
         else
           currk(rank1)%ptr%coords(1)=coords(1)
           currk(rank1)%ptr%coords(2)=coords(2)
           currk(rank1)%ptr%coords(3)=0.0
           currk(rank1)%ptr%coords(4)=0.0
!if(rank==0) then
!print*, rank1,currk(rank1)%ptr%coords
!end if
           num(rank1)=num(rank1)+1 
           allocate(currk(rank1)%ptr%next)
           currk(rank1)%ptr=>currk(rank1)%ptr%next
         end if
      end do

    end do

    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    do i=0,size-1
!       deallocate(pic2d%ful2dsend_head(i)%ptr,currk(i)%ptr)
       deallocate(ful2dsend_head(i)%ptr)
       nullify(currk(i)%ptr)
    end do
   
    deallocate(num) 
   end subroutine para_accept_reject_gaussian1d_ful2d_per_per

   subroutine para_accprej_gaus2d2v_fulgyro_unifield_per_per(ful2d_head,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), intent(inout),pointer :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer ::ful2dsend_head,currk
    real8 :: mean(2),y1,y2,py,pymax,x,x1,py1,gmin(2),gmax(2),mu
    int4,dimension(:),pointer :: num
    class(ful2d_node), pointer :: ful2dtmp
    real8 :: theta, sigma, coords(4),rho,mumin,mumax,vperp,tempt
    int4 :: rank1,numcircle,size,rank,comm
    int4 :: ierr,i,j

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    allocate(num(0:size-1))
    allocate(currk(0:size-1),ful2dsend_head(0:size-1))
    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       currk(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mumin=pic2d%para2d%mumin
    mumax=pic2d%para2d%mumax
    pymax=1.0/sqrt(2.0*pi_)/sigma
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if

!    call random_seed()       
    do i=1,pic2d%para2d%numparticle
       x=1.0
       py=0.0
       do while(x.gt.py) 
         y1=gmin(1)+(gmax(1)-gmin(1))*rand()  
         py=exp(-(y1-mean(1))*(y1-mean(1))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x=pymax*rand()
       end do

       x=1.0
       py=0.0
       do while(x.gt.py)
         y2=gmin(2)+(gmax(2)-gmin(2))*rand()        
         py=exp(-(y2-mean(2))*(y2-mean(2))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x=pymax*rand()
       end do

100    mu=-tempt*log(rand())
      if(mu.gt.mumax) then
        goto 100
      end if

      rho=sqrt(2._f64*mu)

     ful2dtmp=>ful2d_head
 !     vperp=sqrt(2._f64*y3)
     do j=0,numcircle-1
           theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=y1+rho*cos(theta)
           coords(2)=y2+rho*sin(theta)      
           coords(3)=rho*cos(theta+pi_/2.0_f64)
           coords(4)=rho*sin(theta+pi_/2.0_f64) 
!           call coordinates_pointoutbound_per_per(coords(1:2),pic2d) 
           rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(rank1==rank) then
           ful2dtmp%coords(1:4)=coords(1:4)
           allocate(ful2dtmp%next)
           ful2dtmp=>ful2dtmp%next
           num(rank)=num(rank)+1
         else  
           currk(rank1)%ptr%coords(1)=coords(1)
           currk(rank1)%ptr%coords(2)=coords(2)
           currk(rank1)%ptr%coords(3)=coords(3)
           currk(rank1)%ptr%coords(4)=coords(3)
!print*, coords(1:2),rank1
           num(rank1)=num(rank1)+1 
           allocate(currk(rank1)%ptr%next)
           currk(rank1)%ptr=>currk(rank1)%ptr%next
         end if
      end do

    end do

    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    do i=0,size-1
       deallocate(ful2dsend_head(i)%ptr)
!       nullify(pic2d%ful2dsend_head(i)%ptr)
       nullify(currk(i)%ptr)
    end do
    
    deallocate(num)
   end subroutine para_accprej_gaus2d2v_fulgyro_unifield_per_per  


   subroutine para_accprej_gaus1d2v_fulgyro_unifield_per_per(ful2d_head,pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer :: ful2dsend_head, currk
    real8 :: mean(2),y1,y2,py,pymax,x,x1,py1,gmin(2),gmax(2),mu
    int4,dimension(:),pointer :: num
    class(ful2d_node), pointer :: ful2dtmp
    real8 :: theta, sigma, coords(4),rho,mumin,mumax,vperp,tempt
    int4 :: rank1,numcircle,size,rank,comm
    int4 :: ierr,i,j

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    allocate(num(0:size-1))
    allocate(currk(0:size-1),ful2dsend_head(0:size-1))
    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       currk(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mumin=pic2d%para2d%mumin
    mumax=pic2d%para2d%mumax
    pymax=1.0/sqrt(2.0*pi_)/sigma
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if

!    call random_seed()       
    do i=1,pic2d%para2d%numparticle
       y2=gmin(1)+(gmax(1)-gmin(1))*rand()    !generate_random_number()

       x=1.0
       py=0.0
       do while(x.gt.py)
         y1=gmin(2)+(gmax(2)-gmin(2))*rand()        
         py=exp(-(y1-mean(2))*(y1-mean(2))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x=pymax*rand()
       end do

100    mu=-tempt*log(rand())
      if(mu.gt.mumax) then
        goto 100
      end if

      rho=sqrt(2._f64*mu)


!           coords(1)=y1
!           coords(2)=y2      
!           coords(3)=0.0
!           coords(4)=0.0 
!           call coordinates_pointoutbound_per_per(coords(1:2),pic2d) 
!           rank1=compute_process_of_point(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
!               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
!           currk(rank1)%ptr%coords(1)=coords(1)
!           currk(rank1)%ptr%coords(2)=coords(2)
!           currk(rank1)%ptr%coords(3)=coords(3)
!           currk(rank1)%ptr%coords(4)=coords(3)
!           num(rank1)=num(rank1)+1 
!           allocate(currk(rank1)%ptr%next)
!           currk(rank1)%ptr=>currk(rank1)%ptr%next
!

     do j=0,numcircle-1
           theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=y1+rho*cos(theta)
           coords(2)=y2+rho*sin(theta)      
           coords(3)=rho*cos(theta+pi_/2.0_f64)
           coords(4)=rho*sin(theta+pi_/2.0_f64) 
!           call coordinates_pointoutbound_per_per(coords(1:2),pic2d) 
           rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           
         if(rank1==rank) then
           ful2dtmp%coords(1:4)=coords(1:4)
           allocate(ful2dtmp%next)
           ful2dtmp=>ful2dtmp%next
           num(rank)=num(rank)+1           
         else
           currk(rank1)%ptr%coords(1)=coords(1)
           currk(rank1)%ptr%coords(2)=coords(2)
           currk(rank1)%ptr%coords(3)=coords(3)
           currk(rank1)%ptr%coords(4)=coords(4)
           num(rank1)=num(rank1)+1 
           allocate(currk(rank1)%ptr%next)
           currk(rank1)%ptr=>currk(rank1)%ptr%next
         end if
      end do

    end do

    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    do i=0,size-1
 !      deallocate(pic2d%ful2dsend_head(i)%ptr,currk(i)%ptr)
       deallocate(ful2dsend_head(i)%ptr)
       nullify(currk(i)%ptr)
    end do
    
    deallocate(num)
   end subroutine para_accprej_gaus1d2v_fulgyro_unifield_per_per  




end module para_random_sample
