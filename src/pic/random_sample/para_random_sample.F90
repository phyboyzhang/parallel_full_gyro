module para_random_sample
#include "work_precision.h"
  use omp_lib
  use utilities_module, only: gp_error
  use piclayout, only: &
       ful2dsend_node,ful2d_node, gy2dmu_node, gy2d_node, &
       parameters_array_2d
  use paradata_type, only: pic_para_total2d_base

  use constants, only: pi_
  use m_moveparticles, only: mpi2d_alltoallv_box_per_per 
  use paradata_utilities, only: compute_process_of_point_per_per
  use m_picutilities, only: mpi2d_alltoallv_send_particle_2d, &
                            mpi2d_alltoallv_send_particle_2d_gy
  use paradata_utilities, only: coordinates_pointoutbound_per_per
  use piclayout, only: ful2d_node, gy2d_node,ful2dsend_node,  &
                        gy2dsend_node
  use congruence_sampling_seeds, only: congruence, lcrg_anbn, lcrg_evaluate
  use field_initialize, only: test_trigonfun
!  use orbit_data_base, only: ful2d_node,ful2dsend_node
  implicit none

  public:: para_accept_reject_gaussian1d_ful2d_per_per, &
           para_accprej_gaus2d2v_fulgyro_unifield_per_per, &
           para_accprej_gaus1d2v_fulgyro_unifield_per_per, &
           congru_accprej_gaus1d2v_fulgyro_unifield_per_per, &
   !        congru_accprej_flatpert2d2v_fulgyro_unifield_per_per, &
           congru_accprej_2d2v_fulgyro_unifield_per_per, &
           congru_accprej_2d2v_per_per_ful
  
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

    if(.not.associated(ful2d_head)) then
      print*, "pic2d%ful2d_head is not allocated"
      stop
    end if

!    call random_seed()       
    do i=1,pic2d%para2d%numparticles
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
    do i=1,pic2d%para2d%numparticles
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


   subroutine para_accprej_gaus1d2v_fulgyro_unifield_per_per(ful2d_head,gy2dmu_head,pic2d,pamearray)
    class(pic_para_total2d_base), pointer :: pic2d
    class(parameters_array_2d), pointer, intent(in) :: pamearray
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer :: ful2dsend_head, fulcur
!    class(gy2d_node), pointer, intent(inout) :: gy2d_head
    class(gy2dmu_node), dimension(:), pointer,intent(inout) :: gy2dmu_head
    class(gy2dsend_node), dimension(:), pointer :: gy2dsend_head,gycur
    real8 :: mean(2),y(2),yf(2),py,pymax,x,x1,py1,gmin(2),gmax(2),mu
    int4,dimension(:),pointer :: num
    int4,dimension(:), pointer :: num_gy
    class(ful2d_node), pointer :: ful2dtmp
    class(gy2d_node), pointer :: gy2dtmp
    class(gy2dmu_node), dimension(:), pointer :: gy2dmutmp
!    int4, dimension(:), pointer :: munum_partion
    real8 :: theta, sigma, coords(4),rho,mumin,mumax,vperp,tempt
    real8 :: integ
    int4 :: rank1,numcircle,size,rank,comm,mu_num
    int4 :: ierr,i,j,k

    int4 :: sum=0

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    mu_num=pic2d%para2d%mu_num
    allocate(num(0:size-1),num_gy(0:size-1))
    allocate(fulcur(0:size-1),ful2dsend_head(0:size-1))
    allocate(gy2dmutmp(1:mu_num))
!    allocate(munum_partion(1:pic2d%para2d%mu_num))
    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       fulcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    ful2dtmp=>ful2d_head
 
    do i=1,mu_num
       allocate(gy2dmu_head(i)%ptr)
       gy2dmutmp(i)%ptr=>gy2dmu_head(i)%ptr     
    end do

    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mumin=pic2d%para2d%mumin
    mumax=pic2d%para2d%mumax
!    pymax=1.0/sqrt(2.0*pi_)/sigma
    pymax=1.0
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if
!    if(.not.associated(gy2d_head)) then
!      print*, "gy2d_head is not allocated"
!      stop
!    end if

  do j=1,mu_num    
    num_gy=0
    allocate(gy2dsend_head(0:size-1))
    allocate(gycur(0:size-1)) 
    do k=0,size-1
       allocate(gy2dsend_head(k)%ptr)
       gycur(k)%ptr=>gy2dsend_head(k)%ptr
    end do

    do i=1, pamearray%munum_partition(j) 
       y(2)=gmin(1)+(gmax(1)-gmin(1))*rand()    !generate_random_number()

       x=1.0
       py=0.0
       do while(x.gt.py)
         y(1)=gmin(2)+(gmax(2)-gmin(2))*rand()        
    !     py=exp(-(y(1)-mean(2))*(y(1)-mean(2))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         py=exp(-(y(1)-mean(2))*(y(1)-mean(2))/sigma)
         x=pymax*rand()
       end do
       yf=y
      
      rank1=compute_process_of_point_per_per(y,pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)

      if(rank1==rank) then
         gy2dmutmp(j)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
         allocate(gy2dmutmp(j)%ptr%next)
         gy2dmutmp(j)%ptr=>gy2dmutmp(j)%ptr%next
         num_gy(rank)=num_gy(rank)+1

       else 

         gycur(rank1)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
         allocate(gycur(rank1)%ptr%next)
         gycur(rank1)%ptr=>gycur(rank1)%ptr%next
         num_gy(rank1)=num_gy(rank1)+1
       end if

 
     rho=sqrt(2._f64*pamearray%mu_nodes(j))
     do k=0,numcircle-1

          theta=real(k,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=yf(1)+rho*cos(theta)
           coords(2)=yf(2)+rho*sin(theta)      
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
           fulcur(rank1)%ptr%coords(1)=coords(1)
           fulcur(rank1)%ptr%coords(2)=coords(2)
           fulcur(rank1)%ptr%coords(3)=coords(3)
           fulcur(rank1)%ptr%coords(4)=coords(4)
           num(rank1)=num(rank1)+1 
           allocate(fulcur(rank1)%ptr%next)
           fulcur(rank1)%ptr=>fulcur(rank1)%ptr%next
         end if
      end do

    end do


    call mpi2d_alltoallv_send_particle_2d_gy(gy2dmu_head,gy2dsend_head,num_gy,j,pic2d)

    deallocate(gy2dsend_head,gycur)

    end do  !! end loop denoted by j
 
    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    deallocate(ful2dsend_head)
 
   
    deallocate(num,num_gy)
    deallocate(gy2dmutmp)
!    deallocate(munum_partion)
   end subroutine para_accprej_gaus1d2v_fulgyro_unifield_per_per  


   subroutine congru_accprej_gaus1d2v_fulgyro_unifield_per_per(ful2d_head,gy2dmu_head,pic2d,pamearray)
    class(pic_para_total2d_base), pointer :: pic2d
    class(parameters_array_2d), pointer, intent(in) :: pamearray
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer :: ful2dsend_head, fulcur
!    class(gy2d_node), pointer, intent(inout) :: gy2d_head
    class(gy2dmu_node), dimension(:), pointer,intent(inout) :: gy2dmu_head
    class(gy2dsend_node), dimension(:), pointer :: gy2dsend_head,gycur
    real8 :: mean(2),y(2),yf(2),py,pymax,x,x1,py1,gmin(2),gmax(2),mu
    int4,dimension(:),pointer :: num
    int4,dimension(:), pointer :: num_gy
    class(ful2d_node), pointer :: ful2dtmp
    class(gy2d_node), pointer :: gy2dtmp
    class(gy2dmu_node), dimension(:), pointer :: gy2dmutmp
!    int4, dimension(:), pointer :: munum_partion
    real8 :: theta, sigma, coords(4),rho,mumin,mumax,vperp,tempt
    real8 :: integ
    int4 :: rank1,numcircle,size,rank,comm,mu_num
    int4 :: ierr,i,j,k

    int4 :: sum=0

!!!!!! for congurence sampling
    integer, parameter :: numsample=10000000
    integer :: a
    integer, pointer :: an
    integer :: b
    integer, pointer :: bn
    integer :: c
    integer :: error
    integer :: id
    integer :: h,m
    integer :: k_hi
    integer :: u
    integer :: v
    real(8) :: ratio

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    mu_num=pic2d%para2d%mu_num
    allocate(num(0:size-1),num_gy(0:size-1))
    allocate(fulcur(0:size-1),ful2dsend_head(0:size-1))
    allocate(gy2dmutmp(1:mu_num))
    allocate(an,bn)
!    allocate(munum_partion(1:pic2d%para2d%mu_num))
    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       fulcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    ful2dtmp=>ful2d_head
 
    do i=1,mu_num
       allocate(gy2dmu_head(i)%ptr)
       gy2dmutmp(i)%ptr=>gy2dmu_head(i)%ptr     
    end do

    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mumin=pic2d%para2d%mumin
    mumax=pic2d%para2d%mumax
!    pymax=1.0/sqrt(2.0*pi_)/sigma
    pymax=1.0
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if
!    if(.not.associated(gy2d_head)) then
!      print*, "gy2d_head is not allocated"
!      stop
!    end if


!!!!!!+++++++++++++++++++++
!!!! The begining of congurence sampling
    a = 16807
    b = 0
    c = 2147483647
   
    k_hi = size*numsample

    call lcrg_anbn(a, b, c, size,an, bn)

    v = 12345
    h = 1

    do while(h.le.rank)
      u=v
      v=lcrg_evaluate(a,b,c,u)
      h=h+1
    end do

    m=rank+size
   
    

! do while(m.le.k_hi)
  do j=1,mu_num    
    num_gy=0
    allocate(gy2dsend_head(0:size-1))
    allocate(gycur(0:size-1)) 
    do k=0,size-1
       allocate(gy2dsend_head(k)%ptr)
       gycur(k)%ptr=>gy2dsend_head(k)%ptr
    end do

    do i=1, pamearray%munum_partition(j) 
        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is nesseary."
          stop 
        endif
        ratio=real(v,8)/real(c,8)
        y(2)=gmin(2)+(gmax(2)-gmin(2))*ratio    !generate_random_number()
      
       x=1.0
       py=0.0      
       do while(x.gt.py)
         u=v
         v=lcrg_evaluate(an,bn,c,u)
         m=m+size
         if(m.ge.k_hi) then
           print*, "#ERROR: the inital given sampling number is nesseary."
           stop
         endif
         ratio=real(v,8)/real(c,8)         
         y(1)=gmin(1)+(gmax(1)-gmin(1))*ratio        
         py=exp(-(y(1)-mean(1))*(y(1)-mean(1))/sigma)

         u=v
         v=lcrg_evaluate(an,bn,c,u)
         m=m+size
         if(m.ge.k_hi) then
           print*, "#ERROR: the inital given sampling number is nesseary."
           stop
         endif
         ratio=real(v,8)/real(c,8)
         x=pymax*ratio
       end do
       yf=y
      
      rank1=compute_process_of_point_per_per(y,pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)

      if(rank1==rank) then
         gy2dmutmp(j)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
         allocate(gy2dmutmp(j)%ptr%next)
         gy2dmutmp(j)%ptr=>gy2dmutmp(j)%ptr%next
         num_gy(rank)=num_gy(rank)+1

       else 

         gycur(rank1)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
         allocate(gycur(rank1)%ptr%next)
         gycur(rank1)%ptr=>gycur(rank1)%ptr%next
         num_gy(rank1)=num_gy(rank1)+1
       end if

 
     rho=sqrt(2._f64*pamearray%mu_nodes(j))
     do k=0,numcircle-1

          theta=real(k,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=yf(1)+rho*cos(theta)
           coords(2)=yf(2)+rho*sin(theta)      
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
           fulcur(rank1)%ptr%coords(1)=coords(1)
           fulcur(rank1)%ptr%coords(2)=coords(2)
           fulcur(rank1)%ptr%coords(3)=coords(3)
           fulcur(rank1)%ptr%coords(4)=coords(4)
           num(rank1)=num(rank1)+1 
           allocate(fulcur(rank1)%ptr%next)
           fulcur(rank1)%ptr=>fulcur(rank1)%ptr%next
         end if
      end do
!      if(rank==0) then
!        print*, "# The sampling for the", j, "th mu is finished."
!      endif
    end do


    call mpi2d_alltoallv_send_particle_2d_gy(gy2dmu_head,gy2dsend_head,num_gy,j,pic2d)

    deallocate(gy2dsend_head,gycur)

    end do  !! end loop denoted by j
 
    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    deallocate(ful2dsend_head)
 
   
    deallocate(num,num_gy)
    deallocate(gy2dmutmp)
!    deallocate(munum_partion)
  end subroutine congru_accprej_gaus1d2v_fulgyro_unifield_per_per  

  subroutine congru_accprej_2d2v_fulgyro_unifield_per_per(ful2d_head,gy2dmu_head,pic2d,pamearray,eqdistr)
    class(pic_para_total2d_base), pointer :: pic2d
    class(parameters_array_2d), pointer, intent(in) :: pamearray
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer :: ful2dsend_head, fulcur
!    class(gy2d_node), pointer, intent(inout) :: gy2d_head
    class(gy2dmu_node), dimension(:), pointer,intent(inout) :: gy2dmu_head
    character(len=*),intent(in) :: eqdistr

    class(gy2dsend_node), dimension(:), pointer :: gy2dsend_head,gycur
    real8 :: mean(2),y(2),yf(2),py,pymax,x,x1,py1,gmin(2),gmax(2),mu
    int4,dimension(:),pointer :: num
    int4,dimension(:), pointer :: num_gy
    class(ful2d_node), pointer :: ful2dtmp
    class(gy2d_node), pointer :: gy2dtmp
    class(gy2dmu_node), dimension(:), pointer :: gy2dmutmp
!    int4, dimension(:), pointer :: munum_partion
    real8 :: theta, sigma, coords(4),rho,mumin,mumax,vperp,tempt
    real8 :: integ
    int4 :: rank1,numcircle,size,rank,comm,mu_num
    int4 :: ierr,i,j,k

    int4 :: sum=0

!!!!!! for congurence sampling
    integer, parameter :: numsample=10000000
    integer :: a
    integer, pointer :: an
    integer :: b
    integer, pointer :: bn
    integer :: c
    integer :: error
    integer :: id
    integer :: h,m
    integer :: k_hi
    integer :: u
    integer :: v
    real(8) :: ratio

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    mu_num=pic2d%para2d%mu_num
    allocate(num(0:size-1),num_gy(0:size-1))
    allocate(fulcur(0:size-1),ful2dsend_head(0:size-1))
    allocate(gy2dmutmp(1:mu_num))
    allocate(an,bn)
!    allocate(munum_partion(1:pic2d%para2d%mu_num))
    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       fulcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    ful2dtmp=>ful2d_head
 
    do i=1,mu_num
       allocate(gy2dmu_head(i)%ptr)
       gy2dmutmp(i)%ptr=>gy2dmu_head(i)%ptr     
    end do

    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mumin=pic2d%para2d%mumin
    mumax=pic2d%para2d%mumax
!    pymax=1.0/sqrt(2.0*pi_)/sigma
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if

!!!!!!+++++++++++++++++++++
!!!! The begining of congurence sampling
    a = 16807
    b = 0
    c = 2147483647
   
    k_hi = size*numsample

    call lcrg_anbn(a, b, c, size,an, bn)

    v = 12345
    h = 1

    do while(h.le.rank)
      u=v
      v=lcrg_evaluate(a,b,c,u)
      h=h+1
    end do

    m=rank+size 
 
    pymax=1.0+2.0*pic2d%para2d%amp

  do j=1,mu_num    
    num_gy=0
    allocate(gy2dsend_head(0:size-1))
    allocate(gycur(0:size-1)) 
    do k=0,size-1
       allocate(gy2dsend_head(k)%ptr)
       gycur(k)%ptr=>gy2dsend_head(k)%ptr
    end do

    do i=1, pamearray%munum_partition(j) 
!      select case(distri)
!        case("trig")    
!          call congru_sampling_kernel_trigonometry(v,size,k_hi,m,an,bn,c,y,yf,pymax,pic2d)
!        case("flat")
!          call congru_sampling_kernel_flat(v,size,k_hi,m,an,bn,c,y,yf,pymax,pic2d)
!        case default
!          stop
!      end select
          call congru_sampling_kernel_trig_fulgy(v,size,k_hi,m,an,bn,c,y,yf,pymax,mean,eqdistr,pic2d)

        rank1=compute_process_of_point_per_per(y,pic2d%para2d%numproc,pic2d%para2d%gxmin, &
                pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)

        if(rank1==rank) then
          gy2dmutmp(j)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
          allocate(gy2dmutmp(j)%ptr%next)
          gy2dmutmp(j)%ptr=>gy2dmutmp(j)%ptr%next
          num_gy(rank)=num_gy(rank)+1

        else 

          gycur(rank1)%ptr%coords(1:3)=(/y(1),y(2),pamearray%mu_nodes(j)/)
          allocate(gycur(rank1)%ptr%next)
          gycur(rank1)%ptr=>gycur(rank1)%ptr%next
          num_gy(rank1)=num_gy(rank1)+1
        end if

 
        rho=sqrt(2._f64*pamearray%mu_nodes(j))
        do k=0,numcircle-1

           theta=real(k,8)*2.0_f64*pi_/real(numcircle,8)   
           coords(1)=yf(1)+rho*cos(theta)
           coords(2)=yf(2)+rho*sin(theta)      
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
            fulcur(rank1)%ptr%coords(1)=coords(1)
            fulcur(rank1)%ptr%coords(2)=coords(2)
            fulcur(rank1)%ptr%coords(3)=coords(3)
            fulcur(rank1)%ptr%coords(4)=coords(4)
            num(rank1)=num(rank1)+1 
            allocate(fulcur(rank1)%ptr%next)
            fulcur(rank1)%ptr=>fulcur(rank1)%ptr%next
          end if
        end do

    end do


    call mpi2d_alltoallv_send_particle_2d_gy(gy2dmu_head,gy2dsend_head,num_gy,j,pic2d)

    deallocate(gy2dsend_head,gycur)

    end do  !! end loop denoted by j
 
    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    deallocate(ful2dsend_head)
 
   
    deallocate(num,num_gy)
    deallocate(gy2dmutmp)
!    deallocate(munum_partion)
   end subroutine congru_accprej_2d2v_fulgyro_unifield_per_per  


  subroutine congru_accprej_2d2v_per_per_ful(ful2d_head,pic2d,eqdistr,vmin,vmax)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    class(ful2dsend_node), dimension(:),pointer :: ful2dsend_head, fulcur
    character(len=*),intent(in) :: eqdistr
    real8 :: mean(2),yf(4),py,pymax(2),x,x1,py1,gmin(2),gmax(2),mu
    real8,intent(in) :: vmin(2),vmax(2)
    int4,dimension(:),pointer :: num
    class(ful2d_node), pointer :: ful2dtmp
    real8 :: sigma, coords(4),rho,mumin,mumax,vperp,tempt
    real8 :: integ
    int4 :: rank1,numcircle,size,rank,comm
    int4 :: ierr,i,j,k

    int4 :: sum=0

!!!!!! for congurence sampling
    integer, parameter :: numsample=20000000
    integer :: a
    integer, pointer :: an
    integer :: b
    integer, pointer :: bn
    integer :: c
    integer :: error
    integer :: id
    integer :: h,m
    integer :: k_hi
    integer :: u
    integer :: v
    real(8) :: ratio

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    numcircle=pic2d%para2d%numcircle
    sigma=pic2d%para2d%sigma
    tempt=pic2d%para2d%tempt
    allocate(num(0:size-1))
    allocate(fulcur(0:size-1),ful2dsend_head(0:size-1))
    allocate(an,bn)

    do i=0,size-1
       allocate(ful2dsend_head(i)%ptr)
       fulcur(i)%ptr=>ful2dsend_head(i)%ptr
    end do
    ful2dtmp=>ful2d_head
 
    num=0
    gmin=pic2d%para2d%gxmin
    gmax=pic2d%para2d%gxmax
    mean(1)=(pic2d%para2d%gxmax(1)+pic2d%para2d%gxmin(1))/2.0_f64
    mean(2)=(pic2d%para2d%gxmax(2)+pic2d%para2d%gxmin(2))/2.0_f64      

    if(.not.associated(ful2d_head)) then
      print*, "ful2d_head is not allocated"
      stop
    end if

!!!!!!+++++++++++++++++++++
!!!! The begining of congurence sampling
    a = 16807
    b = 0
    c = 2147483647
   
    k_hi = size*numsample

    call lcrg_anbn(a, b, c, size,an, bn)

    v = 12345
    h = 1

    do while(h.le.rank)
      u=v
      v=lcrg_evaluate(a,b,c,u)
      h=h+1
    end do

    m=rank+size 
 
    pymax(1)=1.0+2.0*pic2d%para2d%amp
    pymax(2)=1.0/sqrt(2.0*pi_)

    do i=1, pic2d%para2d%numparticles     

          call congru_sampling_kernel_trig_ful(v,size,k_hi,m,an,bn,c,yf,pymax,mean,vmin,vmax,eqdistr,pic2d)
          rank1=compute_process_of_point_per_per(yf(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)           

          if(rank1==rank) then
            ful2dtmp%coords(1:4)=yf(1:4)
            allocate(ful2dtmp%next)
            ful2dtmp=>ful2dtmp%next
            num(rank)=num(rank)+1           
          else
            fulcur(rank1)%ptr%coords(1:4)=yf(1:4)
            num(rank1)=num(rank1)+1 
            allocate(fulcur(rank1)%ptr%next)
            fulcur(rank1)%ptr=>fulcur(rank1)%ptr%next
          end if
     end do


    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

    deallocate(ful2dsend_head)
 
   
    deallocate(num)
  end subroutine congru_accprej_2d2v_per_per_ful  


   subroutine congru_sampling_kernel_trig_fulgy(v,size,k_hi,m,an,bn,c,y,yf,pmax,mean,eqdistr,pic2d) 
     class(pic_para_total2d_base), intent(in) :: pic2d
     int4, intent(in) :: size,k_hi,an,bn,c
     real8, intent(in) :: mean(2)
     character(len=*), intent(in) :: eqdistr
     int4 :: m,v
     real(8), intent(in) :: pmax
     real(8), intent(inout) :: y(2),yf(2)

     real(8) :: gmin(2),gmax(2),amp,waveone,wavetwo
     int4 :: u
     real(8) :: ratio,x,py

     gmin=pic2d%para2d%gxmin
     gmax=pic2d%para2d%gxmax
     amp=pic2d%para2d%amp
     waveone=pic2d%para2d%waveone
     wavetwo=pic2d%para2d%wavetwo
 
!     select case(eqdistr)
!       case("gaussian")
        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
        endif
        ratio=real(v,8)/real(c,8)
        y(2)=gmin(2)+(gmax(2)-gmin(2))*ratio    !generate_random_number()
      x=1.0
      py=0.0
      do while(x.gt.py)
        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
        endif
        ratio=real(v,8)/real(c,8)
        y(1)=gmin(1)+(gmax(1)-gmin(1))*ratio
        py=test_trigonfun(y,amp,waveone,wavetwo,pic2d%para2d%sigma,mean,eqdistr)

        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
        endif
        ratio=real(v,8)/real(c,8)        
        x=pmax*ratio
      end do
        yf=y


   end subroutine congru_sampling_kernel_trig_fulgy

   subroutine congru_sampling_kernel_flat(v,size,k_hi,m,an,bn,c,y,yf,pmax,pic2d)
     class(pic_para_total2d_base), intent(in) :: pic2d
     int4, intent(in) :: size,k_hi,an,bn,c
     int4 :: m,v
     real(8), intent(in) :: pmax
     real(8), intent(inout) :: y(2),yf(2)

     real(8) :: gmin(2),gmax(2),amp,waveone,wavetwo
     int4 :: u
     real(8) :: ratio,x,py

        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
        endif
        ratio=real(v,8)/real(c,8)
        y(2)=gmin(2)+(gmax(2)-gmin(2))*ratio

        u=v
        v=lcrg_evaluate(an,bn,c,u)
        m=m+size
        if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
        endif
        ratio=real(v,8)/real(c,8)
        y(1)=gmin(1)+(gmax(1)-gmin(1))*ratio

        yf=y
   end subroutine congru_sampling_kernel_flat

   subroutine congru_sampling_kernel_trig_ful(v,size,k_hi,m,an,bn,c,yf,pymax,mean,vmin,vmax,eqdistr,pic2d) 
     class(pic_para_total2d_base), intent(in) :: pic2d
     int4, intent(in) :: size,k_hi,an,bn,c
     real8, intent(in) :: mean(2),vmin(2),vmax(2)
     character(len=*), intent(in) :: eqdistr
     int4 :: m,v
     real(8), intent(in) :: pymax(2)
     real(8), intent(inout) :: yf(4)

     real(8) :: gmin(2),gmax(2),amp,waveone,wavetwo,y(4)
     int4 :: u
     real(8) :: ratio,x,py

     gmin=pic2d%para2d%gxmin
     gmax=pic2d%para2d%gxmax
     amp=pic2d%para2d%amp
     waveone=pic2d%para2d%waveone
     wavetwo=pic2d%para2d%wavetwo
 
     call congru_seed(v,an,bn,c,m,size,k_hi) 
     ratio=real(v,8)/real(c,8)
     y(1)=gmin(2)+(gmax(2)-gmin(2))*ratio    !generate_random_number()

      x=1.0
      py=0.0
      do while(x.gt.py)
        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)
        y(2)=gmin(1)+(gmax(1)-gmin(1))*ratio
        py=test_trigonfun(y(1:2),amp,waveone,wavetwo,pic2d%para2d%sigma,mean,eqdistr)

        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)        
        x=pymax(1)*ratio
      end do

      x=1.0 
      py=0.0
      do while(x.gt.py)
        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)
        y(3)=vmin(1)+(vmax(1)-vmin(1))*ratio
        py=exp(-y(3)**2)/sqrt(2.0*pi_)
        
        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)
        x=pymax(2)*ratio
      end do

      x=1.0
      py=0.0
      do while(x.gt.py)
        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)
        y(4)=vmin(1)+(vmax(1)-vmin(1))*ratio
        py=exp(-y(4)**2)/sqrt(2.0*pi_)

        call congru_seed(v,an,bn,c,m,size,k_hi)
        ratio=real(v,8)/real(c,8)
        x=pymax(2)*ratio
      end do
       
        yf=y


   end subroutine congru_sampling_kernel_trig_ful

   
   subroutine congru_seed(v,an,bn,c,m,size,k_hi)
     int4, intent(in) :: an,bn,c,size,k_hi
     int4, intent(inout) :: v,m
     int4 :: u
     
     u=v
     v=lcrg_evaluate(an,bn,c,u)
     m=m+size
     if(m.ge.k_hi) then
          print*, "#ERROR: the inital given sampling number is not big enough."
          stop
     endif

   end subroutine




end module para_random_sample
