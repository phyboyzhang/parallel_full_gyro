program test_pic_noise
#include "work_precision.h"
use piclayout,only: ful2d_node
use cartesian_mesh,only: cartesian_mesh_1d, &
                         init_para_cartesian_mesh_1d
use m_picutilities, only: singlepart_to_mesh 

  class(cartesian_mesh_1d), pointer :: m_x1, m_x2
  int4   :: nodes_x1, nodes_x2
  real8 :: gmin(2), gmax(2)
  real8 :: delta(2)
  real8 :: x,y1,y2,mu,py,pymax,mean(2),x1(2)
  real8,dimension(:,:),pointer :: density
  real8 :: pi_=3.1415926
  int4 :: numparticle=50000
  real8 :: sigma=1.0,tempt=1.5,mumax=20.0
 
  character(100) :: filename 
  logical :: alive

   nodes_x1=15
  nodes_x2=15
  gmin=(/0.0,0.0/)
  gmax=(/3.0,3.0/)

  allocate(density(nodes_x1,nodes_x2))
  density=0.0_f64

    pymax=1.0/sqrt(2.0*pi_)/sigma
    mean(1)=(gmax(1)+gmin(1))/2.0_f64
    mean(2)=(gmax(2)+gmin(2))/2.0_f64

  delta(1)=(gmax(1)-gmin(1))/real(nodes_x1-1,8)
  delta(2)=(gmax(2)-gmin(2))/real(nodes_x2-1,8)
  
  m_x1=>init_para_cartesian_mesh_1d(nodes_x1,gmin(1),gmax(1),delta(1))
  m_x2=>init_para_cartesian_mesh_1d(nodes_x1,gmin(2),gmax(2),delta(2))

   do i=1,numparticle
       x=1.0
       py=0.0
       do while(x.gt.py)
         y1=gmin(1)+(gmax(1)-gmin(1))*rand()   
         py=exp(-(y1-mean(1))*(y1-mean(1))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
         x=pymax*rand()
       end do

!       x=1.0
!       py=0.0
!       do while(x.gt.py)
!         y2=gmin(2)+(gmax(2)-gmin(2))*rand()
!         py=exp(-(y2-mean(2))*(y2-mean(2))/2.0/sigma/sigma)/sqrt(2.0*pi_)/sigma
!         x=pymax*rand()
!       end do
       y2=gmin(1)+(gmax(1)-gmin(1))*rand()

100    mu=-tempt*log(rand())
      if(mu.gt.mumax) then
        goto 100
      end if

      x1(1)=y1
      x1(2)=y2
      call  singlepart_to_mesh(density,x1,delta,gmin)
     
    end do

  filename="/PARA/blsc950/electrostatic_exp/run/density_series.txt"
  inquire(file=trim(filename),exist=alive)
    if(.not.alive) then
       open(10,file=trim(filename),status='new')
    else
       open(10,file=trim(filename),status='replace')
    end if   
    do j=1,nodes_x2
       y2=gmin(2)+real(j-1,8)*delta(2)
       do i=1,nodes_x1
          y1=gmin(1)+real(i-1,8)*delta(1)
          if(i==1)  then
             write(10,*)
          end if
             write(10,"(3F19.9)") y2,density(i,j),y1    
         !    print*, y2,density(i,j),y1 
       end do
    end do



end program test_pic_noise
