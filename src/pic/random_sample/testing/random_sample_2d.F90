program random_sample_2d
#include "work_precision.h"
!  use utilities_module, only: gp_error
!  use piclayout, only: &
!       ful2dsend_node, &
!       ful2d_node
!
 implicit none


    real8 :: gmin=-3, gmax=3
    real8 :: mean,py,pymax,x,delta,length,y,x1,py1
    real8 :: theta,sigma,coords
    real8 :: pi=3.1415926
    int4  :: ierr,number
    real8, dimension(:,:), allocatable :: statis
    int4 :: cell=200
    int4 :: i,j

    allocate(statis(cell,2)) 
    mean=(gmax+gmin)/2.0
    delta=(gmax-gmin)/real(cell,8)
    do i=1, cell
       statis(i,1)=gmin+real(i-1,8)*delta
    end do
    sigma=1.0  
    number=10

    do i=1, number
      pymax=1.0/sqrt(2.0*pi)/sigma
      x=1.0
      py=0.0
      do while(x.gt.py)
         y=gmin+(gmax-gmin)*rand()
         py=exp(-(y-mean)*(y-mean)/2.0/sigma/sigma)/sqrt(2.0*pi)/sigma
         x=pymax*rand()
print*, "x,py=",x,py
         if(x.le.py) then
           goto 100
         end if
      end do
100   x1=1.0
      py1=0.0
      do while(x1.gt.py1)
         y=gmin+(gmax-gmin)*rand()
         py1=exp(-(y-mean)*(y-mean)/2.0/sigma/sigma)/sqrt(2.0*pi)/sigma
         x1=pymax*rand()
print*,"x,py",x,py, "x1,py1", x1,py1
         if(x1.le.py1) then
             goto 150
         end if
      end do     
150     end do


!   length=gmin
!      j=0
!      do while(length.le.y)
!         length=length+delta
!         j=j+1   
!      end do
!      statis(j-1,2)=statis(j-1,2)+1.0
!    end do
!
!    open(10,file="/Users/zhang/program/debug/gnuplot/gauss.txt",status="replace")
!    do i=1, cell
!       write(10,"(2F19.10)") statis(i,:)
!    end do
!
!    close(10)

 end program random_sample_2d
