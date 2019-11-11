program random_sample_mu
#include "work_precision.h"
!  use utilities_module, only: gp_error
!  use piclayout, only: &
!       ful2dsend_node, &
!       ful2d_node
!
 implicit none


    real8 :: mean,py,pymax,x,delta,length,y
    real8 :: theta,coords
    real8 :: pi=3.1415926
    int4  :: ierr,number
    real8, dimension(:,:), allocatable :: statis
    int4 :: cell=200
    int4 :: i,j,fileitem=10
    character(100) :: filename
    logical :: alive    
    real8 :: mu,mumax=30.0,mumin=0.0
    real8 :: tempt=2.0



    allocate(statis(cell,2)) 
    delta=(mumax-mumin)/real(cell,8)
    do i=1,cell
      statis(i,1)=mumin+real(i-1,8)*delta
!    print*, statis(i,1)
    end do

    number=300000

    do i=1, number
100      mu=-tempt*log(rand())
      if(mu.gt.mumax) then
        goto 100
      end if 
!      print*, "mu=",mu     
      length=mumin
      j=0
      do while(length.le.mu)
         length=length+delta
         j=j+1   
      end do
      statis(j,2)=statis(j,2)+1.0
    end do

    filename="/PARA/blsc950/electrostatic_exp/run/gauss.txt"
      
    inquire(file=trim(filename),exist=alive)
    if(.not.alive) then
       open(fileitem,file=trim(filename),status='new')
    else
       open(fileitem,file=trim(filename),status='replace')
    end if

    do i=1, cell
       write(fileitem,"(2F19.10)") statis(i,1),statis(i,2)
       print*, statis(i,1)
    end do
 
    close(fileitem)

 end program random_sample_mu
