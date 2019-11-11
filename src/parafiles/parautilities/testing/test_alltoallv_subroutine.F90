program test_alltoallv_subroutine
#include "work_precision.h"
use m_parautilities, only: mpi2d_alltoallv_box_per_per
implicit none
include "mpif.h"

  int4 :: row,myid,numproc(2),globalsz
  real8,dimension(:,:),pointer :: box,rw,re,rn,rs,rsw,rse,rnw,rne
  int4 :: boxindex(4),size
  int4 :: ierr,  num1,num2,num
  int4 :: n1,n2,i,j
  call MPI_Init(ierr)
  call mpi_comm_size(mpi_comm_world,size,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
!  size=9
  row=1
  globalsz=9
!  numproc=NINT(sqrt(real(size,8)))
!  num=3
  numproc(1)=floor(sqrt(real(size,8)+0.05))
  numproc(2)=size/numproc(1)
  print*, numproc
  n1=modulo(myid,numproc(1))
   n2=(myid-n1)/numproc(1)
  num=globalsz/numproc(1)

  boxindex(1)=n1*num+1
  boxindex(2)=n1*num+num
  boxindex(3)=n2*num+1
  boxindex(4)=n2*num+num
  allocate(box(num,num),rw(num,row),re(num,row),rn(row,num),rs(row,num),rsw(row,row),&
           rse(row,row),rnw(row,row),rne(row,row))

  do i=1,boxindex(2)-boxindex(1)+1
     do j=1,boxindex(4)-boxindex(3)+1
         box(i,j)=i+j
     end do
  end do
 call mpi2d_alltoallv_box_per_per(row,mpi_comm_world,myid,numproc, &
                                       box,rw,re,rn,rs,rsw,rse,rnw,rne,boxindex)

 if(myid==2) then
    print*, "rs=",rs
    print*, "rn=",rn
    print*, "rne=",rne
    print*, "re:", re
 end if
  call MPI_FINALIZE(IERR)
  print*, 28

 end program

