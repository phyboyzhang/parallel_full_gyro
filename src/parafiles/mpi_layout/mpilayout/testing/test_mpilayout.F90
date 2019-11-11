program test_mpilayout
#include "work_precision.h"
use m_mpilayout, only: initialize_layout_with_distributed_2d_array, &
                       split_domain_into_box1d
  implicit none
  include "mpif.h"

  int4 :: ierr,myid,nproc
  int4 :: numproc1=4,numproc2
  int4, dimension(:,:), allocatable :: intervals_x1,intervals_x2

  int4 :: global_npx1=9
  int4 :: i

  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!  numproc1=floor(sqrt(real(nproc,8)+0.01))
   numproc1=nproc
   ALLOCATE( intervals_x1(0:numproc1-1,2), stat=ierr )
!  ALLOCATE( intervals_x2(0:numproc2-1,2), stat=ierr )
!  if(myid==0) then
!print*, 1
!  end if
  intervals_x1(0:numproc1-1,1:2) = &
         split_domain_into_box1d( 1, global_npx1, numproc1 )
   if(myid==0) then
       print*,numproc1
   end if
  if(myid==0) then
  do i=0, numproc1-1
    print*, "i=",i, intervals_x1(i,:)
  end do
  end if
!if(myid==0) then
!  print*,
!end if
call mpi_finalize(ierr)

end program
