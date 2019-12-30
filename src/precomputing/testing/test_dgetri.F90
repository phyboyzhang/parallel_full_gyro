program test_dgetri
#include "work_precision.h"
implicit none

 int4 :: NUMDIM,LDA,INFO,LWORK
 int4, dimension(:), allocatable :: ipiv
 real8, dimension(:), allocatable :: work
 real8, dimension(:,:), allocatable :: buf,buf1
 int4 :: i,j

 NUMDIM=1000
 LDA=numdim
 LWORK=numdim*64
 allocate(IPIV(numdim))
 allocate(work(lwork))
 allocate(buf(numdim,numdim))

do i=1,numdim
  do j=1,numdim
    if(i==j) then
      buf(i,j)=4.0
    else
      buf(i,j)=0.0
    endif
  enddo
enddo

! buf=reshape((/2.0, 3.0, 4.0,5.0, &
!           6.0, 7.0, 8.0,9.0, &
!           10.0, 4.0, 2.3, 3.4, &
!           7.0, 4.0, 5.0, 3.0/) ,shape(buf))
!buf1=buf
 call dgetrf(numdim,numdim,buf,LDA,IPIV,INFO)

 call dgetri(numdim,buf,LDA,IPIV,WORK,LWORK,INFO)

!  print*,"buf=", buf
!
! buf=matmul(buf,buf1)

!  print*, "buf=",buf

  do i=1, numdim
    do j=1, numdim
      if(i==j) then
        print*, buf(i,j)
      endif 
    enddo
  enddo

end program

