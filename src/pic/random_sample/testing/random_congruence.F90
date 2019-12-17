program  para_random_congruence
!use congruence_sampling_seeds, only: congruence, lcrg_anbn, lcrg_evaluate
implicit none
  include "mpif.h"


    integer, parameter :: num=500000, numcel=300, &
                          pnum=10000 
    real(8), parameter :: minn=0.0,maxx=10.0

    integer :: a
    integer, pointer :: an
    integer :: b
    integer, pointer :: bn
    integer :: c
    integer :: error
    integer :: id
    integer :: j
    integer :: k
    integer :: k_hi
    integer :: p
    integer :: u
    integer :: v
    real(8), dimension(:), pointer :: array
    real(8), dimension(:,:), pointer :: summ

    real(8) :: dx=maxx/real(numcel,8)
    character(200) :: filepath="/home/qmlu/zsx163/parallel_full_gyro/run/distr.txt"
    integer :: i,ierr,lowind

    real(8) :: pi=3.14125, sigma=1.0,mean=5.0
    real(8) :: x,y,py,pymax,delta,ratio,lbound,ubound    
    real(8), dimension(:), pointer :: buf

    call MPI_Init(error)
    call MPI_Comm_size(MPI_COMM_WORLD, p, error)
    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

    allocate(an,bn) 
    allocate(array(numcel+1), summ(2,numcel+1),buf(numcel+1))

    array=0.0

    do i=1, numcel
      summ(1,i)=dx*real(i-1,8)
    end do

    pymax=1.0/sqrt(2.0*pi)/sigma


    a = 16807
    b = 0
    c = 2147483647

    if(id==0) then
       print*,
       print*, " LCRG parameters:"
       print*, 
       print*, " A = ", a
       print*, " B = ", b
       print*, " C = ", c
 
    end if

    k_hi = p*num
     
    call lcrg_anbn(a, b, c, p,an, bn)

    if(id==0) then
      print*, 
      print*, "LCRG parameters for P processor:"
      print*,
      print*, "AN = ", an
      print*, "BN = ", bn
      print*, "C = ", c
      print*, 
      print*, "Have all the processors participate in computing"
      print*, "the same random number sequence."
      print*, 
      print*, "   K        ID                Input           Output"
      print*
    endif

    v = 12345
    j=1
    do while(j.le.id)
      u=v
      v=lcrg_evaluate(a,b,c,u)
      j=j+1
    end do
    k=id


!    call interp(array,v,c,maxx,dx)
    k=id+p
    j=0
    do while(k.le.k_hi)
      if(j.ge.pnum) then
        exit
      else
        u=v
        v=lcrg_evaluate(an,bn,c,u)
        ratio=real(v,8)/real(c,8)      
        y=minn+(maxx-minn)*ratio
        py=exp(-(y-mean)*(y-mean)/2.0/sigma/sigma)/sqrt(2.0*pi)/sigma
        u=v
        v=lcrg_evaluate(an,bn,c,u)
        ratio=real(v,8)/real(c,8)
        x=pymax*ratio
        if(x<=py) then
          j=j+1
          lowind= floor(y/dx)
          x=y/dx-real(lowind,8)
          array(lowind+1)=array(lowind+1)+1.0-x
          array(lowind+2)=array(lowind+2)+x          
        endif
      endif
      k=k+2*p   
      if(k.ge.k_hi.and.j.lt.pnum) then
         print*, "#ERROR: more sampling number is nesseary."
         stop 
      endif   
    end do

    call MPI_Reduce(array, buf, numcel+1, MPI_DOUBLE_PRECISION, &
         mpi_sum, 0,mpi_comm_world,error)
   
    summ(2,:)=buf(:)
    if(id==0) then
      open(10,file=trim(filepath),status='replace')
        do i=1,numcel+1
          write(10,'(2F18.10)') summ(:,i)
        end do
      close(10)
!   print*, array    !summ(2,:)
    end if


    call mpi_finalize(error)
    
    print*, "sampling is finished"
  contains

  Subroutine congruence(a,b,c,cong,error)
    integer, intent(in) :: a,b,c
    integer :: error
    integer :: cong
    
    integer :: N_MAX=100
    integer ::  a_copy
    integer ::  a_mag
    integer ::  a_sign
    integer ::  b_copy
    integer ::  b_mag
    integer ::  b_sign
    integer ::  c_copy
    integer ::  g
    integer ::  k
    integer ::  n
    real(8) ::  norm_new
    real(8) ::  norm_old
    integer, dimension(:), allocatable ::  q
    integer ::  swap
    integer ::  temp
    integer ::  x
    integer ::  xnew
    integer ::  y
    integer ::  ynew
    integer ::  z

    allocate(q(N_MAX))    
    error=0
    x=0 
    y=0
    
    if(a==0.and.b==0.and.c==0) then
      x=0
      cong=x
      goto 100
    else if(a==0.and.b==0.and.c.ne.0) then
      error=1
      x=0
      cong=x
      goto 100
    else if(a==0.and.b.ne.0.and.c==0) then
      x=0
      cong=x
      goto 100
    else if(a==0.and.b.ne.0.and.c.ne.0) then
      x=0
      if(mod(c,b).ne.0) then
        error=2
      endif
      cong=x
      goto 100
    else if(a/=0.and.b==0.and.c==0) then
      x=0
      cong=x     !!!!?????
      goto 100
    else if(a.ne.0.and.b==0.and.c.ne.0) then
      x=c/a
      if(mod(c,a).ne.0) then
        error=3
        cong=x
        goto 100
      endif
    else if(a.ne.0.and.b.ne.0.and.c==0) then
      x=0
      cong=x
      goto 100
    endif

   g=i4_gcd(a,b)

   if(mod(c,g).ne.0) then
      error=4
      cong=x
      goto 100
    endif

    a_copy = a/g
    b_copy = b/g
    c_copy = c/g

    a_mag = abs(a_copy)
    a_sign = i4_sign(a_copy)
    b_mag = abs(b_copy)
    b_sign = i4_sign(b_copy)

    if(a_mag==1) then
      x=a_sign*c_copy
      cong=x
      goto 100
    else if(b_mag==1) then
      x=0
      cong=x
      goto 100
    endif

   if(b_mag<=a_mag) then
      swap=0
      q(0)=a_mag
      q(1)=b_mag
    else
      swap=1;
      q(0)=b_mag
      q(1)=a_mag
    endif

    n=3

    do while(.true.)
      q(n-1)=mod(q(n-3),q(n-2))
      if(q(n-1)==1) then
        exit
      endif
  
      n=n+1
      if(N_MAX<n) then
        error=1
        print*
        print*, "CONGRUENCE - Fatal error" 
        print*, "Exceeded number of iterations."
        stop
      endif
    enddo

    y=0
    k=n
    do while(2.le.k)
      x=y
      y=(1-x*q(k-2))/q(k-1)
      k=k-1
    end do

    if(swap==1) then
      z=x
      x=y
      y=z
    endif

    x=x*a_sign

    x=x*c_copy

    x=mod(x,b)

    if(x<0) then
      x=x+b
    endif

    cong=x

100  end subroutine congruence


  function i4_gcd(i,j)
    integer :: i,j
    integer :: i4_gcd
    
    integer :: ip
    integer :: iq
    integer :: ir

    if(i==0) then
       i4_gcd=i4_max(1,abs(j))
    else if(j==0) then
       i4_gcd=i4_max(1,abs(i))
    endif

    ip=i4_max(abs(i),abs(j))
    iq=i4_min(abs(i),abs(j))

    do while(.true.)
      ir=mod(ip,iq)
      if(ir==0) then
        exit
      endif
      ip=iq
      iq=ir
    enddo  

    i4_gcd=iq

  end function i4_gcd


  function i4_max(i1,i2)
    integer :: i1, i2,i4_max

    integer :: value

    if(i2<i1) then
      value=i1
    else
      value=i2
    end if
    i4_max=value
    
  end function i4_max

  function i4_min(i1, i2)
    integer :: i1,i2,i4_min
    integer :: value

    if(i1<i2) then
      value=i1
    else
      value=i2
    endif
      i4_min=value

  end function i4_min


  function i4_sign(i)
    integer :: i, i4_sign

    integer :: value

    if(i<0) then
      value=-1
    else
      value=1
    endif
      i4_sign=value

  end function i4_sign

  subroutine lcrg_anbn(a,b,c,n,an,bn)
    integer :: a,b,c,n
    integer, pointer :: an,bn

    integer :: am1
    integer :: anm1tb
    integer :: ierr

    if(n<0) then
      print*
      print*, "LCRG_ANBN - Fatal error!"
      print*, "Illegal input value of N = ", n
      stop
    endif

    if(c<=0) then
      print*
      print*, "LCRG_ANBN - Fatal error!"
      print*, "Illegal input value of C = ", c
      stop
    endif

    if(n==0) then
      an=1
      bn=0
    else if(n==1) then
      an=a
      bn=b

    else

      an=power_mod(a,n,c)

      am1= a-1

      anm1tb=(an-1)*b

      call congruence(am1,c,anm1tb,bn,ierr);
 
      if( ierr .ne.  0) then
        print*,
        print*,"LCRG_ANBN - Fatal error!"
        print*,"An error occurred in the CONGRUENCE routine."
        stop
      endif
    
     end if
  end subroutine


  function lcrg_evaluate(a,b,c,x)
    integer :: a,b,c,x,lcrg_evaluate
    integer(8) :: a8
    integer(8) :: b8
    integer(8) :: c8
    integer(8) :: x8
    integer :: y
    integer(8) :: y8

    a8=int(a,8)
    b8=int(b,8)
    c8=int(c,8)
    x8=int(x,8)

    y8=mod(a8*x8+b8,c8)
    y=int(y8, 4)
    
    if(y<0) then
      y=y+c
    endif
 
    lcrg_evaluate=y
!  print*, "lcrg=",y 
  end function lcrg_evaluate


  function power_mod(a,nn,m)
    integer, intent(in) :: a, nn, m
    integer :: power_mod
    integer(8) :: a_square2
    integer :: d
    integer(8) :: m2
    integer :: x,n
    integer(8) :: x2

    n=nn
    if(a<0) then
       power_mod=-1
    endif
    if(m<=0) then
       power_mod=-1
    endif
    if(n<0) then
       power_mod=-1
    endif

    a_square2=int(a,8)
    m2=int(m,8)
    x2=int(1,8)

   do while(0<n)
      d=mod(n,2)
      if(d==1) then
        x2=mod(x2*a_square2,m2)
      endif
      a_square2=mod(a_square2*a_square2,m2)
      n=(n-d)/2
   end do

   if(x2 < 0) then
      x2=x2+m2
    endif

    x=int(x2,4)

    power_mod = x

      
  end function power_mod

  
  subroutine interp(array,v,c,maxx,dx)
    real(8),dimension(:),pointer :: array
    integer,intent(in) :: v,c
    real(8),intent(in) :: maxx,dx
    integer :: lowind
    real(8) :: x

    lowind= floor(real(v,8)/real(c,8)*maxx/dx)
    x=real(v,8)/real(c,8)*maxx/dx-real(lowind,8)
    array(lowind+1)=array(lowind+1)+1.0-x
    array(lowind+2)=array(lowind+2)+x

  end subroutine


end program para_random_congruence

