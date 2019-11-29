module utilities_module
#include "work_precision.h"

  use constants,only: pi_
  implicit none

  public :: &
       polar_to_cartesian, &
       cartesian_to_polar, &
       derivative_polar_to_cartesian, &
       gp_error,  &
       f_is_power_of_two, &
       s_int2string,      &
       f_is_even, &
       muarray_euler_maclaurin_choice, &
       muarray_eulermaclaurin

contains
  
  subroutine polar_to_cartesian(x_polar,x_cart)
    real8,intent(in) :: x_polar(2)
    real8,intent(inout) :: x_cart(2)

    x_cart(1)=x_polar(1)*cos(x_polar(2))
    x_cart(2)=x_polar(1)*sin(x_polar(2))

  end subroutine polar_to_cartesian

  subroutine cartesian_to_polar(x_cart,x_polar)
    real8,intent(in) :: x_cart(2)
    real8,intent(inout) :: x_polar(2)

    x_polar(1)=sqrt(x_cart(1)**2+x_cart(2)**2)
    x_polar(2)=atan2(x_cart(2),x_cart(1))
  end subroutine cartesian_to_polar


  subroutine derivative_polar_to_cartesian(x_polar,deriv_polar,deriv_cart,delta)
    real8, intent(in) :: x_polar(2),deriv_polar(2),delta(2)
    real8, intent(inout) :: deriv_cart(2)
    real8 :: r,theta
    real8 :: drdx,drdy,dthdx,dthdy
    
    r=x_polar(1)
    theta=x_polar(2)
    if(x_polar(1)==0.0_f64) then
       print*, "the radial poistion is zero,r=", x_polar(1)
       stop
    end if

    if(abs(theta).le.delta(2)**2.or.abs(theta-pi_).le.delta(2)**2) then
       drdx=1._f64/cos(theta)
       dthdy=1._f64/(r*cos(theta))
       deriv_cart(1)=drdx*deriv_polar(1)
       deriv_cart(2)=dthdy*deriv_polar(2)
    else if(abs(theta-pi_/2.0_f64).le.delta(2)**2.or.abs(theta-pi_/2.0_f64).le.delta(2)**2) then
       dthdx=-1._f64/(r*sin(theta))
       drdy=1._f64/sin(theta)
       deriv_cart(1)=dthdx*deriv_polar(2)
       deriv_cart(2)=drdy*deriv_polar(1)

    else
       drdx=1._f64/cos(theta)
       dthdy=1._f64/(r*cos(theta))
       dthdx=-1._f64/(r*sin(theta))
       drdy=1._f64/sin(theta)
       deriv_cart(1)=drdx*deriv_polar(1)+dthdx*deriv_polar(2)
       deriv_cart(2)=drdy*deriv_polar(1)+dthdy*deriv_polar(2)       

    end if
  end subroutine derivative_polar_to_cartesian


  subroutine gp_error(err_code,file_name)
    integer :: err_code
!    character(len=*),, intent(in) :: descriptor
    character(len=*), intent(in) :: file_name
!    integer, optional,intent(in) :: line_number

    if(err_code .ne. 0) then
       write(*, '(a,i8)') ' Triggered in FILE '//file_name// &
            ''
       stop 'ERROR: gp_error(): exiting program'
    end if       
  
  end subroutine gp_error

  function f_is_power_of_two( n )
    int8, intent(in) :: n
    logical               :: f_is_power_of_two

    intrinsic :: not, iand

    if( (n>0) .and. (0 .eq. (iand(n,(n-1)))) ) then
       f_is_power_of_two = .true.
    else
       f_is_power_of_two = .false.
    end if
  end function f_is_power_of_two


  subroutine s_int2string( istep, cstep )
    integer         , intent(in ) :: istep   !< input integer
    character(len=*), intent(out) :: cstep   !< output string

    integer          :: l
    character(len=8) :: str_fmt

    l = len(cstep)

    if ( istep >= 0  .and. istep < 10**l) then
       str_fmt="(I0"//char(l+48)//"."//char(l+48)//")"
       write(cstep,str_fmt) istep
    else
       !      SLL_WARNING( 's_int2string', 'index is negative or too big' )
       write(*, '(a)'), 's_int2string, index is negative or too big'
       print*, 'index =', istep, ' cstep length = ', l
       cstep = 'xxxx'
    end if

  end subroutine s_int2string


  function f_is_even( n )
    int4, intent(in) :: n
    logical               :: f_is_even

    intrinsic :: modulo

    if( modulo(n,2) .eq. 0 ) then
       f_is_even = .true.
    else
       f_is_even = .false.
    end if
  end function f_is_even    

 subroutine deri_forward_coef(p,d,coef)
  integer,intent(in) :: p,d  ! p is the order of the accuracy and d denotes   the dth order of th1 derivativ2
  real(8),dimension(p+d),intent(out) :: coef
  integer :: i,j,k,imax,imin
  real(8), dimension(:,:), allocatable :: buffer,buffer1,buffer2
  real(8), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv
  real(8) :: info

  imax=p+d-1
  imin=0
  allocate(buffer(imax+1,imax+1))
  allocate(buffer1(imax+1,imax+1))
   allocate(buffer2(imax+1,imax+1))
  allocate(work(imax+1),ipiv(imax+1))

  do i=1,imax+1
     do j=1,imax+1
        buffer(i,j)=real(j-1,8)**(i-1)
     end do
  end do

  buffer1=buffer
  buffer2=0.0

  call DGETRF(imax+1,imax+1,buffer,imax+1,ipiv,info)

  if(info.ne.0) then
    stop "matrix is numerically singulari"
  end if

   call DGETRI(imax+1,buffer,imax+1,ipiv,work,imax+1,info)

   if(info.ne.0) then
    stop "matrix inversion failed"
   end if

 do i=1,imax+1
   coef(i)=buffer(i,d+1)
 end do
 deallocate(buffer)
 deallocate(buffer1)
 deallocate(buffer2)

end subroutine deri_forward_coef


 subroutine deri_central_coef(p,d,coef)
  integer,intent(in) :: p,d  ! p is the order of the accuracy and d denotes the dth order of the derivative
  real(8),dimension(:),intent(inout) :: coef
  integer :: i,j,k,imax,imin
  real(8), dimension(:,:), allocatable :: buffer
  real(8), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv
  integer :: info, im

  imax=(p+d-1)/2
  imin=-(p+d-1)/2
  if(modulo(p+d,2)==0) then
     im=p+d-1
  else
     im=p+d
  end if
  allocate(buffer(im,im))
  allocate(work(im),ipiv(im))

  do i=1,im
     do j=1,im
        buffer(i,j)=real(imin+j-1,8)**(i-1)
     end do
  end do
!  buffer1=buffer
  call DGETRF(im,im,buffer,im,ipiv,info)

  if(info.ne.0) then
    stop "matrix is numerically singulari"
  end if

 call DGETRI(im,buffer,im,ipiv,work,im,info)

   if(info.ne.0) then
    stop "matrix inversion failed"
   end if

 do i=1,im
   coef(i)=buffer(i,d+1)
 end do

end subroutine deri_central_coef

 subroutine muarray_euler_maclaurin_choice(mu_bound,mu_num,mus,muweight,scheme)
!!! It seems that this function has some problems
   integer :: pmax=8, d=6  ! d is the order of accuracy;pmax is the order derivative
  real(8),intent(in) :: mu_bound
  integer,intent(in) :: mu_num,scheme
  real(8), dimension(:),pointer, intent(inout) :: mus,muweight
  real(8), dimension(:,:), allocatable :: coefmat
  real(8), dimension(:), allocatable :: coef,bernumvec

  real(8) :: binorm, bernum
  real(8) :: integ
  real(8) :: mul,mul1,mul2,mul3
  integer :: i,j,k,l
  real(8) :: B0
  real(8) :: dvperp,vperp
  real(8),dimension(:),allocatable ::trapeze(:),vptrapeze(:),g(:)

  allocate(coefmat(pmax/2,pmax+d))
  allocate(coef(pmax+d))
  allocate(bernumvec(pmax+4))


  allocate(trapeze(0:mu_num-1),vptrapeze(0:mu_num-1),g(0:mu_num-1))

  select case (scheme)

  case (1)   ! 1: forward finite difference
    do i=1,pmax/2
     coef=0.0
     call deri_forward_coef(d,i*2-1,coef)
        coefmat(i,:)=coef(:)
    end do
  case (2)  ! 2: central finite difference
     do i=1,pmax/2
     coef=0.0
     call deri_central_coef(d,i*2-1,coef)
         coefmat(i,:)=coef(:)
     end do
  case (3)  !3: michel's scheme
     do i=1,pmax/2
       coef=0.0
       call deri_central_coef(d,2*i,coef)
       coefmat(i,:)=coef(:)
     end do
     case default
        print*, "input correct scheme, scheme1=", scheme
        stop
     end select


!  print*, coefmat(2,:)

  !!Bernu number

  B0=1.0
  do i=1,pmax+4
     if(i==1)then
        bernumvec(1)=-0.5
     else
    ! compute the binormal factor
     mul1=1
     do j=1,i+1
        mul1=mul1*j
     end do
        bernum=0.0
        do k=0,i-1
           mul2=1
           mul3=1
           if(k==0) then
              mul2=1
           else
              do j=1,k
                 mul2=mul2*j
              end do
           end if
           do j=1,i+1-k
              mul3=mul3*j
           end do
           ! compute the bernulli number
           binorm=real(mul1,8)/(real(mul2,8)*real(mul3,8))
           if(k==0) then
              bernum=bernum+binorm*B0
           else
              bernum=bernum+binorm*bernumvec(k)
           end if
        end do
         bernumvec(i)=-bernum/real(i+1,8)
     end if
!           print*, "bernum",i,bernumvec(i)
  end do

  dvperp=sqrt(2.0*mu_bound)/real(mu_num-1,8)


  select case (scheme)

     case (1)
  !!! forward difference
  do i=0,mu_num-1
     if(i==0) then
        vptrapeze(i)=0.5
        do j=1,pmax/2
              vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,1)/real(2*j,8)
        end do
     else if(i.ge.1.and.i.le.(pmax+d-2)) then
        vptrapeze(i)=1.0
        do j=1,pmax/2
           if((2*j+d-2).ge.i) then
            vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,i+1)/real(2*j,8)
           end if
         end do
      else if(i==mu_num-1) then
         vptrapeze(mu_num-1)=1.0/2.0
      else
         vptrapeze(i)=1.0    !vptrapeze(i)+1.0
      end if
   end do

    case (2)
!!!!central difference

   do i=0,mu_num-1
     if(i==0) then
        vptrapeze(i)=0.5
        do j=1,pmax/2
           vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,(2*j+d+1)/2+1)/real(2*j,8)
!           print*, bernumvec(2*j), coefmat(j,(2*j+d+1)/2+1)
        end do
     else if(i.ge.1.and.i.le.(pmax+d-1)/2) then
        vptrapeze(i)=1.0
        do j=1,pmax/2
           if((2*j+d-2)/2.ge.i) then
              vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*(-coefmat(j,(2*j+d)/2-i)+ &
                            coefmat(j,(2*j+d)/2+i))/real(2*j,8)
           end if
         end do
      else if(i==mu_num-1) then
         vptrapeze(mu_num-1)=1.0/2.0
      else
         vptrapeze(i)=1.0
      end if
   end do

    case (3)
  !!!!Michel's scheme: central difference
   do i=0,mu_num-1
     if(i==0) then
        vptrapeze(i)=dvperp*bernumvec(2)/2.0
        do j=1,pmax/2
           vptrapeze(i)=vptrapeze(i)+dvperp*bernumvec(2*(j+1))*coefmat(j,(2*j+d-1)/2+1)/real(2*(j+1),8)
        end do
     else if(i.ge.1.and.i.le.(pmax+d-1)/2) then
        vptrapeze(i)=real(i,8)*dvperp
        do j=1,pmax/2
           if((2*j+d-1)/2.ge.i) then
              vptrapeze(i)=vptrapeze(i)+dvperp*bernumvec(2*j+2)*(coefmat(j,(2*j+d-1)/2+1-i)+ &
                            coefmat(j,(2*j+d-1)/2+1+i))/real(2*j+2,8)
           end if
         end do
      else if(i==mu_num-1) then
         vptrapeze(mu_num-1)=1.0/2.0*real(i,8)*dvperp
      else
         vptrapeze(i)=real(i,8)*dvperp
      end if
   end do

   case default
      print*, "input the correct scheme,scheme2=", scheme
      stop
   end select

 !  print*,vptrapeze(:)

  select case (scheme)

  case (1)
     do i=1,mu_num
      vperp=real(i-1,8)*dvperp
      muweight(i)=vptrapeze(i-1)*vperp*dvperp
      mus(i)=vperp**2/2.0_f64
   end do

  case (2)
   do i=1,mu_num
      vperp=real(i-1,8)*dvperp
      muweight(i)=vptrapeze(i-1)*vperp*dvperp
      mus(i)=vperp**2/2.0_f64
   end do

  case (3)
    do i=1,mu_num
      vperp=real(i-1,8)*dvperp
      muweight(i)=vptrapeze(i-1)*dvperp
      mus(i)=vperp**2/2.0_f64
   end do
   case default
      print*, "input the correct scheme,scheme3=", scheme
      stop
   end select

 end subroutine muarray_euler_maclaurin_choice

subroutine muarray_eulermaclaurin(mu_bound,mu_num,mus,muweight)
  integer :: pmax=8, d=6  ! d is the order of accuracy;pmax is the order derivative
  real(8),intent(in) :: mu_bound
  integer,intent(in) :: mu_num
  real(8), dimension(:,:), allocatable :: coefmat
  real(8), dimension(:), allocatable :: coef,bernumvec
  real(8), dimension(:), pointer,intent(inout) :: mus,muweight

  real(8) :: binorm, bernum
  real(8) :: integ
  integer:: mul,mul1,mul2,mul3
  integer :: i,j,k,l

  real(8) :: B0

  real(8) :: dvperp,vperp,vperpmin=0.0,sum
  real(8), dimension(:),allocatable ::trapeze(:),vptrapeze(:),g(:)

  allocate(coefmat(pmax/2,pmax+d))
  allocate(coef(pmax+d))
  allocate(bernumvec(pmax))
  allocate(trapeze(0:mu_num-1),vptrapeze(0:mu_num-1),g(0:mu_num-1))

  do i=1,pmax/2
     coef=0.0_f64
     call deri_forward_coef(d,i*2-1,coef)
 !    call deri_central_coef(d,i*2-1,coef)
     coefmat(i,:)=coef(:)
  end do

 ! print*, coefmat(2,:)

  !!Bernu number

  B0=1.0_f64
  do i=1,pmax
     if(i==1)then
        bernumvec(1)=-0.5_f64
     else
    ! compute the binormal factor
     mul1=1
     do j=1,i+1
        mul1=mul1*j
     end do
        bernum=0.0_f64
        do k=0,i-1
           mul2=1
           mul3=1
           if(k==0) then
              mul2=1
           else
              do j=1,k
                 mul2=mul2*j
              end do
           end if
           do j=1,i+1-k
              mul3=mul3*j
           end do
           ! compute the bernulli number
           binorm=real(mul1,f64)/(real(mul2,f64)*real(mul3,f64))
           if(k==0) then
              bernum=bernum+binorm*B0
           else
 !             print*,bernum
              bernum=bernum+binorm*bernumvec(k)
 !             print*, binorm, bernumvec(1)
           end if
        end do
         bernumvec(i)=-bernum/real(i+1,8)
     end if
 !          print*, "bernum",i,bernumvec(i)
  end do

  dvperp=sqrt(2.0_f64*mu_bound)/real(mu_num-1,8)
  !!! forward difference
  do i=0,mu_num-1
     if(i==0) then
        muweight(i+1)=0.5_f64
        do j=1,pmax/2
              muweight(i+1)=muweight(i+1)+bernumvec(2*j)*coefmat(j,1)/real(2*j,8)
        end do
     else if(i.ge.1.and.i.le.(pmax+d-2)) then
        muweight(i+1)=1.0
        do j=1,pmax/2
           if((2*j+d-2).ge.i) then
            muweight(i+1)=muweight(i+1)+bernumvec(2*j)*coefmat(j,i+1)/real(2*j,8)
           end if
         end do
      else if(i==mu_num-1) then
         muweight(mu_num)=1.0/2.0
      else
         muweight(i+1)=1.0
      end if
   end do

   do i=1,mu_num
      vperp=real(i-1,8)*dvperp
      muweight(i)=muweight(i)*vperp*dvperp
      mus(i)=vperp**2/2.0_f64
   end do
   deallocate(coefmat)
   deallocate(bernumvec)
   deallocate(coef)
 end subroutine muarray_eulermaclaurin

end module utilities_module
