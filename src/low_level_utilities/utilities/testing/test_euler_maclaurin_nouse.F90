program   test_euler_maclaurin
#include "work_precision.h"
use utilities_module, only : &
      muarray_euler_maclaurin_choice
!      muarray_eulermaclaurin
  implicit none

  int4 :: pmax=8, d=6  ! d is the order of accuracy;pmax is the order derivative
  real8 :: mu_bound=20,dmu
  int4 ::  mu_num(4)   !
  real8, dimension(:,:,:),pointer :: muweight  ! 1st mesh of each  mu_num; 2nd,mu_num; 3nd scheme;
  real8, dimension(:,:), pointer :: mus
!  sll_real64 :: integ
  int4 :: mul,mul1,mul2,mul3
  int4 :: i,j,k,l,h

  int4 :: n=10000000
  real8 :: integ1,integ2,integ3,bessel,integ4,integ5
!  real8 :: int1,int2,int3,int4,int5
  real8 :: a(3)
  real8 :: integ(9,4)  ! 1st index for forward, central, derivative_reduction, trapzoidal
  real8 :: integanal(2)
  real8 :: dvperp,vperp,vperpmin=0.0,sum
!  sll_real64,dimension(:),allocatable ::trapeze(:),vptrapeze(:),g(:)
  int4 :: scheme(3)
  character(20) ::  number
  character(100) :: filename="/home/u2/zhang/Selalib_rep/selalib/data_multimu/",filename1
  logical :: alive 
 
  real8, dimension(:), pointer :: musbuf,muweightbuf 
  
  mu_num=(/17,33,65,129/)
  allocate(muweight(mu_num(4),4,4))
  allocate(mus(mu_num(4),4))    ! 1st mu number with redundancy 2nd: a(3)
  a=(/0.5_f64,1.0_f64,2.2_f64/)
  scheme=(/1,2,3/)  ! 1: forward 2: central  3: derivative_reduction

  do k=1,4
     integ(1,k)=mu_num(k)

  end do

  do j=1,3  ! for a
!!!!! this part is for the analytical value
  dvperp=sqrt(2.0*mu_bound)/real(n,8)
  integanal(1)=0.0_f64
  integanal(2)=0.0_f64
  integ3=0.0_f64
  integ4=0.0
  integ5=0.0
  bessel=0.0_f64

  do i=1,n
     vperp=0.0_f64+real(i-1,8)*dvperp
     bessel=bessel_j0(vperp)
     integanal(1)=integanal(1)+exp(-vperp**2/(2.0_f64*a(j)))*vperp*dvperp
     integanal(2)=integanal(2)+bessel**2*exp(-vperp**2/(2.0*a(j)))*vperp*dvperp
     bessel=bessel_j0(vperp)
     integ2=integ2+bessel*exp(-vperp**2/(2.0*a(j)))*vperp*dvperp
     integ4=integ4+vperp**2/2.0*exp(-vperp**2/(2.0*a(j)))*vperp*dvperp
     integ5=integ5+exp(-abs(vperp-1.0)-abs(vperp-1.02)-abs(vperp-1.04)-abs(vperp-1.06)-abs(vperp-1.08))*vperp*dvperp
   end do

   print*, integanal(:)

  do k=1,4 ! for mu_num
        do h=2,9
           integ(h,k)=0.0
        end do
  end do

 !!!! This part is for euler-maclaurin numuerical value
     do k=1,4  ! The four mu numbers
        allocate(musbuf(mu_num(k)))
        allocate(muweightbuf(mu_num(k)))
        dvperp=sqrt(2.0_f64*mu_bound)/real(mu_num(k)-1,8)
        do i=1,3  ! scheme
!          call muarray_euler_maclaurin_choice(mu_bound,mu_num(k),mus(:,k),muweight(:,k,i),scheme(i))
          call muarray_euler_maclaurin_choice(mu_bound,mu_num(k),musbuf,muweightbuf,scheme(i)) 
           do l=1,mu_num(k)
           vperp=real(l-1,8)*dvperp
           bessel=bessel_j0(vperp)
              integ(1+2*(i-1)+1,k)=integ(1+2*(i-1)+1,k)+muweight(l,k,i)*exp(-mus(l,k)/a(j))
    !          integ(1+2*(i-1)+2,k)=integ(1+2*(i-1)+2,k)+muweight(l,k,i)*bessel**2*exp(-vperp*vperp/(2.0*a(j)))
              integ(1+2*(i-1)+2,k)=integ(1+2*(i-1)+2,k)+muweight(l,k,i)*bessel**2*exp(-mus(l,k)/a(j))
            end do
        end do  !!! i index
!!! The trapzoidal scheme
        muweight(1,k,4)=dvperp*0.5
        muweight(mu_num(k),k,4)=dvperp*0.5
        do l=2,mu_num(k)-1
           muweight(l,k,4)=dvperp
        end do

        do l=1, mu_num(k)
           vperp=real(l-1,8)*dvperp
           bessel=bessel_j0(vperp)
           integ(8,k)=integ(8,k)+muweight(l,k,4)*exp(-vperp*vperp/(2.0*a(j)))*vperp
           integ(9,k)=integ(9,k)+muweight(l,k,4)*bessel**2*exp(-vperp*vperp/(2.0*a(j)))*vperp
        end do
        mus(:,k)=musbuf
        muweight(:,k,i)=muweightbuf
  end do  !!!! k index

    !!! compute the precision
     do k=1,4
        do i=1,4  ! scheme
           integ(1+2*(i-1)+1,k)= abs(integ(1+2*(i-1)+1,k)-integanal(1))
           integ(1+2*(i-1)+2,k)= abs(integ(1+2*(i-1)+2,k)-integanal(2))
        end do
     end do


     write(number,'(I1)') j
!     filename1=trim(filename)//"/eumaclaurin_"//trim(number)//".data"
!
!     open(j,file=filename1,status="replace")
!     do k=1,4
!       write(j,*) integ(:,k)
!     end do
  end do

end program test_euler_maclaurin
