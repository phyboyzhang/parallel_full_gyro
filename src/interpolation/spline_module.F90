module spline_module
#include "work_precision.h"
use utilities_module,only: polar_to_cartesian, &
                     cartesian_to_polar, &
                     derivative_polar_to_cartesian, &
                     gp_error

implicit none
!include "mpif.h"

  public :: &
  s_splcoefper1d0old,&
  s_compute_splines_coefs_matrix_per_1d, &
  s_localize_per_per, &
  s_contribution_spl, &
  compute_spl2d_firstorder_derivative_point_per_per, &
  compute_spl2d_firstorder_derivative_point_nat_per, &
  compute_spl2d_firstorder_derivative_point_polar_nat_per, &    
  splcoefper1d0old, &
!  compute_splines_coefs_matrix_per_1d , &
  compute_spl2D_double_per_weight, &
  compute_spl2D_nat_per_weight, &
  s_compute_splines_coefs_matrix_nat_1d_new, &
  s_localize_nat_per_new, &
  compute_spl2d_field_point_per_per, &
  compute_spl2d_field_point_nat_per, &
  s_contribution_spl_1d, &
  compute_spl2d_field_point_nat_1d, &
  compute_spl2d_field_point_per_1d, &
  compute_D_spl2D_nat_per_noblock, &
  compute_D_spl2D_per_per_noblock, &
  s_localize_new
  
 contains
  
  subroutine s_splcoefper1d0old(dper,lper,mper,N)
    int4,intent(in)::N
    real8,dimension(0:N-1),intent(out)::dper,lper,mper
    int4::i
    
    dper(0)=4.0d0
    mper(0)=0.25d0
    do i=0,N-2
      lper(i)=1.0d0/dper(i)
      dper(i+1)=4.0d0-lper(i)
      mper(i+1)=-mper(i)/dper(i+1)
    enddo
    dper(N-1)=dper(N-1)-(lper(N-2)+2.0d0*mper(N-2))  
    do i=0,N-1
      dper(i)=1.0d0/dper(i)
    enddo
  end subroutine s_splcoefper1d0old

  subroutine s_splcoefnat1d0old_new(dper,lper,mper,N)
    int4,intent(in)::N
    real8,dimension(0:N-1),intent(out)::dper,lper,mper
    int4::i
    
    dper(0)=6._f64
    mper(0)=0._f64
    do i=0,N-2
      lper(i)=1._f64/dper(i)
      dper(i+1)=4._f64-lper(i)
      mper(i+1)=-mper(i)/dper(i+1)
    enddo
    dper(N-1)=dper(N-1)-(lper(N-2)+2._f64*mper(N-2))  
    do i=0,N-1
      dper(i)=1._f64/dper(i)
    enddo
  end subroutine s_splcoefnat1d0old_new  


  subroutine s_compute_splines_coefs_matrix_nat_1d_new(mat,N)
    int4,intent(in)::N
    real8,dimension(:,:),intent(inout) :: mat
    real8,dimension(:),allocatable :: dnat,a,b,d
    real8,dimension(:,:), allocatable :: v,u
    int4 :: i,j,l
    int4 :: low,up

    low=lbound(mat,1)
    up=ubound(mat,1)
 
    allocate(v(low:up,low:up))
    allocate(u(low:up,low:up))
    allocate(dnat(low:up),a(low:up),b(low:up),d(low:up))

    d(low)=6.0_f64
    a(low)=0.0_f64
    b(low)=0.0_f64
    do i=low+1,up-1
       d(i)=4.0_f64
       a(i)=1.0_f64
       b(i)=1.0_f64
    end do
    d(up)=6.0_f64
    a(up)=0.0_f64
    b(up)=0.0_f64    
   
    d=d/6.0_f64
    a=a/6.0_f64
    b=b/6.0_f64
    
    dnat(low)=d(low)
    do i=low+1,up
       dnat(i)=d(i)-b(i)*a(i-1)/dnat(i-1)
    end do

    do i=low,up
       v(i,i)=1.0_f64
       u(i,i)=1.0_f64/dnat(i)
    end do
    
    do j=low+1,up
       do i=j-1,1,-1
          v(i,j)=0._f64
          u(i,j)=-a(i)*u(i+1,j)/dnat(i)
       end do       
    end do

    do i=low+1,up
       do j=1,i-1
          u(i,j)=0._f64
          v(i,j)=-b(i)*v(i-1,j)/dnat(i-1)
       end do
    end do

    mat=0.0
    do i=low,up
       do j=low,up
          do l=max(i,j),up
             mat(i,j)=mat(i,j)+u(i,l)*v(l,j)
          end do
       end do
    end do

    deallocate(u)
    deallocate(v)
  end subroutine s_compute_splines_coefs_matrix_nat_1d_new    

  subroutine s_compute_splines_coefs_matrix_per_1d(mat,dper,lper,mper,N)
    int4,intent(in)::N
    real8,dimension(0:N-1,0:N-1),intent(inout) :: mat
    real8,dimension(0:N+2),intent(in)::dper,lper,mper
    int8 :: i,j

    mat = 0.0d0
    
    do j = 0,N-1
      mat(j,j) = 6.0d0
    enddo

    do i = 1,N-1
      do j = 0,N-1
      mat(i,j) = mat(i,j)-lper(i-1)*mat(i-1,j)
      enddo
    enddo

    do i = 0,N-2
      do j = 0,N-1
      mat(N-1,j) = mat(N-1,j)-mper(i)*mat(i,j)
      enddo
    enddo
    
    do j = 0,N-1
      mat(N-1,j) = mat(N-1,j)*dper(N-1) 
    enddo
    
    do j = 0,N-1
      mat(N-2,j) = dper(N-2)*(mat(N-2,j)-(1.0d0-mper(N-3))*mat(N-1,j))
    enddo

    do i=N-3,1,-1
      do j = 0,N-1
        mat(i,j) = dper(i)*(mat(i,j)-mat(i+1,j)+mper(i-1)*mat(N-1,j))
      enddo
    enddo
        
    do j = 0,N-1
      mat(0,j) = dper(0)*(mat(0,j)-mat(1,j)-mat(N-1,j))
    enddo
  
  end subroutine s_compute_splines_coefs_matrix_per_1d
  

  subroutine s_localize_per_per(x,eta_min,eta_max,ii,eta,N,flag)
    real8,intent(in)::x(2),eta_min(2),eta_max(2)
    int4,intent(out)::ii(2)
    int4,intent(in)::N(2)
    real8,intent(out)::eta(2)
    int4,optional :: flag
    eta=x
    call localize_per(ii(1),eta(1),eta_min(1),eta_max(1),N(1))
    call localize_per(ii(2),eta(2),eta_min(2),eta_max(2),N(2))

  end subroutine s_localize_per_per

  subroutine s_localize_new(x,eta_min,eta_max,ii,eta,N,flag)
    real8,intent(in)::x(2),eta_min(2),eta_max(2)
    int4,intent(out)::ii(2)
    int4,intent(in)::N(2)
    real8,intent(out)::eta(2)
    int4,optional :: flag
    eta=x
    call localize_new(ii(1),eta(1),eta_min(1),eta_max(1),N(1))
    call localize_new(ii(2),eta(2),eta_min(2),eta_max(2),N(2))

  end subroutine s_localize_new


  subroutine localize_per(i,x,xmin,xmax,N)
    int4,intent(out)::i
    real8,intent(inout)::x
    real8,intent(in)::xmin,xmax
    int4,intent(in)::N
    x=(x-xmin)/(xmax-xmin)
    x=x-real(floor(x),8)
    x=x*real(N,8)
    i=floor(x)
    x=x-real(i,8)
    if(i==N)then
      i=0
      x=0.0d0
    endif
  end subroutine localize_per

 subroutine localize_new(i,x,xmin,xmax,N) !!!Here, x is within xmin and xmax
    int4,intent(out)::i
    real8,intent(inout)::x
    real8,intent(in)::xmin,xmax
    int4,intent(in)::N
    real8 :: xbar, delta
    x=(x-xmin)/(xmax-xmin)
    x=x*real(N,8)
    i=floor(x)
    x=x-real(i,8)
    i=i+1
  end subroutine localize_new

  
  subroutine localize_nat(i,x,xmin,xmax,N,flag)  ! N is the number of cells
    int4,intent(out)::i
    real8,intent(inout)::x
    real8,intent(in)::xmin,xmax
    int4,intent(in)::N
    int4,optional :: flag
    x=(x-xmin)/(xmax-xmin)
    x=x*real(N,8)
    if(x>=real(N,8))then
      x=real(N,8)
    endif
    if(x<0.0d0)then
       x=0.0d0
       flag=0
    endif
 !   if(x<0._f64)then  !! with this judgement, the points out of the boudnary are ignored.
 !      flag=0
 !   end if
    i=floor(x)   !! here indicates when i==0, we choose i=0
    x=x-real(i,8)
    if(i==N)then
      i=N-1
      x=1.0d0
   endif
   
 end subroutine localize_nat

  subroutine s_localize_nat_per_new(x,eta_min,eta_max,ii,eta,N,flag)
    real8,intent(in)::x(2),eta_min(2),eta_max(2)
    int4,intent(out)::ii(2)
    int4,intent(in)::N(2)
    real8,intent(out)::eta(2)
    int4,optional :: flag
    eta=x
    call localize_nat(ii(1),eta(1),eta_min(1),eta_max(1),N(1),flag)

    call localize_per(ii(2),eta(2),eta_min(2),eta_max(2),N(2))
  end subroutine s_localize_nat_per_new
 
   subroutine s_contribution_spl(x,val)
    real8,intent(in)::x(0:1)
    real8,intent(out)::val(-1:2,-1:2)
    int4::s
    real8::w(-1:2,0:1)
    do s=0,1
      w(-1,s)=(1.0_f64/6._f64)*(1._f64-x(s))*(1._f64-x(s))*(1._f64-x(s));
      w(0,s)=1._f64/6._f64+0.5_f64*(1._f64-x(s))*(-(1._f64-x(s))*&
       (1._f64-x(s))+(1._f64-x(s))+1._f64);
      w(1,s)=1._f64/6._f64+0.5_f64*x(s)*(-x(s)*x(s)+x(s)+1._f64);
      w(2,s)=(1._f64/6._f64)*x(s)*x(s)*x(s);
    enddo
    
    
    val(-1,-1)=w(-1,0)*w(-1,1)  
    val(-1,0)=w(-1,0)*w(0,1)  
    val(-1,1)=w(-1,0)*w(1,1)  
    val(-1,2)=w(-1,0)*w(2,1)  

    val(0,-1)=w(0,0)*w(-1,1)  
    val(0,0)=w(0,0)*w(0,1)  
    val(0,1)=w(0,0)*w(1,1)  
    val(0,2)=w(0,0)*w(2,1)  
    
    val(1,-1)=w(1,0)*w(-1,1)  
    val(1,0)=w(1,0)*w(0,1)  
    val(1,1)=w(1,0)*w(1,1)  
    val(1,2)=w(1,0)*w(2,1)  
    
    val(2,-1)=w(2,0)*w(-1,1)  
    val(2,0)=w(2,0)*w(0,1)  
    val(2,1)=w(2,0)*w(1,1)  
    val(2,2)=w(2,0)*w(2,1)  

  end subroutine s_contribution_spl

  subroutine s_contribution_spl_1d(x,val)
    real8,intent(in) :: x
    real8,intent(out) :: val(-1:2)
      val(-1)=(1.0_f64/6._f64)*(1._f64-x)*(1._f64-x)*(1._f64-x)
      val(0)=1._f64/6._f64+0.5_f64*(1._f64-x)*(-(1._f64-x)*&
       (1._f64-x)+(1._f64-x)+1._f64)
      
      val(1)=1._f64/6._f64+0.5_f64*x*(-x*x+x+1._f64)
      val(2)=(1._f64/6._f64)*x*x*x
  end subroutine
 
  subroutine compute_spl_derivative_spline_1d(x,val_deriv)
    real8, intent(in) :: x
    real8, intent(inout) :: val_deriv(-1:2)
  end subroutine compute_spl_derivative_spline_1d

  subroutine compute_spl_derivative_spline_2d(x, val_first_deriv,val_second_deriv,eta_delta)
    real8,intent(in)::x(0:1)
!    real8,intent(out)::val(-1:2,-1:2)
    real8,intent(inout) :: val_first_deriv(-1:2,-1:2)
    real8,intent(inout) :: val_second_deriv(-1:2,-1:2)    
 !   real8,intent(out) :: val_deriv_two(-1:2,-1:2)
    real8, dimension(0:1), intent(in) :: eta_delta   
    int4 :: s, j
    real8 :: w(-1:2,0:1)
    real8 :: w_deriv(-1:2,0:1)
    do s=0,1
      w(-1,s)=(1.0_f64/6._f64)*(1._f64-x(s))*(1._f64-x(s))*(1._f64-x(s))
      w(0,s)=1._f64/6._f64+0.5_f64*(1._f64-x(s))*(-(1._f64-x(s))*   &
       (1._f64-x(s))+(1._f64-x(s))+1._f64)
      w(1,s)=1._f64/6._f64+0.5_f64*x(s)*(-x(s)*x(s)+x(s)+1._f64)
      w(2,s)=(1._f64/6._f64)*x(s)*x(s)*x(s)
    enddo

   do s=0,1
      w_deriv(-1,s)=-1.0_f64/(2._f64*eta_delta(s))*(1._f64-x(s))**2
      w_deriv(0,s)=-1._f64/(2._f64*eta_delta(s))-1._f64/eta_delta(s)*(1._f64-x(s)) &
           +3._f64/(2._f64*eta_delta(s))*(1._f64-x(s))**2
      w_deriv(1,s)=1._f64/(2._f64*eta_delta(s))+1._f64*x(s)/eta_delta(s)-3._f64*x(s)**2/(2._f64*eta_delta(s))
      w_deriv(2,s)=1._f64/(2._f64*eta_delta(s))*x(s)**2
   enddo

!   w_deriv=-1.0_f64*w_deri

!!!===> The imterpolation of the first coordinate at x
   
    val_first_deriv(-1,-1)=w_deriv(-1,0)*w(-1,1)  
    val_first_deriv(-1,0)=w_deriv(-1,0)*w(0,1)  
    val_first_deriv(-1,1)=w_deriv(-1,0)*w(1,1)  
    val_first_deriv(-1,2)=w_deriv(-1,0)*w(2,1)  

    val_first_deriv(0,-1)=w_deriv(0,0)*w(-1,1)  
    val_first_deriv(0,0)=w_deriv(0,0)*w(0,1)  
    val_first_deriv(0,1)=w_deriv(0,0)*w(1,1)  
    val_first_deriv(0,2)=w_deriv(0,0)*w(2,1)  
    
    val_first_deriv(1,-1)=w_deriv(1,0)*w(-1,1)  
    val_first_deriv(1,0)=w_deriv(1,0)*w(0,1)  
    val_first_deriv(1,1)=w_deriv(1,0)*w(1,1)  
    val_first_deriv(1,2)=w_deriv(1,0)*w(2,1)  
    
    val_first_deriv(2,-1)=w_deriv(2,0)*w(-1,1)  
    val_first_deriv(2,0)=w_deriv(2,0)*w(0,1)  
    val_first_deriv(2,1)=w_deriv(2,0)*w(1,1)  
    val_first_deriv(2,2)=w_deriv(2,0)*w(2,1)

!!!===> The imterpolation of the second coordinate at x  
    
    val_second_deriv(-1,-1)=w(-1,0)*w_deriv(-1,1)  
    val_second_deriv(-1,0)=w(-1,0)*w_deriv(0,1)  
    val_second_deriv(-1,1)=w(-1,0)*w_deriv(1,1)  
    val_second_deriv(-1,2)=w(-1,0)*w_deriv(2,1)  

    val_second_deriv(0,-1)=w(0,0)*w_deriv(-1,1)  
    val_second_deriv(0,0)=w(0,0)*w_deriv(0,1)  
    val_second_deriv(0,1)=w(0,0)*w_deriv(1,1)  
    val_second_deriv(0,2)=w(0,0)*w_deriv(2,1)  
    
    val_second_deriv(1,-1)=w(1,0)*w_deriv(-1,1)  
    val_second_deriv(1,0)=w(1,0)*w_deriv(0,1)  
    val_second_deriv(1,1)=w(1,0)*w_deriv(1,1)  
    val_second_deriv(1,2)=w(1,0)*w_deriv(2,1)  
    
    val_second_deriv(2,-1)=w(2,0)*w_deriv(-1,1)  
    val_second_deriv(2,0)=w(2,0)*w_deriv(0,1)  
    val_second_deriv(2,1)=w(2,0)*w_deriv(1,1)  
    val_second_deriv(2,2)=w(2,0)*w_deriv(2,1)    
   

  end subroutine compute_spl_derivative_spline_2d


  subroutine compute_spl2d_firstorder_derivative_point_per_per(x, &
       ep_weight,  &
       eta_min,    &
       eta_max,    &
       eta_delta,  &
       NC,         &
       deri_firstorder)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    real8,dimension(:,:),pointer,intent(in) :: ep_weight
    real8,dimension(:), intent(inout) :: deri_firstorder
    real8,dimension(:), intent(in) :: eta_min(2),eta_max(2), eta_delta(2)
    int4, intent(in) :: NC(2)
    real8 :: eta_star(2)
    real8 :: val_1st_deri(-1:2,-1:2), val_2nd_deri(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2   
    int4 :: flag
    
    deri_firstorder(1)=0.0_f64
    deri_firstorder(2)=0.0_f64
    call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,NC,flag)
    call compute_spl_derivative_spline_2d(eta_star,val_1st_deri,val_2nd_deri,eta_delta)
    
    do ell_2=-1,2
       ind(2)=modulo(ii(2)+ell_2,NC(2))
       do ell_1=-1,2
          ind(1)=modulo(ii(1)+ell_1,NC(1))         
          deri_firstorder(1)=deri_firstorder(1)+val_1st_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
          deri_firstorder(2)=deri_firstorder(2)+val_2nd_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
       end do
    end do
!print*, "per_per,deri_", deri_firstorder(1:2) 
  end subroutine compute_spl2d_firstorder_derivative_point_per_per
  
   subroutine compute_spl2d_firstorder_derivative_point_nat_per(x, & 
       ep_weight,  &
       eta_min,    & 
       eta_max,    &
       eta_delta,  &
       NC,         &
       deri_firstorder)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    real8,dimension(:,:),pointer,intent(in) :: ep_weight
    real8,dimension(:), intent(inout) :: deri_firstorder
    real8,dimension(:), intent(in) :: eta_min(2),eta_max(2), eta_delta(2)
    int4, intent(in) :: NC(2)
    real8 :: eta_star(2)
    real8 :: val_1st_deri(-1:2,-1:2), val_2nd_deri(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2   
    int4 :: flag
    
    deri_firstorder(1)=0.0_f64
    deri_firstorder(2)=0.0_f64
    call s_localize_nat_per_new(x,eta_min,eta_max,ii,eta_star,NC,flag)
    call compute_spl_derivative_spline_2d(eta_star,val_1st_deri,val_2nd_deri,eta_delta)
    do ell_2=-1,2
       ind(2)=modulo(ii(2)+ell_2,NC(2))
       do ell_1=-1,2
          ind(1)=ii(1)+ell_1
        if(ind(1)<0) then
          ind(1)=0
        else if(ind(1).ge.(Nc(1)-1)) then
          ind(1)=Nc(1)-1
        end if          
          deri_firstorder(1)=deri_firstorder(1)+val_1st_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
          deri_firstorder(2)=deri_firstorder(2)+val_2nd_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
       end do
    end do

 !  print*, "nat_per, deri_firstorder=",deri_firstorder 
  end subroutine compute_spl2d_firstorder_derivative_point_nat_per 


   subroutine compute_spl2d_firstorder_derivative_point_polar_nat_per(x_polar, & 
       ep_weight,  &
       eta_min,    & 
       eta_max,    &
       eta_delta,  &
       NC,         &
       deriv_cart)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x_polar(2)
    real8,dimension(:,:),pointer,intent(in) :: ep_weight
    real8,dimension(:), intent(inout) :: deriv_cart
    real8,dimension(:), intent(in) :: eta_min(2),eta_max(2), eta_delta(2)
    int4, intent(in) :: NC(2)
    real8 :: eta_star(2)
    real8 :: val_1st_deri(-1:2,-1:2), val_2nd_deri(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2   
    int4 :: flag
    real8 :: deriv_polar(2)
    
    deriv_polar(1)=0.0_f64
    deriv_polar(2)=0.0_f64
    call s_localize_nat_per_new(x_polar,eta_min,eta_max,ii,eta_star,NC,flag)
    call compute_spl_derivative_spline_2d(eta_star,val_1st_deri,val_2nd_deri,eta_delta)
    
    do ell_2=-1,2
       ind(2)=modulo(ii(2)+ell_2,NC(2))
       do ell_1=-1,2
         ind(1)=ii(1)+ell_1
       if(ind(1)<0) then
        ind(1)=0
       else if(ind(1).ge.(Nc(1)-1)) then
          ind(1)=Nc(1)-1
       end if          
          deriv_polar(1)=deriv_polar(1)+val_1st_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
          deriv_polar(2)=deriv_polar(2)+val_2nd_deri(ell_1,ell_2)*ep_weight(ind(1)+1,ind(2)+1)
       end do
    end do

    call derivative_polar_to_cartesian(x_polar,deriv_polar,deriv_cart,eta_delta)
 
  end subroutine compute_spl2d_firstorder_derivative_point_polar_nat_per   

   function compute_spl2d_field_point_per_per(x, &
       field_weight,  &
       eta_min,    &
       eta_max,    &
       eta_delta,  &
       NC) result(field_value)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    real8,dimension(:,:),pointer,intent(in) :: field_weight
    real8,dimension(:), intent(in) :: eta_min(2),eta_max(2), eta_delta(2)
    int4, intent(in) :: NC(2)
    real8 :: field_value
    real8 :: eta_star(2)
    real8 :: val(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2   
    int4 :: flag
    
    call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,NC,flag)
    call  s_contribution_spl(eta_star, val)

    field_value = 0._f64
    do ell_2=-1,2
       ind(2)=modulo(ii(2)+ell_2,NC(2))
       do ell_1=-1,2
          ind(1)=modulo(ii(1)+ell_1,NC(1))         
          field_value=field_value+val(ell_1,ell_2)*field_weight(ind(1)+1,ind(2)+1)         
       end do
    end do
 
  end function compute_spl2d_field_point_per_per

  function compute_spl2d_field_point_per_1d(x, &
       field_weight,  &
       eta_min,    &
       eta_max,    &
       eta_delta,  &
       NC) result(field_value)
    real8,intent(inout) :: x
    real8,dimension(:),pointer,intent(in) :: field_weight
    real8,intent(in) :: eta_min,eta_max, eta_delta
    int4, intent(in) :: NC
    real8 :: field_value
    real8 :: val(-1:2)
    int4 :: ii,ind,ell
    int4 :: flag

    call localize_per(ii,x,eta_min,eta_max,NC)
    call s_contribution_spl_1d(x,val)
    
    field_value=0._f64
    do ell=-1,2
      ind=modulo(ii+ell,NC)
      field_value=field_value+val(ell)*field_weight(ind+1)
    end do

  end function



   function compute_spl2d_field_point_nat_per(x, &
       field_weight,  &
       eta_min,    &
       eta_max,    &
       eta_delta,  &
       NC) result(field_value)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    real8,dimension(:,:),pointer,intent(in) :: field_weight
    real8,dimension(:), intent(in) :: eta_min(2),eta_max(2), eta_delta(2)
    int4, intent(in) :: NC(2)
    real8 :: field_value
    real8 :: eta_star(2)
    real8 :: val(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2   
    int4 :: flag
    
    call s_localize_nat_per_new(x,eta_min,eta_max,ii,eta_star,NC,flag)
    call  s_contribution_spl(eta_star, val)

    field_value = 0._f64
    do ell_2=-1,2
       ind(2)=modulo(ii(2)+ell_2,NC(2))
       do ell_1=-1,2           
          ind(1)=ii(1)+ell_1
                
          if(ind(1)<0) then
             ind(1)=0
          else if(ind(1).ge.(Nc(1)-1)) then
             ind(1)=Nc(1)-1
          end if        
          field_value=field_value+val(ell_1,ell_2)*field_weight(ind(1)+1,ind(2)+1)         
       end do
    end do
 
  end function compute_spl2d_field_point_nat_per 

  function compute_spl2d_field_point_nat_1d(x, &
       field_weight,  &
       eta_min,    &
       eta_max,    &
       eta_delta,  &
       NC) result(field_value)
    
    real8 :: x
    real8,dimension(:),pointer,intent(in) :: field_weight
    real8,intent(in) :: eta_min,eta_max, eta_delta
    int4, intent(in) :: NC
    real8 :: field_value
    real8 :: val(-1:2)
    int4 :: ii,ind,ell
    int4 :: flag

    call localize_nat(ii,x,eta_min,eta_max,NC)
    call s_contribution_spl_1d(x,val)

    field_value=0._f64
    do ell=-1,2
      ind=ii+ell
      if(ind<0) then
             ind=0
      else if(ind.ge.(NC-1)) then
             ind=NC-1
      end if            
      field_value=field_value+val(ell)*field_weight(ind+1)
    end do

  end function

  subroutine splcoefper1d0old(dper,lper,mper,N)
    int4,intent(in)::N
    real8,dimension(0:N-1),intent(out)::dper,lper,mper
    int4::i
    
    dper(0)=4._f64
    mper(0)=0.25_f64
    do i=0,N-2
      lper(i)=1._f64/dper(i)
      dper(i+1)=4._f64-lper(i)
      mper(i+1)=-mper(i)/dper(i+1)
    enddo
    dper(N-1)=dper(N-1)-(lper(N-2)+2._f64*mper(N-2))  
    do i=0,N-1
      dper(i)=1._f64/dper(i)
    enddo
  end subroutine splcoefper1d0old


!!$subroutine compute_splines_coefs_matrix_per_1d(mat,dper,lper,mper,N)
!!$    int4,intent(in)::N
!!$    int4,dimension(0:N-1,0:N-1),intent(inout) :: mat
!!$    real8,dimension(0:N+2),intent(in)::dper,lper,mper
!!$    int4 :: i,j
!!$
!!$    mat = 0._f64
!!$    
!!$    do j = 0,N-1
!!$      mat(j,j) = 6._f64
!!$    enddo
!!$
!!$    do i = 1,N-1
!!$      do j = 0,N-1
!!$      mat(i,j) = mat(i,j)-lper(i-1)*mat(i-1,j)
!!$      enddo
!!$    enddo
!!$
!!$    do i = 0,N-2
!!$      do j = 0,N-1
!!$      mat(N-1,j) = mat(N-1,j)-mper(i)*mat(i,j)
!!$      enddo
!!$    enddo
!!$    
!!$    do j = 0,N-1
!!$      mat(N-1,j) = mat(N-1,j)*dper(N-1) 
!!$    enddo
!!$    
!!$    do j = 0,N-1
!!$      mat(N-2,j) = dper(N-2)*(mat(N-2,j)-(1._f64-mper(N-3))*mat(N-1,j))
!!$    enddo
!!$
!!$    do i=N-3,1,-1
!!$      do j = 0,N-1
!!$        mat(i,j) = dper(i)*(mat(i,j)-mat(i+1,j)+mper(i-1)*mat(N-1,j))
!!$      enddo
!!$    enddo
!!$        
!!$    do j = 0,N-1
!!$      mat(0,j) = dper(0)*(mat(0,j)-mat(1,j)-mat(N-1,j))
!!$    enddo
!!$  
!!$  end subroutine compute_splines_coefs_matrix_per_1d  
  

  subroutine compute_spl2D_double_per_weight( &   !==> compute the spline weight at the mesh points
    field_init, &
    mesh_field_weight, &
    num_cells_x, &
    num_cells_y )
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:),pointer, intent(in) :: field_init
    real8, dimension(:,:), pointer, intent(inout) :: mesh_field_weight
    real8, dimension(:),allocatable,target::dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
    real8, dimension(:),pointer::pointer_dper_x,pointer_lper_x,pointer_mper_x
    real8, dimension(:),pointer::pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8, dimension(:,:),allocatable,target:: mat_per_x,mat_per_y
    real8, dimension(:,:),pointer::pointer_mat_per_x,pointer_mat_per_y
    int4  :: j
    int4  :: k,i,l

    real8, dimension(:,:),allocatable :: buffer
    real8, dimension(:), allocatable :: buf1_field,buf2_field
    
    Nx = num_cells_x
    Ny = num_cells_y


    ALLOCATE(dper_x(0:Nx-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(lper_x(0:Nx-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mper_x(0:Nx-1),stat=ierr) 
call gp_error(ierr,"spline_module")   
    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
    ALLOCATE(mat_per_x(0:Nx-1,0:Nx-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr) 
call gp_error(ierr,"spline_module")   
    allocate(buffer(Ny*Nx,Ny*Nx),stat=ierr)
    call gp_error(ierr,"spline_module")
    allocate(buf1_field(Ny*Nx),stat=ierr)
    call gp_error(ierr,"spline_module")
    allocate(buf2_field(Ny*Nx),stat=ierr)
    call gp_error(ierr,"spline_module")


    pointer_dper_y => dper_y
    pointer_lper_y => lper_y
    pointer_mper_y => mper_y    
    pointer_mat_per_y => mat_per_y
 !   pointer_mat_spl2D_circ => mat_spl2D_circ

    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)
    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
    

    pointer_dper_x => dper_x
    pointer_lper_x => lper_x
    pointer_mper_x => mper_x    
    pointer_mat_per_x=> mat_per_x
   
    call s_splcoefper1d0old(dper_x,lper_x,mper_x,Nx)
    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_x,pointer_dper_x,pointer_lper_x,pointer_mper_x,Nx)

    do j=1,Ny
       do k=1,Nx
          do i=1,Ny
             do l=1,Nx
                buffer((j-1)*Nx+k,(i-1)*Nx+l)=mat_per_y(j-1,i-1)*mat_per_x(k-1,l-1)
             end do
          end do
       end do
    enddo

    do i=1, Ny
       do j=1,Nx
          buf1_field((i-1)*Nx+j)=field_init(j,i)
       end do
    end do

    buf2_field=0.0_f64
    do i=1,Nx*Ny
       do j=1,Nx*Ny
          buf2_field(i)=buf2_field(i)+buffer(i,j)*buf1_field(j)
       end do
    end do
    

    do i=1, Ny
       do j=1,Nx
          mesh_field_weight(j,i)=buf2_field((i-1)*Nx+j)
       end do
    end do   

 
    deallocate(buffer)
    deallocate(buf1_field)
    deallocate(buf2_field)
    deallocate(dper_x)
    deallocate(lper_x)
    deallocate(mper_x)
    deallocate(dper_y)
    deallocate(lper_y)
    deallocate(mper_y)
    deallocate(mat_per_x)
    deallocate(mat_per_y)
!  print*, "#mat_spl2d", mat_spl2d(:,:,1)
  end subroutine compute_spl2D_double_per_weight



  subroutine compute_spl2D_nat_per_weight( &   !==> compute the spline weight at the mesh points
    field_init, &
    mesh_field_weight, &
    num_cells_x, &
    num_cells_y )
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:), pointer,intent(in) :: field_init
    real8, dimension(:,:), pointer,intent(inout) :: mesh_field_weight
    real8, dimension(:),allocatable,target::dnat,lnat,dper_y,lper_y,mper_y
    real8, dimension(:),pointer::pointer_dnat,pointer_lnat
    real8, dimension(:),pointer::pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8, dimension(:,:),allocatable, target :: mat_nat
    real8, dimension(:,:),pointer :: pointer_mat_nat
    real8, dimension(:,:),allocatable,target:: mat_per_y
    real8, dimension(:,:),pointer::pointer_mat_per_x,pointer_mat_per_y
    int4  :: j
    int4  :: k,i,l
 !   comp8, dimension(:), allocatable :: fft_array
 !   real8, dimension(:), allocatable :: buf_fft
    real8, dimension(:,:),allocatable :: buffer
    real8, dimension(:), allocatable :: buf1_field,buf2_field
    
    Nx = num_cells_x
    Ny = num_cells_y


    ALLOCATE(dnat(0:Nx),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(lnat(0:Nx),stat=ierr) 
call gp_error(ierr,"spline_module")  
    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mat_nat(0:Nx,0:Nx),stat=ierr)
call gp_error(ierr,"spline_module")
    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr)
call gp_error(ierr,"spline_module")    
!    ALLOCATE(buf_fft(4*Ny+100),stat=ierr)
!    ALLOCATE(fft_array(Ny),stat=ierr)
    allocate(buffer(Ny*(Nx+1),Ny*(Nx+1)),stat=ierr)
    call gp_error(ierr,"buffer")
    allocate(buf1_field(Ny*(Nx+1)),stat=ierr)
    call gp_error(ierr,"spline_module")
    allocate(buf2_field(Ny*(Nx+1)),stat=ierr)
    call gp_error(ierr,"spline_module")
    pointer_dnat => dnat
    pointer_lnat => lnat
    pointer_mat_nat => mat_nat
    pointer_dper_y => dper_y
    pointer_lper_y => lper_y
    pointer_mper_y => mper_y    
    pointer_mat_per_y => mat_per_y
 !   pointer_mat_spl2D_circ => mat_spl2D_circ

    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)

    call s_compute_splines_coefs_matrix_nat_1d_new(pointer_mat_nat,Nx+1)    
    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)

    do j=1,Ny
       do k=1,Nx+1
          do i=1,Ny
             do l=1,Nx+1
                buffer((j-1)*(Nx+1)+k,(i-1)*(Nx+1)+l)=mat_per_y(j-1,i-1)*mat_nat(k-1,l-1)
             end do
          end do
       end do
    enddo

    do i=1, Ny
       do j=1,Nx+1
          buf1_field((i-1)*(Nx+1)+j)=field_init(j,i)
       end do
    end do

    buf2_field=0.0_f64
    do i=1,(Nx+1)*Ny
       do j=1,(Nx+1)*Ny
          buf2_field(i)=buf2_field(i)+buffer(i,j)*buf1_field(j)
       end do
    end do
    

    do i=1, Ny
       do j=1,Nx+1
          mesh_field_weight(j,i)=buf2_field((i-1)*(Nx+1)+j)
       end do
    end do   

 
    deallocate(buffer)
    deallocate(buf1_field)
    deallocate(buf2_field)
    deallocate(lnat)
    deallocate(dnat)
    deallocate(dper_y)
    deallocate(lper_y)
    deallocate(mper_y)
    deallocate(mat_nat)
    deallocate(mat_per_y)
!  print*, "#mat_spl2d", mat_spl2d(:,:,1)
  end subroutine compute_spl2D_nat_per_weight


subroutine compute_D_spl2D_per_per_noblock( &
    num_cells_x, &
    num_cells_y, &
    mat_spl2D)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:), intent(inout) :: mat_spl2D
!    real8, dimension(:,:,:),allocatable :: buf_3d
    real8, dimension(:),allocatable,target :: dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
    real8, dimension(:),pointer :: pointer_dper_x,pointer_lper_x,pointer_mper_x
    real8, dimension(:),pointer :: pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8, dimension(:,:),allocatable,target :: mat_per_x,mat_per_y
    real8, dimension(:,:),pointer :: pointer_mat_per_x,pointer_mat_per_y
    real8, dimension(:,:),pointer :: buf_x,buf_y
    int4 :: k,i,h,j
    real8 :: testm(num_cells_x,num_cells_y),test(4,4)
        
    Nx = num_cells_x
    Ny = num_cells_y


    ALLOCATE(dper_x(0:Nx-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(lper_x(0:Nx-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mper_x(0:Nx-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")    
    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mat_per_x(0:Nx-1,0:Nx-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(buf_3d(1:size(mat_spl2d,3),1:size(mat_spl2d,1),1:size(mat_spl2d,2)),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(buf_x(0:Nx,0:Nx),stat=ierr)
!    ALLOCATE(buf_y(0:Ny,0:Ny),stat=ierr)


    dper_x=0.0d0
    lper_x=0.0d0
    mper_x=0.0d0
    dper_y=0.0d0
    lper_y=0.0d0
    mper_y=0.0d0
    mat_per_x=0.0d0
    mat_per_y=0.0d0

    pointer_dper_y => dper_y
    pointer_lper_y => lper_y
    pointer_mper_y => mper_y    
    pointer_mat_per_y => mat_per_y
 !   pointer_mat_spl2D_circ => mat_spl2D_circ

    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)
    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
    

    pointer_dper_x => dper_x
    pointer_lper_x => lper_x
    pointer_mper_x => mper_x    
    pointer_mat_per_x => mat_per_x
    
    call s_splcoefper1d0old(dper_x,lper_x,mper_x,Nx)
    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_x,pointer_dper_x,pointer_lper_x,pointer_mper_x,Nx)

 do h=0,Ny-1
    do j=0,Ny-1
       do i=0,Nx-1
          do k=0,Nx-1
             mat_spl2d(h*(Nx)+i+1,j*(Nx)+k+1)=mat_per_y(h,j)*mat_per_x(i,k)
          end do
       end do
    enddo
 end do   

    deallocate(dper_x)
    deallocate(lper_x)
    deallocate(mper_x)
    deallocate(dper_y)
    deallocate(mper_y)
    deallocate(lper_y)
    deallocate(mat_per_x)
    deallocate(mat_per_y)
!    deallocate(buf_x)
!    deallocate(buf_y)
    
  end subroutine compute_D_spl2D_per_per_noblock



 subroutine compute_D_spl2D_nat_per_noblock( &
    num_cells_x, &
    num_cells_y, &
    mat_spl2D)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:), intent(inout) :: mat_spl2D
    real8,dimension(:),allocatable,target::dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
    real8,dimension(:),pointer::pointer_dper_x,pointer_lper_x,pointer_mper_x
    real8,dimension(:),pointer::pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8,dimension(:,:),allocatable,target:: mat_nat_x,mat_per_y
    real8,dimension(:,:),pointer::pointer_mat_nat_x,pointer_mat_per_y
    int4 :: k,i,h,j     
    real8, dimension(:,:),pointer :: buf_y

    
    Nx = num_cells_x
    Ny = num_cells_y
   
    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mat_nat_x(0:Nx,0:Nx),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr)    
call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(buf_y(0:Ny,0:Ny))
    pointer_dper_y => dper_y
    pointer_lper_y => lper_y
    pointer_mper_y => mper_y    
    pointer_mat_per_y => mat_per_y

    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
 
    pointer_mat_nat_x=> mat_nat_x
    
    call s_compute_splines_coefs_matrix_nat_1d_new(pointer_mat_nat_x,Nx+1)


 do h=0, Ny-1 
    do j=0,Ny-1
       do i=0,Nx
          do k=0,Nx
             mat_spl2D(h*(Nx+1)+i+1,j*(Nx+1)+k+1)=mat_per_y(h,j)*mat_nat_x(i,k)
          end do
       end do
    enddo
 end do

    deallocate(mat_nat_x)
    deallocate(mat_per_y)
    deallocate(dper_y)
    deallocate(mper_y)
    deallocate(lper_y)
!    deallocate(buf_y)    

  end subroutine compute_D_spl2D_nat_per_noblock

  end module spline_module


  
