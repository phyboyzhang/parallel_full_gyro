module boris_rk_orbit
#include "work_precision.h"
  use gyroaverage_2d, only: gyroaverage_2d_plan
  use spline_module, only: &
       compute_spl2d_firstorder_derivative_point_per_per, &
       compute_spl2d_firstorder_derivative_point_polar_nat_per, &     
       compute_spl2d_firstorder_derivative_point_nat_per, & 
       compute_spl2d_field_point_per_per, &
       compute_spl2d_field_point_nat_per, &
       compute_spl2d_field_point_per_per
  
  use cartesian_mesh, only: cartesian_mesh_1d

  use field_2d_mesh, only: field_2d_plan
  use utilities_module, only: &
       polar_to_cartesian, &
       cartesian_to_polar
  
  implicit none

  public :: borissolve_2d_per_per, &
            borissolve_2d_nat_per, &
            fulrksolve_per_per, &
            fulrksolve_nat_per, &
            gyrorksolve,  &
            gyrorksolve_2ndorder
       
contains
 

subroutine borissolve_2d_per_per(x,v,&
       ep_weight, &
       Bf_weight, &
       eta_min, &
       eta_max, &
       delta_eta, &
       NC, &
       dtreal, &
       nt,geometry)
       
    real8, dimension(:,:), intent(in),pointer :: ep_weight, Bf_weight
!    real8, dimension(:,:), intent(in) :: NC  !!!==> number cells
    real8, dimension(:), intent(in) :: delta_eta, eta_min, eta_max
    real8, dimension(:), intent(inout) :: x,v
    real8, intent(in) :: dtreal
    int4,  intent(in) :: NC(2)
    int4,  intent(in) :: nt   ! the iteration number
    character(len=*),intent(in) :: geometry
    real8 :: magf(3)
    real8 :: pomat(3,3),nomat(3,3),xco(3),vel(3),vec(3,1)
    real8 :: deri_firstorder(2)
    real8 :: elef(3), x1(2)
!    int4, intent(in) :: flag
    
    int4 :: N=3, NRHS=1, LDA=3, LDB=3,INFO
    int4 :: IPIV(3)
    int4 i,j

          magf(3)= compute_spl2d_field_point_per_per(x(1:2), &
               Bf_weight, &
               eta_min,   &
               eta_max,   &
               delta_eta, &
               NC )
          magf(2)=0.0_f64
          magf(1)=0.0_f64

!   print*, "per_per,magf(3)=",magf(3)
   x1(1:2)=x(1:2)
   call obtain_interpolation_elefield_per_per(x1,    &
                                      elef,  &
                                      ep_weight,  &
                                      eta_min,    &
                                      eta_max,    &
                                      delta_eta,  &
                                      NC)
   
!   do i=1,NC(1)
!     print*, "ep_weight,i", ep_weight(i,:)
!   END DO 
   pomat(1,1)=0.0_f64
   pomat(1,2)=-magf(3)*dtreal/2.0_f64
   pomat(1,3)=magf(2)*dtreal/2.0_f64
   pomat(2,1)=magf(3)*dtreal/2.0_f64
   pomat(2,2)=0.0_f64
   pomat(2,3)=-magf(1)*dtreal/2.0_f64
   pomat(3,1)=-magf(2)*dtreal/2.0_f64
   pomat(3,2)=magf(1)*dtreal/2.0_f64
   pomat(3,3)=0.0_f64
   
   nomat(:,:)=-pomat(:,:)
   nomat(1,1)=nomat(1,1)+1.0_f64
   nomat(2,2)=nomat(2,2)+1.0_f64
   nomat(3,3)=nomat(3,3)+1.0_f64

   pomat(1,1)=pomat(1,1)+1.0_f64
   pomat(2,2)=pomat(2,2)+1.0_f64
   pomat(3,3)=pomat(3,3)+1.0_f64

   do i=1,3,1
    xco(i)=x(i)
    vel(i)=v(i)
   end do

   do i=1,3,1
     vec(i,1)=0.0_f64
     do j=1,3,1
     vec(i,1)=vec(i,1)+nomat(i,j)*vel(j)
     end do
     vec(i,1)=vec(i,1)+dtreal*elef(i)
   end do

   call dgesv(N,NRHS,pomat,LDA,IPIV,vec,LDB,INFO)

   do i=1,3,1
     xco(i)=xco(i)+vec(i,1)*dtreal
   end do

   do i=1,3,1
    x(i)=xco(i)
    v(i)=vec(i,1)
   end do
   return
 end subroutine borissolve_2d_per_per


subroutine borissolve_2d_nat_per(x,v,&
       ep_weight, &
       Bf_weight, &
       eta_min, &
       eta_max, &
       delta_eta, &
       NC, &
       dtreal, &
       nt,geometry)
       
    real8, dimension(:,:), intent(in),pointer :: ep_weight, Bf_weight
!    real8, dimension(:,:), intent(in) :: NC  !!!==> number cells
    real8, dimension(2), intent(in) :: delta_eta, eta_min, eta_max
    real8, dimension(3), intent(inout) :: x,v
    real8, intent(in) :: dtreal
    int4,  intent(in) :: NC(2)
    int4,  intent(in) :: nt   ! the iteration number
    character(len=*),intent(in) :: geometry
    real8 :: magf(3)
    real8 :: pomat(3,3),nomat(3,3),xco(3),vel(3),vec(3,1)
    real8 :: deri_firstorder(2)
    real8 :: elef(3), polar(3)
  !  int4, intent(in) :: flag
    
    int4 :: N=3, NRHS=1, LDA=3, LDB=3,INFO
    int4 :: IPIV(3)
    int4 i,j
    real8 :: x1(2)
    if(geometry=="polar") then
      call cartesian_to_polar(x(1:2),polar(1:2))
    else
      polar(1:2)=x(1:2)
    end if
          magf(3)= compute_spl2d_field_point_nat_per(polar(1:2),&
               Bf_weight, &
               eta_min,   &
               eta_max,   &
               delta_eta, &
               NC )
          magf(2)=0.0_f64
          magf(1)=0.0_f64

   x1(1:2)=x(1:2)
   call obtain_interpolation_elefield_nat_per(polar(1:2),    & 
                                      elef,  &
                                      ep_weight,  &
                                      eta_min,    &
                                      eta_max,    &
                                      delta_eta,  &
                                      NC,geometry)
 
!   do i=1,NC(1)
!     print*, "ep_weight,i", ep_weight(i,:)
!   END DO  
!  if(geometry=="nat_per") then
!   write(11,"(3f12.8)") elef
!  end if
 !  print*, "nat_per,elef=",elef                                
   pomat(1,1)=0.0_f64
   pomat(1,2)=-magf(3)*dtreal/2.0_f64
   pomat(1,3)=magf(2)*dtreal/2.0_f64
   pomat(2,1)=magf(3)*dtreal/2.0_f64
   pomat(2,2)=0.0_f64
   pomat(2,3)=-magf(1)*dtreal/2.0_f64
   pomat(3,1)=-magf(2)*dtreal/2.0_f64
   pomat(3,2)=magf(1)*dtreal/2.0_f64
   pomat(3,3)=0.0_f64
   
   nomat(:,:)=-pomat(:,:)
   nomat(1,1)=nomat(1,1)+1.0_f64
   nomat(2,2)=nomat(2,2)+1.0_f64
   nomat(3,3)=nomat(3,3)+1.0_f64

   pomat(1,1)=pomat(1,1)+1.0_f64
   pomat(2,2)=pomat(2,2)+1.0_f64
   pomat(3,3)=pomat(3,3)+1.0_f64

   do i=1,3,1
    xco(i)=x(i)
    vel(i)=v(i)
   end do

   do i=1,3,1
     vec(i,1)=0.0_f64
     do j=1,3,1
     vec(i,1)=vec(i,1)+nomat(i,j)*vel(j)
     end do
     vec(i,1)=vec(i,1)+dtreal*elef(i)
   end do

   call dgesv(N,NRHS,pomat,LDA,IPIV,vec,LDB,INFO)

   do i=1,3,1
     xco(i)=xco(i)+vec(i,1)*dtreal
   end do

   do i=1,3,1
    x(i)=xco(i)
    v(i)=vec(i,1)
   end do
   return
 end subroutine borissolve_2d_nat_per
 

subroutine fulrkfunc_f_per_per(f, &
       Bf_weight, &
       vec, &
       elef, &
       eta_min, &
       eta_max, &
       delta_eta, &
       NC,   &
       time)  ! It includes the magnetic field perturbation and the electric field perturbation
   real8, dimension(:,:), pointer, intent(in) :: Bf_weight
   real8, dimension(:), intent(in) :: elef, eta_min, eta_max,delta_eta
   real8, dimension(:), intent(in) :: vec
   real8, dimension(:), intent(inout) :: f
   real8, intent(in) :: time
   int4, dimension(:),intent(in) :: NC
!   character(len=*), intent(in) :: geometry
!   int4,  intent(in) :: flag
   
!   real(8) :: x1(3),v1(3),vperp,rholength,gc(3),time,k
   int4 :: i
   real8 x(3),v(3),mainb(3),minb(3),magf(3),k
   real8 :: polar(2)
 
   x(1)=vec(1)
   x(2)=vec(2)
   x(3)=vec(3)
   v(1)=vec(4)
   v(2)=vec(5)
   v(3)=vec(6)

         magf(3)= compute_spl2d_field_point_per_per(x(1:2),  &
               Bf_weight, &
               eta_min,   &
               eta_max,   &
               delta_eta, &
               NC )
          magf(2)=0.0_f64
          magf(1)=0.0_f64

   
   f(1)=vec(4)
   f(2)=vec(5)
   f(3)=vec(6)
   f(4)=(elef(1)+v(2)*magf(3)-v(3)*magf(2))
   f(5)=(elef(2)+v(3)*magf(1)-v(1)*magf(3))
   f(6)=(elef(3)+v(1)*magf(2)-v(2)*magf(1))

   return
 end subroutine fulrkfunc_f_per_per

subroutine fulrkfunc_f_nat_per(f,  &
       Bf_weight, &
       vec, &
       elef, &
       eta_min, &
       eta_max, &
       delta_eta, &
       NC,  &
       time,geometry)  ! It includes the magnetic field perturbation and the electric field perturbation
   real8, dimension(:,:), pointer, intent(in) :: Bf_weight
   real8, dimension(:),intent(in) :: eta_min, eta_max, delta_eta
   real8,intent(in) :: elef(3)
   real8, dimension(:), intent(in) :: vec
   real8, dimension(:), intent(inout) :: f
   int4, dimension(:), intent(in) :: NC
   real8, intent(in) :: time
   character(len=*),intent(in) :: geometry
   real8 :: polar(2)
   
   real8 x(3),v(3),mainb(3),minb(3),magf(3),k
 
   x(1)=vec(1)
   x(2)=vec(2)
   x(3)=vec(3)
   v(1)=vec(4)
   v(2)=vec(5)
   v(3)=vec(6)
       
       if(geometry=="polar") then
          call cartesian_to_polar(x(1:2),polar(1:2))
       else
         polar(1:2)=x(1:2)
       end if
          magf(3)= compute_spl2d_field_point_nat_per(polar(1:2),&
               Bf_weight, &
               eta_min,   &
               eta_max,   &
               delta_eta, &
               NC )
          magf(2)=0.0_f64
          magf(1)=0.0_f64
   
   f(1)=vec(4)
   f(2)=vec(5)
   f(3)=vec(6)
   f(4)=(elef(1)+v(2)*magf(3)-v(3)*magf(2))
   f(5)=(elef(2)+v(3)*magf(1)-v(1)*magf(3))
   f(6)=(elef(3)+v(1)*magf(2)-v(2)*magf(1))

   return
 end subroutine fulrkfunc_f_nat_per


!subroutine fulrkfunc_f_nat_per(f,  &
!       Bf_weight, &
!       vec, &
!       elef, &
!       eta_min, &
!       eta_max, &
!       delta_eta, &
!       NC,  &
!       time, geometry)  ! It includes the magnetic field perturbation and the electric field perturbation
!   real8, dimension(:,:), pointer, intent(in) :: Bf_weight
!   real8, dimension(:),intent(in) :: eta_min, eta_max, delta_eta
!   real8,intent(in) :: elef(3)
!   real8, dimension(:), intent(in) :: vec
!   real8, dimension(:), intent(inout) :: f
!   int4, dimension(:), intent(in) :: NC
!   real8, intent(in) :: time
!   character(len=*),intent(in) :: geometry
!   real8 :: polar(2)
!   
!!   real(8) :: x1(3),v1(3),vperp,rholength,gc(3),time,k
!   int4 :: i
!   real8 x(3),v(3),mainb(3),minb(3),magf(3),k
!!   real8 :: polar(2)
! 
!   x(1)=vec(1)
!   x(2)=vec(2)
!   x(3)=vec(3)
!   v(1)=vec(4)
!   v(2)=vec(5)
!   v(3)=vec(6)
!      
!        if(geometry="polar")
!          call cartesian_to_polar(x(1:2),polar(1:2))
!        else
!          polar(1:2)=x(1:2)
!        end if  
!          magf(3)= compute_spl2d_field_point_nat_per(polar(1:2),&
!               Bf_weight, &
!               eta_min,   &
!               eta_max,   &
!               delta_eta, &
!               NC )
!          magf(2)=0.0_f64
!          magf(1)=0.0_f64
!   
!   f(1)=vec(4)
!   f(2)=vec(5)
!   f(3)=vec(6)
!   f(4)=(elef(1)+v(2)*magf(3)-v(3)*magf(2))
!   f(5)=(elef(2)+v(3)*magf(1)-v(1)*magf(3))
!   f(6)=(elef(3)+v(1)*magf(2)-v(2)*magf(1))
!
!   return
! end subroutine fulrkfunc_f_polar_nat_per


 subroutine fulrksolve_per_per(vec, &
       field_2d, &
       m_x1,  &
       m_x2,  &
       dtreal,&
       time)  ! It includes the magnetic field perturbation and the electric field perturbation
   class(field_2d_plan),intent(in), pointer :: field_2d
   class(cartesian_mesh_1d), intent(in),pointer :: m_x1,m_x2
   real8, dimension(:), intent(inout) :: vec
   real8, intent(in) :: dtreal
   real8, intent(in) :: time
!   character(len=*), intent(in) :: geometry
!   int4,  intent(in) :: flag
   real8 :: f1(6), f2(6),f3(6),f4(6),x1(2)
   real8 :: elef(3), eta_min(2),eta_max(2),deri_firstorder(2),delta_eta(2)
   int4  :: NC(2)
   int4 :: i

   real8 :: vec1(6)

    vec1(:)=vec(:)        !for magnetic and electric perturbation
    eta_min(1)=m_x1%eta_min
    eta_min(2)=m_x2%eta_min
    eta_max(1)=m_x1%eta_max
    eta_max(2)=m_x2%eta_max
    NC(1)=m_x1%num_cells
    NC(2)=m_x2%num_cells
    delta_eta(1)=m_x1%delta_eta
    delta_eta(2)=m_x2%delta_eta

  x1(1:2)=vec(1:2)    
  call obtain_interpolation_elefield_per_per(x1,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC)                                              
 
  write(12, "(3f19.9)") elef  
  call fulrkfunc_f_per_per(f1,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta,NC,time)
  
  vec1(:)=vec(:)+0.5_f64*dtreal*f1(:)
  x1(1:2)=vec1(1:2)    
  call obtain_interpolation_elefield_per_per(x1,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC)                                                  
  call fulrkfunc_f_per_per(f2,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta,NC,time+dtreal*0.5_f64)
  
  vec1(:)=vec(:)+0.5_f64*dtreal*f2(:)
  x1(1:2)=vec1(1:2)    
  call obtain_interpolation_elefield_per_per(x1,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC)   
  call fulrkfunc_f_per_per(f3,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta,NC,time+dtreal*0.5_f64)

  
  vec1(:)=vec(:)+dtreal*f3(:)
  x1(1:2)=vec1(1:2)    
  call obtain_interpolation_elefield_per_per(x1,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC)   
  call fulrkfunc_f_per_per(f4,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta,NC,time+dtreal)

   vec(:)=vec(:)+(dtreal/6.0_f64)*(f1(:)+2.0_f64*f2(:)+2.0_f64*f3(:)+f4(:))


 end subroutine fulrksolve_per_per
 

 subroutine fulrksolve_nat_per(vec, &
       field_2d, &
       m_x1,  &
       m_x2,  &
       dtreal,&
       time,geometry)  ! It includes the magnetic field perturbation and the electric field perturbation
   class(field_2d_plan),pointer,intent(in) :: field_2d
   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
   real8, dimension(:), intent(inout) :: vec
   real8, intent(in) :: dtreal
   real8, intent(in) :: time
   character(len=*), intent(in) :: geometry
 !  int4,  intent(in) :: flag
   real8 :: f1(6), f2(6),f3(6),f4(6),x1(2),polar(2)
   real8 :: elef(3), eta_min(2),eta_max(2),deri_firstorder(2),delta_eta(2)
   int4  :: NC(2)
   int4 :: i

   real8 :: vec1(6)
    vec1(:)=vec(:)        !for magnetic and electric perturbation
    eta_min(1)=m_x1%eta_min
    eta_min(2)=m_x2%eta_min
    eta_max(1)=m_x1%eta_max
    eta_max(2)=m_x2%eta_max
    NC(1)=m_x1%num_cells
    NC(2)=m_x2%num_cells
    delta_eta(1)=m_x1%delta_eta
    delta_eta(2)=m_x2%delta_eta
   
  if(geometry=="polar") then 
    call cartesian_to_polar(vec(1:2),polar(1:2))
  else 
    polar(1:2)=vec(1:2)
  end if
  call obtain_interpolation_elefield_nat_per(polar,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC,geometry)                                              
   call fulrkfunc_f_nat_per(f1,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta,NC,time,geometry) 
  vec1(:)=vec(:)+0.5_f64*dtreal*f1(:)
  if(geometry=="polar") then 
    call cartesian_to_polar(vec1(1:2),polar(1:2))
  else 
    polar(1:2)=vec1(1:2)
  end if
  call obtain_interpolation_elefield_nat_per(polar,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC,geometry)                                                  
  call fulrkfunc_f_nat_per(f2,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta, &
                                 NC,time+dtreal*0.5_f64,geometry)
  
  vec1(:)=vec(:)+0.5_f64*dtreal*f2(:)
  if(geometry=="polar") then 
    call cartesian_to_polar(vec1(1:2),polar(1:2))
  else 
    polar(1:2)=vec1(1:2)
  end if
call obtain_interpolation_elefield_nat_per(polar,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC,geometry)   
  call fulrkfunc_f_nat_per(f3,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta, &
                                 NC,time+dtreal*0.5_f64,geometry)

  
  vec1(:)=vec(:)+dtreal*f3(:)
  if(geometry=="polar") then 
    call cartesian_to_polar(vec1(1:2),polar(1:2))
  else 
    polar(1:2)=vec(1:2)
  end if
call obtain_interpolation_elefield_nat_per(polar,    &
                                     elef,  &
                                     field_2d%ep_weight_init,  &
                                     eta_min,    &
                                     eta_max,    &
                                     delta_eta,  &
                                     NC,geometry)   
  call fulrkfunc_f_nat_per(f4,field_2d%Bf_weight_3rd,vec1,elef,eta_min,eta_max,delta_eta, &
                                 NC,time+dtreal,geometry)

   vec(:)=vec(:)+(dtreal/6.0_f64)*(f1(:)+2.0_f64*f2(:)+2.0_f64*f3(:)+f4(:))


 end subroutine fulrksolve_nat_per
 
 subroutine gyrorksolve_f(x,&
      mu, &
      Bf_weight, &
      gyro_weight, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,     &
      f,      &
      dtreal, &
      geometry,boundary)
   real8, intent(in) :: mu
   real8, dimension(:),intent(inout) :: x
   real8,dimension(:,:),pointer,intent(in) :: gyro_weight, Bf_weight
   real8,dimension(:),intent(in) :: eta_min,eta_max, delta_eta
   int4, intent(in) :: NC(2)
   real8,dimension(:),intent(inout) :: f
   real8,intent(in) :: dtreal
   character(len=*),intent(in) :: geometry,boundary
 !  int4, intent(in) :: flag
   real8 :: deri_ep(2), deri_Bf(2), magf(3)
   real8 :: x_polar(2)


   select case (geometry)
   case ("cartesian")
   select case (boundary)
   case ("double_per")
   call  compute_spl2d_firstorder_derivative_point_per_per(x, &
        gyro_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_ep)

   call  compute_spl2d_firstorder_derivative_point_per_per(x, &
        Bf_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_Bf)

      magf(3)=compute_spl2d_field_point_per_per(x, &
        Bf_weight,  &
        eta_min,    &
        eta_max,    &
        delta_eta,  &
        NC)

   f(1)=(-deri_ep(2)-mu*deri_Bf(2))/magf(3)
   f(2)=(deri_ep(1)+mu*deri_Bf(1))/magf(3)     
!   f(1)=-1.0_f64*deri_ep(2)
!   f(2)=deri_ep(1)
!print*, "gyro,per_per,deri=", deri_ep
    case ("nat_per")
    call  compute_spl2d_firstorder_derivative_point_nat_per(x_polar, &
        gyro_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_ep)

    call  compute_spl2d_firstorder_derivative_point_nat_per(x_polar, &
        Bf_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_Bf)

    magf(3)=compute_spl2d_field_point_nat_per(x_polar, &
        Bf_weight,  &
        eta_min,    &
        eta_max,    &
        delta_eta,  &
        NC)

   f(1)=(-deri_ep(2)-mu*deri_Bf(2))/magf(3)
   f(2)=(deri_ep(1)+mu*deri_Bf(1))/magf(3)
      case default    
        print*, "boundary is not right, boundary=", boundary
      end select

   case ("polar")


   case default
      print*, "geometry is not right, geometry=", geometry
   end select

 end subroutine gyrorksolve_f


subroutine gyrorksolve_f_2ndorder(x,&
      mu, &
      Bf_weight, &
      gyro_weight, &
      driftsquare_weight, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,     &
      f,      &
      dtreal, &
      geometry, boundary)
   real8, intent(in) :: mu
   real8, dimension(:),intent(inout) :: x
   real8,dimension(:,:),intent(in),pointer :: gyro_weight,Bf_weight,driftsquare_weight
   real8,dimension(:),intent(in) :: eta_min,eta_max, delta_eta
   int4, intent(in) :: NC(2)
   real8,dimension(:),intent(inout) :: f
   real8,intent(in) :: dtreal
   character(len=*),intent(in) :: geometry,boundary
 !  int4, intent(in) :: flag
   real8 :: deri_ep(2), deri_Bf(2), magf(3), deri_driftsquare(2)
   real8 :: x_polar(2)


   select case (geometry)
   case ("cartesian")
   select case (boundary)
   case ("double_per")
   call  compute_spl2d_firstorder_derivative_point_per_per(x, &
        gyro_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_ep)

 call  compute_spl2d_firstorder_derivative_point_per_per(x, &
        Bf_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_Bf)

   call compute_spl2d_firstorder_derivative_point_per_per(x, &
        driftsquare_weight, &
        eta_min,  &
        eta_max,  &
        delta_eta, &
        NC,       &
        deri_driftsquare)

    magf(3)=compute_spl2d_field_point_per_per(x_polar, &
        Bf_weight,  &
        eta_min,    &
        eta_max,    &
        delta_eta,  &
        NC)

!   f(1)=-1.0_f64*deri_ep(2)-0.5_f64*deri_driftsquare(2)
!   f(2)=deri_ep(1)+0.5_f64*deri_driftsquare(1)
  f(1)=(-deri_ep(2)-mu*deri_Bf(2))/magf(3)-0.5_f64*deri_driftsquare(2)/magf(3)
   f(2)=(deri_ep(1)+mu*deri_Bf(1))/magf(3)+0.5_f64*deri_driftsquare(1)/magf(3)
!print*, "deri_driftsquare=", deri_driftsquare
   case ("nat_per")
    call  compute_spl2d_firstorder_derivative_point_nat_per(x, &
        gyro_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_ep)

    call  compute_spl2d_firstorder_derivative_point_nat_per(x, &
        Bf_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_Bf)

   call compute_spl2d_firstorder_derivative_point_nat_per(x, &
        driftsquare_weight, &
        eta_min,  &
        eta_max,  &
        delta_eta, &
        NC,       &
        deri_driftsquare)

  magf(3)=compute_spl2d_field_point_nat_per(x, &
        Bf_weight,  &
        eta_min,    &
        eta_max,    &
        delta_eta,  &
        NC)



   f(1)=(-deri_ep(2)-mu*deri_Bf(2))/magf(3)-0.5_f64*deri_driftsquare(2)/magf(3)
   f(2)=(deri_ep(1)+mu*deri_Bf(1))/magf(3)+0.5_f64*deri_driftsquare(1)/magf(3)


   case default
     stop
   end select

   case ("polar")
    call cartesian_to_polar(x,x_polar)
    call  compute_spl2d_firstorder_derivative_point_polar_nat_per(x_polar, &
        gyro_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_ep)

    call  compute_spl2d_firstorder_derivative_point_polar_nat_per(x_polar, &
        Bf_weight,  &
        eta_min,      &
        eta_max,      &
        delta_eta,    &
        NC,           &
        deri_Bf)

  call compute_spl2d_firstorder_derivative_point_polar_nat_per(x_polar, &
        driftsquare_weight, &
        eta_min,  &
        eta_max,  &
        delta_eta, &
        NC,       &
        deri_driftsquare)

   magf(3)=compute_spl2d_field_point_nat_per(x_polar, &
        Bf_weight,  &
        eta_min,    &
        eta_max,    &
        delta_eta,  &
        NC)



   f(1)=(-deri_ep(2)-mu*deri_Bf(2))/magf(3)-0.5_f64*deri_driftsquare(2)/magf(3)
   f(2)=(deri_ep(1)+mu*deri_Bf(1))/magf(3)+0.5_f64*deri_driftsquare(1)/magf(3)

   case default
     print*, "geometry should be chosen as polar when flag=0"
   end select

 end subroutine gyrorksolve_f_2ndorder


 subroutine gyrorksolve(x,&
   !         mu, &
            Bf_weight, &
            gyro, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry,boundary)
 !  real8, intent(in) :: mu
   real8, dimension(:,:), pointer,intent(in) :: Bf_weight
   class(cartesian_mesh_1d),intent(in),pointer :: m_x1,m_x2
   class(gyroaverage_2d_plan), intent(in),pointer :: gyro
   int4, intent(in) :: n  !!!  n=40000
   !   int4, intent(in) :: flag
   character(len=*), intent(in) :: geometry,boundary
   real8, dimension(:) :: x
   real8,intent(in) :: mu
   real8,intent(in) :: dtreal
   real8 :: eta_min(2), eta_max(2), delta_eta(2)
   int4 :: NC(2)
   int4 :: i

   real8 :: x1(2),f1(2),f2(2),f3(2),f4(2)

   real8 :: deri_firstorder(2)

   eta_min(1)=m_x1%eta_min
   eta_min(2)=m_x2%eta_min
   eta_max(1)=m_x1%eta_max
   eta_max(2)=m_x2%eta_max
   delta_eta(1)= m_x1%delta_eta
   delta_eta(2)= m_x2%delta_eta
   NC(1)=m_x1%num_cells
   NC(2)=m_x2%num_cells
 
     x1=x
     call gyrorksolve_f(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f1,&
      dtreal, &
      geometry,boundary)
     x1=x+0.5_f64*dtreal*f1
 !    print*, "f1=", f1
     call gyrorksolve_f(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f2,&
      dtreal, &
      geometry,boundary)
     x1=x+0.5_f64*dtreal*f2
  !    print*, "x1=", x1, "f2=",f2
  
     call gyrorksolve_f(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f3,&
      dtreal, &
      geometry,boundary)
     x1=x+dtreal*f3
    !  print*, "x1=", x1, "f3=",f3
     
       call gyrorksolve_f(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f4,&
      dtreal, &
      geometry,boundary)
   !   print*, "x1=", x1, "f4=",f4

      x(:)=x(:)+(dtreal/6.0_f64)*(f1(:)+2.0_f64*f2(:)+2.0_f64*f3(:)+f4(:))
   
    end subroutine gyrorksolve
    
    
    subroutine gyrorksolve_2ndorder(x,&
   !         mu, &
            Bf_weight, &
            gyro, &
            driftsquare_weight, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry, &
            boundary)
 !  real8, intent(in) :: mu
   real8, dimension(:,:), pointer,intent(in) :: Bf_weight, driftsquare_weight
   class(cartesian_mesh_1d),intent(in),pointer :: m_x1,m_x2
   class(gyroaverage_2d_plan), intent(in),pointer :: gyro
   int4, intent(in) :: n  !!!  n=40000
   !   int4, intent(in) :: flag
   character(len=*), intent(in) :: geometry, boundary
   real8, dimension(:) :: x
   real8,intent(in) :: mu
   real8,intent(in) :: dtreal
   real8 :: eta_min(2), eta_max(2), delta_eta(2)
   int4 :: NC(2)
   int4 :: i

   real8 :: x1(2),f1(2),f2(2),f3(2),f4(2)

   real8 :: deri_firstorder(2)

   eta_min(1)=m_x1%eta_min
   eta_min(2)=m_x2%eta_min
   eta_max(1)=m_x1%eta_max
   eta_max(2)=m_x2%eta_max
   delta_eta(1)= m_x1%delta_eta
   delta_eta(2)= m_x2%delta_eta
   NC(1)=m_x1%num_cells
   NC(2)=m_x2%num_cells
   x1(1)=2.1
   x1(2)=2.1

     x1=x
     call gyrorksolve_f_2ndorder(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      driftsquare_weight, &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f1,&
      dtreal, &
      geometry, &
      boundary)
     x1=x+0.5_f64*dtreal*f1

    call gyrorksolve_f_2ndorder(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      driftsquare_weight,  &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f2,&
      dtreal, &
      geometry, boundary)
     x1=x+0.5_f64*dtreal*f2

     call gyrorksolve_f_2ndorder(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      driftsquare_weight,  &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f3,&
      dtreal, &
      geometry, boundary)
     x1=x+dtreal*f3

    call gyrorksolve_f_2ndorder(x1,&
      mu, &
      Bf_weight, &
      gyro%ep_weight_gyro, &
      driftsquare_weight,  &
      eta_min,&
      eta_max,&
      delta_eta, &
      NC,&
      f4,&
      dtreal, &
      geometry,boundary)
      x(:)=x(:)+(dtreal/6.0_f64)*(f1(:)+2.0_f64*f2(:)+2.0_f64*f3(:)+f4(:))

    end subroutine gyrorksolve_2ndorder


    subroutine obtain_interpolation_elefield_per_per(x1,    &
                                             elef,  &
                                             ep_weight,  &
                                             eta_min,    &
                                             eta_max,    &
                                             delta_eta,  &
                                             NC)
                                             
      real8, dimension(2),intent(in) :: x1,eta_min, eta_max, delta_eta
      real8, dimension(:,:),pointer,intent(in) :: ep_weight
      real8, dimension(3), intent(inout) :: elef
      int4,  dimension(2), intent(in) :: NC
      real8 :: deri_firstorder(2)
      
       call compute_spl2d_firstorder_derivative_point_per_per(x1, &
            ep_weight, &
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            deri_firstorder)
       elef(1)=-deri_firstorder(1)
       elef(2)=-deri_firstorder(2)
       elef(3)= 0.0_f64
     
     end subroutine obtain_interpolation_elefield_per_per
     
     subroutine obtain_interpolation_elefield_nat_per(x_polar,    &
                                             elef,  &
                                             ep_weight,  &
                                             eta_min,    &
                                             eta_max,    &
                                             delta_eta,  &
                                             NC,geometry)
                                             
      real8, dimension(2),intent(in) :: x_polar,eta_min, eta_max, delta_eta
      real8, dimension(:,:),intent(in),pointer :: ep_weight
      real8, dimension(3), intent(inout) :: elef
      int4,  dimension(2), intent(in) :: NC
      character(len=*),intent(in) :: geometry
      real8 :: deri_cart(2)
      if(geometry=="polar") then
       call compute_spl2d_firstorder_derivative_point_polar_nat_per(x_polar, &
            ep_weight, &
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            deri_cart)
      else
       call compute_spl2d_firstorder_derivative_point_nat_per(x_polar, &
            ep_weight, &
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            deri_cart)
      end if
       elef(1)=-deri_cart(1)
       elef(2)=-deri_cart(2)
       elef(3)= 0.0_f64

     end subroutine obtain_interpolation_elefield_nat_per    
    

end module boris_rk_orbit

