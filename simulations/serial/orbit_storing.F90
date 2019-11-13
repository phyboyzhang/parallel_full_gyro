module orbit_storing
#include "work_precision.h"
  use boris_rk_orbit, only: &
            borissolve_2d_per_per, &
            borissolve_2d_nat_per, &
            fulrksolve_per_per, &
            fulrksolve_nat_per, &
            gyrorksolve,  &
            gyrorksolve_2ndorder  

  
!  use spline_module, only: compute_spl2d_firstorder_derivatve_point
  use cartesian_mesh, only: cartesian_mesh_1d
  use field_2d_mesh, only: field_2d_plan
  use gyroaverage_2d, only: gyroaverage_2d_plan
  use field_initialize, only: &
       efield2d_anal, &
       magfield_2d_interpolation_point
       

  implicit none
  public :: borisrun, &
           fulrkorbit, &
           gyrorkorbit,&
           borisrun_anal 

contains 
 subroutine borisrun(x,v,field_2d,m_x1,m_x2,n,dtreal,geometry,boundary)
   class(cartesian_mesh_1d),pointer, intent(in) :: m_x1,m_x2
   class(field_2d_plan), pointer,intent(in) :: field_2d
   int4, intent(in) :: n  !!!  n=40000
   character(len=*), intent(in) :: geometry,boundary
!   int4, intent(in) :: flag
   real8, dimension(:), intent(inout) :: x,v
   real8, intent(in) :: dtreal
   real8 :: eta_min(2), eta_max(2), delta_eta(2)
   int4 :: NC(2)
   
   real(8) :: x1(3),v1(3),time,k
   int4 :: i
   character(99) :: char_a,char_b,char_c,char_d,char_f

   real8 :: vperp, rholength, gc(3)

   write(unit=char_a,fmt=*) n
   char_b="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/borisorbit_be_"
   char_f="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/borisorbit_gc_"
 !   char_b="/Users/zhang/program/debug/electrostatic_exp/run/borisorbit_be_"
 !  char_f="/Users/zhang/program/debug/electrostatic_exp/run/borisorbit_gc_"
   
   x1=x
   v1=v
   eta_min(1)=m_x1%eta_min
   eta_min(2)=m_x2%eta_min
   eta_max(1)=m_x1%eta_max
   eta_max(2)=m_x2%eta_max
   delta_eta(1)=m_x1%delta_eta
   delta_eta(2)=m_x2%delta_eta
   NC(1)=m_x1%num_cells
   NC(2)=m_x2%num_cells

!   print*, "boris init", x, v
   !!!! open(30,file='/Users/zhang/program/electrostatic_exp/run/boris_orbit_nbe.txt',status='replace')
   open(40,file=trim(char_b)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
   open(20,file=trim(char_f)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
   
 do i=1,n

!!$   call borissolve_be_no_fieldinterpolation_2d(x,v,&
!!$       field_2d%ep_weight_init, &
!!$       eta_min, &
!!$       eta_max, &
!!$       delta_eta, &
!!$       NC, &
!!$       dtreal, &
!!$       i, &
!!$       flag)
    select case (geometry)
    case ("cartesian")
      select case(boundary)
       case("double_per")
       call  borissolve_2d_per_per(x,v,&
            field_2d%ep_weight_init, &
            field_2d%bf_weight_3rd, &    
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            dtreal, &
            i,geometry)
       case ("nat_per")

       call  borissolve_2d_nat_per(x,v,&
            field_2d%ep_weight_init, &
            field_2d%bf_weight_3rd, &    
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            dtreal, &
            i,geometry)
       case default
            stop
       end select
     case("polar")
      call  borissolve_2d_nat_per(x,v,&
            field_2d%ep_weight_init, &
            field_2d%bf_weight_3rd, &    
            eta_min, &
            eta_max, &
            delta_eta, &
            NC, &
            dtreal, &
            i,geometry)
     case default
       stop
     end select  

      vperp=sqrt(v(1)**2+v(2)**2)
      rholength=vperp
      
      gc(1)=x(1)+rholength*v(2)/vperp
      gc(2)=x(2)-rholength*v(1)/vperp

!=====> Here, judge whether the particle reaches out of the boundary
      if(x(1).le.eta_min(1).or.x(1).ge.eta_max(1).or.x(2).le.eta_min(2)  &
           .or.x(2).ge.eta_max(2)) then
         print*,"i=",i, "boris,particle reaches out of the boundary"
         goto 200
      else
         write(20,'(I6,2F18.6)') i,gc(1),gc(2)
         write(40,'(I6,6f18.6)') i,x(1),x(2),x(3),v(1),v(2),v(3)
     end if
 end do

200 close(40)
    close(20)
    
 end subroutine borisrun


subroutine borisrun_anal(x,v, &
           wam, wave_one,wave_two,contr,dtreal,n,geometry)
  real8,dimension(:) :: x,v
  real8, intent(in) :: wam,wave_one,wave_two, contr,dtreal
  int4, intent(in) :: n ! n is the iterative number
  character(len=*), intent(in) :: geometry
  character(99) :: char_d, char_a
  int4 :: flag=0
  int4 :: i

   write(unit=char_a,fmt=*) n
   char_d="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/borisorbit_anal_"  
!   char_d="/Users/zhang/program/debug/electrostatic_exp/run/borisorbit_anal_"  
   open(30,file=trim(char_d)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
  do i=1, n
     call borissolve_be_2d_anal(x,v,dtreal, &
          wam,  &
          wave_one,      &
          wave_two,      &
          contr,   &
          flag)

     !      vperp=sqrt(v1(1)**2+v1(2)**2)
!      rholength=vperp/(2.0*pi*bamp(x1,k))
!      
!      gc(1)=x1(1)+rholength*v(2)/vperp
!      gc(2)=x1(2)-rholength*v(1)/vperp
!      gc(3)=x1(3)
!      write(50,'(I6,3F12.6)') i,gc(1),gc(2),gc(3)
     write(30,'(I6,6f12.6)') i,x(1),x(2),x(3),v(1),v(2),v(3)     
  end do
  close(30)
  
end subroutine borisrun_anal

 subroutine fulrkorbit(x,v, &
            field_2d, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry,boundary)
   class(cartesian_mesh_1d),pointer,intent(in) :: m_x1,m_x2
   class(field_2d_plan), pointer,intent(in) :: field_2d
   int4, intent(in) :: n  !!!  n=40000
   !  int4, intent(in) :: flag
   character(len=*),intent(in) :: geometry,boundary
   real8, dimension(:), intent(inout) :: x,v
   real8,intent(in) :: mu
   real8,intent(in) :: dtreal
   
   real8 :: eta_min(2), eta_max(2), delta_eta(2), elef(3)
   real8 :: deri_firstorder(2), x1(2)
   int4 :: NC(2)
   int4 :: i

   real8 vec(6),f1(6),f2(6),f3(6),f4(6),vec_1(6)
   character(99) :: char_a,char_b,char_c,char_d,char_f

   real8 :: vperp, rholength, gc(3), time 

   write(unit=char_a,fmt=*) n
   char_b="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/fulrkorbit_be_"
   char_f="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/fulrkorbit_gc_"
!    char_b="/Users/zhang/program/debug/electrostatic_exp/run/fulrkorbit_be_"
!   char_f="/Users/zhang/program/debug/electrostatic_exp/run/fulrkorbit_gc_"
    
    time=0.0
    vec(1)=x(1)
    vec(2)=x(2)
    vec(3)=x(3)
    vec(4)=v(1)
    vec(5)=v(2)
    vec(6)=v(3)
    vec_1(:)=vec(:)

   eta_min(1)=m_x1%eta_min
   eta_min(2)=m_x2%eta_min
   eta_max(1)=m_x1%eta_max
   eta_max(2)=m_x2%eta_max
   delta_eta(1)=m_x1%delta_eta
   delta_eta(2)=m_x2%delta_eta
   NC(1)=m_x1%num_cells
   NC(2)=m_x2%num_cells

      open(10,file=trim(char_b)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
      open(70,file=trim(char_f)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
!    open(40,file='/Users/zhang/program/electrostatic_exp/run/gcful_be.txt',status='replace')

    x1(1)=x(1)
    x1(2)=x(2)
!  print*, "vec", vec 
   do i=1,n
       select case (geometry)
        case("cartesian")
        select case(boundary)
          case("double_per")
             call fulrksolve_per_per(vec, &
             field_2d, &
             m_x1, &
             m_x2, &
             dtreal, &
             time)
!print*, "field_2d%ep_weight_init(10,:)=",field_2d%ep_weight_init(10,:)
!print*, "field_2d%ep_weight_init(:,10)=",field_2d%ep_weight_init(:,10)
          case ("nat_per")      
             call  fulrksolve_nat_per(vec, &
             field_2d, &
             m_x1, &
             m_x2, &
             dtreal, &
             time,geometry)
          case default
            print*, "boundary is not right, boundary=", boundary
            stop
          end select
          
       case("polar")
             call  fulrksolve_nat_per(vec, &
             field_2d, &
             m_x1, &
             m_x2, &
             dtreal, &
             time,geometry)
        case default
          stop
        end select

      if(vec(1).le.eta_min(1).or.vec(1).ge.eta_max(1).or.vec(2).le.eta_min(2)  &
           .or.vec(2).ge.eta_max(2)) then
         print*,"i=",i, "fulrk4, particle reaches out of the boundary"
         goto 150
      else
      vperp=sqrt(vec(4)**2+vec(5)**2)
      rholength=vperp
      
      gc(1)=vec(1)+rholength*vec(5)/vperp
      gc(2)=vec(2)-rholength*vec(4)/vperp
      
       write(70,'(I6,2F18.6)') i,gc(1),gc(2)    
   
       write(10,'(I6,6f18.6)') i,vec(1),vec(2),vec(3),vec(4),vec(5),vec(6) 
    end if 
  
  end do


150  close(10)
   close(70)
   

 end subroutine fulrkorbit


 subroutine fulrkorbit_anal(vec, &
           wam, wave_one,wave_two,contr,dtreal,n,geometry)
  real8, dimension(:) :: vec
  real8, intent(in) :: wam,wave_one,wave_two, contr,dtreal
  int4, intent(in) :: n
  character(len=*), intent(in) :: geometry
  character(99) :: char_d, char_a
  real8 :: elef(3), x(2)
  int4 :: flag=0
  real8 :: time,time_1
  int4 :: i


   write(unit=char_a,fmt=*) n 
   char_d="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/fulrkorbit_anal_"
! char_d="/Users/zhang/program/debug/electrostatic_exp/run/fulrkorbit_anal_"

   open(20,file=trim(char_d)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')
   time=0._f64
   time=0._f64
   do i=1,n
      x(1:2)=vec(1:2)
      time=real(i-1,8)*dtreal
      time_1=time
      call  efield2d_anal(elef,  &
           wam,  &
           wave_one, &
           wave_two, &
           x)

      call fulrksolve_noninterpolation( &
           vec,    &
           elef,   &
           dtreal, &
           time_1,   &
           flag )
    write(20,'(I6,6f18.6)') i,vec(1),vec(2),vec(3),vec(4),vec(5),vec(6) 
 end do

 close(20)
 
 end subroutine fulrkorbit_anal

 subroutine gyrorkorbit(x,&
            Bf_weight, &
            gyro, &
            driftsquare_weight, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry,boundary,orbitorder)
            
   real8,dimension(:,:),pointer,intent(in) :: Bf_weight,driftsquare_weight
   class(cartesian_mesh_1d),pointer,intent(in) :: m_x1,m_x2
   class(gyroaverage_2d_plan),pointer, intent(in) :: gyro
   int4, intent(in) :: n  !!!  n=40000
   character(len=*), intent(in) :: geometry,boundary
   int4, intent(in) :: orbitorder
!   int4, intent(in) :: flag
   real8, dimension(:) :: x
   real8,intent(in) :: mu
   real8,intent(in) :: dtreal
   real8 :: eta_min(2), eta_max(2), delta_eta(2)
   int4 :: NC(2)
   int4 :: i

   real8 :: x1(2),f1(2),f2(2),f3(2),f4(2)
   character(99) :: char_a,char_b,char_c

   real8 :: deri_firstorder(2)

   
   write(unit=char_a,fmt=*) n
   char_b="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/gyro_rkorbit_"
!  char_b="/Users/zhang/program/debug/electrostatic_exp/run/gyro_rkorbit_"

   eta_min(1)=m_x1%eta_min
   eta_min(2)=m_x2%eta_min
   eta_max(1)=m_x1%eta_max
   eta_max(2)=m_x2%eta_max
   delta_eta(1)= m_x1%delta_eta
   delta_eta(2)= m_x2%delta_eta
   NC(1)=m_x1%num_cells
   NC(2)=m_x2%num_cells
   x1(1)=2.1   !!!!! The way to assian the initial position is not good! 
   x1(2)=2.2

   open(60,file=trim(char_b)//trim(geometry)//trim(adjustl(char_a))//".txt",status='replace')   
   do i=1,N
      
      if(orbitorder==1) then
         call  gyrorksolve(x,&
            Bf_weight, &
            gyro, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry,boundary)
      else
        call  gyrorksolve_2ndorder(x,&
            Bf_weight, &
            gyro, &
            driftsquare_weight, &
            m_x1, &
            m_x2, &
            mu,   &
            n,    &
            dtreal,&
            geometry,boundary)
      end if 
      if(x(1).le.eta_min(1).or.x(1).ge.eta_max(1).or.x(2).le.eta_min(2)  &
           .or.x(2).ge.eta_max(2)) then
         print*,"i=",i, "gyrork4,particle reaches out of the boundary"
         goto 100
      else    
         write(60, '(I4,8F12.6)') i, x(1), x(2)
 !     print*, i, x(1),x(2)
      end if

    end do

 100  close(60)
   
 end subroutine gyrorkorbit



  subroutine borissolve_nbe_no_fieldinterpolation_2d(x,v, &
       magf, &
       eta_min, &
       eta_max, &
       eta_delta, &
       NC, &
       dtreal, &
       nt,     &
       flat)
       
    !   real8, dimension(:,:), intent(in) :: ep_weight
    real8, intent(in) :: magf(3)
 !   real8, dimension(:,:), intent(in) :: NC  !!!==> number cells
    real8, dimension(2), intent(in) :: eta_delta, eta_min, eta_max
    real8, dimension(3), intent(inout) :: x,v
    real8, intent(in) :: dtreal
    int4,  intent(in) :: NC(2)
    int4,  intent(in) :: nt   ! the iteration number
    int4,  intent(in) :: flat
    
    real8 :: minb(3)
    real8  pomat(3,3),nomat(3,3),xco(3),vel(3),vec(3,1),vec1(3,1)
    real8 :: deri_firstorder(2)


    int4 :: N=3, NRHS=1, LDA=3, IPIV=3, LDB=3,INFO
    int4 i,j
    
 !  call magfield_2d_interpolation_point(x,magf,flat)
!   call minbfield(x,minb,k)

 !  call compute_spl2d_firstorder_derivatve_point(x, val_first_deriv,val_second_deriv,eta_delta)   

   pomat(1,1)=0._f64
   pomat(1,2)=-magf(3)*dtreal/2.0_f64
   pomat(1,3)=magf(2)*dtreal/2.0_f64
   pomat(2,1)=magf(3)*dtreal/2.0_f64
   pomat(2,2)=0._f64
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

 !  rmat=pomat.x.nomat

   do i=1,3,1
     vec(i,1)=0.0_f64
     do j=1,3,1
     vec(i,1)=vec(i,1)+nomat(i,j)*vel(j)
     end do
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
 end subroutine borissolve_nbe_no_fieldinterpolation_2d


 subroutine borissolve_be_2d_anal(x,v,dtreal,wp,wn_one,wn_two,contr,flag)
  real8, dimension(:) :: x,v
  real8, intent(in) :: dtreal,wp,wn_one,wn_two,contr
  int4,  intent(in) :: flag
  real8 :: elef(3)
  real8 :: pomat(3,3),nomat(3,3),xco(3),vel(3),vec(3,1),magf(3)
  
    int4 :: N=3, NRHS=1, LDA=3, LDB=3,INFO
    int4 :: IPIV(3)
    int4 i,j

   if(flag==0) then
      call magfield_2d_interpolation_point(x,magf,flag)
   end if
   !  call minbfield(x,minb
   !  magf(:)=mainb(:)
   call efield2d_anal(elef,wp,wn_one,wn_two,x)
   
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
 end subroutine borissolve_be_2d_anal


subroutine fulrksolve_noninterpolation( &
       vec, &
       elef, &
       dtreal,&
       time, &
       flag)  ! It includes the magnetic field perturbation and the electric field perturbation   real8 :: elef(3)
  real8, dimension(:), intent(inout) :: vec
  real8, dimension(:), intent(in) :: elef
   real8, intent(in) :: dtreal
   real8, intent(in) :: time
!   character(len=*), intent(in) :: geometry
   int4,  intent(in) :: flag
   real8 :: f1(6), f2(6),f3(6),f4(6)
   int4 :: i

   real8 :: vec1(6)

    vec1(:)=vec(:)        !for magnetic and electric perturbation

 !   print*, "vec1",vec1
   call fulrkfunc_be_step(vec1,elef,f1,time,flag)
   vec1(:)=vec(:)+0.5_f64*dtreal*f1(:)
!   print*, "vec1", vec1(:), f1(:)
   call fulrkfunc_be_step(vec1,elef,f2,time+dtreal*0.5_f64,flag)
   vec1(:)=vec(:)+0.5_f64*dtreal*f2(:)

   call fulrkfunc_be_step(vec1,elef,f3,time+dtreal*0.5_f64,flag)
   vec1(:)=vec(:)+dtreal*f3(:)

   call fulrkfunc_be_step(vec1,elef,f4,time+dtreal,flag)
   vec(:)=vec(:)+(dtreal/6.0_f64)*(f1(:)+2.0_f64*f2(:)+2.0_f64*f3(:)+f4(:))

!   print*, "vec", vec

 end subroutine fulrksolve_noninterpolation


 subroutine fulrkfunc_be_step( &
       vec, &
       elef, &
       f, &
       time, &
       flag)  ! It includes the magnetic field perturbation and the electric field perturbation
   real8,intent(in) :: elef(3)
   real8, dimension(:), intent(in) :: vec
   real8, dimension(:), intent(inout) :: f
   real8, intent(in) :: time
   real8 :: polar(2)
   int4, intent(in) :: flag
   
!   real(8) :: x1(3),v1(3),vperp,rholength,gc(3),time,k
   int4 :: i
   real8 x(3),v(3),mainb(3),minb(3),magf(3),k
 
   x(1)=vec(1)
   x(2)=vec(2)
   x(3)=vec(3)
   v(1)=vec(4)
   v(2)=vec(5)
   v(3)=vec(6)

    if(flag==0) then
      call magfield_2d_interpolation_point(x,magf,flag)
!!$          magf(2)=0.0_f64
!!$          magf(1)=0.0_f64
    end if
   
   f(1)=vec(4)
   f(2)=vec(5)
   f(3)=vec(6)
   f(4)=(elef(1)+v(2)*magf(3)-v(3)*magf(2))
   f(5)=(elef(2)+v(3)*magf(1)-v(1)*magf(3))
   f(6)=(elef(3)+v(1)*magf(2)-v(2)*magf(1))

   return
 end subroutine fulrkfunc_be_step
 
 
end module orbit_storing
