module field_initialize
#include "work_precision.h"
  use constants, only: pi_
  use field_2d_mesh, only: field_2d_plan
  use cartesian_mesh, only: cartesian_mesh_1d
  use paradata_type, only: pic_para_total2d_base
  use paradata_utilities, only: get_coords_from_processrank, &
                                coords_from_localind, &
                                copy_boundary_value_per_per

  implicit none
  include "mpif.h"
 
  public :: epotential_2d, &
       epotential_2d_analytical, &
       efield2d_anal, &
       compute_field_2d_mesh, &
       magfield_2d_interpolation_point, &
       para_initialize_field_2d_mesh

contains

   function epotential_2d_analytical(wapone, &
                                    wavenumber_one,  &
                                    wavenumber_two,  &
                                    xmin, &
                                    xmax, &
                                    x) result(epo)   ! the fist kind of electric potential. k is the wave vector needing to be scanned
  real8, intent(in) :: x(2)
  real8 :: epo
  real8, intent(in) :: wapone, wavenumber_one, wavenumber_two 
  real8, intent(in) :: xmin(2),xmax(2)
  real8 :: xmid,y
!  call wavearray(input)
 !    epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
 !    efi(2)=input(2,2)*sin(input(1,2)*x(1))*cos(input(1,2)*x(2)) 
!     epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
! epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))
!<<<<<<< HEAD
!  epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+waveamplitude*(x(1))/4.0
!=======
!  epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+waveamplitude*(x(2)+x(1))/4.0
!>>>>>>> solverk4
!       epo=waveamplitude*x(2)
  !    epo=waveamplitude

   xmid=(xmax(2)+xmin(2))/2.0_F64
   if(x(2).le.xmid) then
     y=x(2)-xmin(2)
     epo=wapone*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+wapone*y/4.0
   else 
     y=xmid-(x(2)-xmid)
     epo=wapone*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+wapone*y/4.0
   end if 

     return
   end function epotential_2d_analytical


   function epotential_2d(wapone,wapeq, &
                                    wavenumber_one,  &
                                    wavenumber_two,  &
                                    xmin, &
                                    xmax, &
                                    x,epkind) result(epo)   ! the fist kind of electric potential. k is the wave vector needing to be scanned
        real8, intent(in) :: x(2)
        character(len=*) :: epkind
        real8 :: epo
        real8, intent(in) :: wapone,wapeq, wavenumber_one, wavenumber_two 
        real8, intent(in) :: xmin(2),xmax(2)
        real8 :: xmid,y
        !  call wavearray(input)
        !    epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
        !    efi(2)=input(2,2)*sin(input(1,2)*x(1))*cos(input(1,2)*x(2)) 
        !     epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
! epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))
        !<<<<<<< HEAD
        !  epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+waveamplitude*(x(1))/4.0
        !=======
        !  epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+waveamplitude*(x(2)+x(1))/4.0
        !>>>>>>> solverk4
!       epo=waveamplitude*x(2)
        !    epo=waveamplitude
        if(epkind=="eqp_linear")  then
          xmid=(xmax(2)+xmin(2))/2.0_F64
            if(x(2).le.xmid) then
              y=x(2)-xmin(2)
              epo=wapone*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+wapeq*y/4.0
            else 
              y=xmid-(x(2)-xmid)
              epo=wapone*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+wapeq*y/4.0
            end if 
        else if(epkind=="eqp_sin") then
          xmid=(xmax(2)+xmin(2))/2.0_F64
          epo=wapone*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+wapeq*sin(pi_/xmid*(x(2)-xmid))
        end if

        return
   end function epotential_2d




   subroutine efield2d_anal(efi,waveamplitude, &
                                    wavenumber_one,  &
                                    wavenumber_two,  &
                                    x)   ! the fist kind of electric potential. k is the wave vector needing to be scanned
  real8, intent(in) :: x(2)
  real8, intent(inout) :: efi(3)
  real8, intent(in) :: waveamplitude, wavenumber_one, wavenumber_two
!  call wavearray(input)
 !    epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
 !    efi(2)=input(2,2)*sin(input(1,2)*x(1))*cos(input(1,2)*x(2)) 
 !    epo=waveamplitude*cos(wavenumber_two*x(2))
!  efi(2)=-waveamplitude
  efi(2)=-waveamplitude*sin(wavenumber_two*x(2))
  efi(1)=0.0_f64
  efi(3)=0.0_f64
    return
  end subroutine efield2d_anal


  function Bfield_third(contr,x)
    real8, intent(in) :: contr
    real8 :: Bfield_third
    real8  x(2)

    Bfield_third=1.0_f64 +1.0_f64/(1.0_f64+contr)

  end function Bfield_third


  subroutine compute_field_2d_mesh( &
      field_2d, &
      m_x1,  &
      m_x2,  &
       amp,&
       wave_one, &
       wave_two, &
       contr,    &
      geometry)
      class(field_2d_plan), pointer,intent(inout) :: field_2d 
      class(cartesian_mesh_1d), pointer,intent(in) :: m_x1,m_x2
      real8, intent(in) :: amp,wave_one,wave_two
      real8, intent(in) :: contr
 !     character(len=*), intent(in) :: boundary
      character(len=*), intent(in) :: geometry
      real8 :: x(2),x1(2),xmin(2),xmax(2)
      int4 :: i, j

    xmin(1:2)=(/m_x1%eta_min,m_x2%eta_min/)
    xmax(1:2)=(/m_x1%eta_max,m_x2%eta_max/)
    select case (geometry)
      case ("cartesian")
         
!!$      select case (boundary)  ! for double_per, Bf=1.0
!!$      case ("double_per")    
       x=0.0_f64
       do i=0,m_x1%nodes-1
         x(1)=x(1)+m_x1%delta_eta
         x(2)=0.0_f64
          do j=0,m_x2%nodes-1
            x(2)=x(2)+m_x2%delta_eta
!           field_2d%epotential_init(i+1,j+1) = epotential_2d_analytical(amp,wave_one, &
!                wave_two,xmin,xmax, x)
           field_2d%epotential_init(i+1,j+1) = 1.0
 

          field_2d%Bf_init_3rd(i+1,j+1)=1.0
          end do
       end do
!print*, field_2d%epotential_init(10,:)
!!$    case ('nat_per')
!!$        do i=0,m_x1%nodes
!!$         x(1)=x(1)+m_x1%delta_eta
!!$         x(2)=0.0_f64
!!$          do j=0,m_x2%nodes-1
!!$            x(2)=x(2)+m_x2%delta_eta
!!$            field_2d%epotential_init(i+1,j+1) = epotential_2d_analytical(amp,wave_one, &
!!$                 wave_two, x)
!!$            field_2d%Bf_init_3rd(i+1,j+1)=Bfield_third(contr,x)
!!$          end do
!!$         end do       
!!$
!!$       case default
!!$          print*, "input the correct boudnary, boundary=", boundary
!!$         stop  
!!$       end select
 
    case ("polar")
        do i=0,m_x1%nodes
          x(1)=x(1)+m_x1%delta_eta
          x(2)=0.0_f64
            do j=0,m_x2%nodes-1
               x(2)=x(2)+m_x2%delta_eta
               x1(1)=x(1)*cos(x(2))
               x1(2)=x(2)*sin(x(2))
               field_2d%epotential_init(i+1,j+1) = epotential_2d_analytical(amp,wave_one, &
                    wave_two,xmin,xmax, x1)
               field_2d%Bf_init_3rd(i+1,j+1)=Bfield_third(contr,x1)
            end do
         end do
         
    case default
       print*, "input the correct geometry, geometry=", geometry
       stop
    end select

  end subroutine compute_field_2d_mesh

  subroutine magfield_2d_interpolation_point(x,mag_field,flag) !=> For flags.ne.flat, it needs the interpolation
    real8, dimension(2),intent(in) :: x               !=> of the magnetic field 
    real8, dimension(3),intent(inout) :: mag_field
    int4, intent(in) :: flag

    if(flag==0) then
       mag_field(1)=0.0_f64
       mag_field(2)=0._f64
       mag_field(3)=1._f64
  

    end if
    
  end subroutine magfield_2d_interpolation_point  
  

  subroutine para_initialize_field_2d_mesh(amp,amp_eq,wave_one, wave_two, pic2d)
      class(pic_para_total2d_base), pointer,intent(inout) :: pic2d
      real8, intent(in) :: amp,amp_eq,wave_one,wave_two
!      real8, intent(in) :: contr
      character(25) :: geometry,boundary
      real8 :: x(2),x1(2),delta(2),xmin(2),xmax(2)
      int4 :: ND(2),rank,global_sz(2),comm,coords(2),numproc(2)
      int4 :: i, j
  
      character(25) :: epkind="eqp_sin"

      comm=pic2d%layout2d%collective%comm
      numproc=pic2d%para2d%numproc
      coords=get_coords_from_processrank(rank,numproc)
      rank=pic2d%layout2d%collective%rank     
      geometry=pic2d%para2d%geometry
      boundary=pic2d%para2d%boundary
      global_sz(1:2)=(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/)
      ND(1)=pic2d%para2d%m_x1%nodes
      ND(2)=pic2d%para2d%m_x2%nodes  
      delta(1)=pic2d%para2d%m_x1%delta_eta
      delta(2)=pic2d%para2d%m_x2%delta_eta
      xmin=pic2d%para2d%gxmin
      xmax=pic2d%para2d%gxmax
         
    select case (geometry)
      case ("cartesian")
        select case (boundary)  ! for double_per, Bf=1.0
          case ("double_per")   
            if(coords(1)==numproc(1)-1.and.coords(2).ne.numproc(2)-1) then
              do j=1, ND(2)
                 do i=1, ND(1)-1
                    x(1:2)=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)
                    pic2d%field2d%ep(i,j)=epotential_2d(amp,amp_eq,wave_one,wave_two,xmin,xmax, x,epkind)
                 end do
             end do              
            else if(coords(1).ne.numproc(1)-1.and.coords(2)==numproc(2)-1) then
              do j=1, ND(2)-1
                 do i=1, ND(1)
                    x(1:2)=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)
                    pic2d%field2d%ep(i,j)=epotential_2d(amp,amp_eq,wave_one,wave_two,xmin,xmax,x,epkind)              
                 end do
              end do
            else if(coords(1)==numproc(1)-1.and.coords(2)==numproc(2)-1)   then
               do j=1,ND(2)-1
                  do i=1,Nd(1)-1
                    x(1:2)=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)
                    pic2d%field2d%ep(i,j)=epotential_2d(amp,amp_eq,wave_one,wave_two,xmin,xmax,x,epkind)
                  end do
               end do
            else
               do j=1,ND(2)
                  do i=1,Nd(1)
                    x(1:2)=coords_from_localind((/i,j/),rank,pic2d%para2d%gboxmin,delta)
                    pic2d%field2d%ep(i,j)=epotential_2d(amp,amp_eq,wave_one,wave_two,xmin,xmax,x,epkind)
                  end do
               end do
            end if                           

         call copy_boundary_value_per_per(pic2d%field2d%ep,rank,pic2d%para2d%numproc,pic2d%layout2d)
        
            do j=1,Nd(2)+1
               do i=1, Nd(1)+1
                  pic2d%field2d%bf03(i,j)=1.0_f64
               enddo
            end do
 
         case("nat_per")
            print*, "error: The current version doesn't include the nat_per boundary"
            stop
         case default
            stop
         end select

    case ("polar")
        print*, "error: The current version doesn't include the polar geometry"        
        stop
    case default
       print*, "input the correct geometry, geometry=", geometry
       stop
    end select

  end subroutine para_initialize_field_2d_mesh




 

end module field_initialize
