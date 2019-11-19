module gyroaverage_2d
#include "work_precision.h"
  use utilities_module, only: gp_error
  use cartesian_mesh, only: cartesian_mesh_1d
  
  use gyroaverage_utilities, only: s_compute_shape_circle
  
use spline_module, only: s_splcoefper1d0old,&
                         s_compute_splines_coefs_matrix_per_1d, &
                         s_compute_splines_coefs_matrix_nat_1d_new, &
                         s_localize_per_per, &
                         s_contribution_spl, &
                         s_splcoefper1d0old, &
                         s_localize_nat_per_new, &
                         compute_spl2D_double_per_weight, &
                         compute_spl2D_nat_per_weight,  &
                         compute_spl2d_field_point_per_per, &
                         compute_spl2d_field_point_nat_per, &
                         compute_spl2d_firstorder_derivative_point_per_per, &
                         compute_spl2d_firstorder_derivative_point_nat_per, &
                         compute_spl2d_firstorder_derivative_point_polar_nat_per

use field_2d_mesh, only: field_2d_plan 
use constants, only: pi_              
use utilities_module, only: gp_error

   implicit none

   public :: inverse_cubic_splines_matrix , &
             matrix_solve_polar_splines,    &
             initialize_gyroaverage_2d_plan, &
             compute_spl2d_driftsquare_weight_per_per, &
             compute_spl2d_driftsquare_weight_nat_per

 type gyroaverage_2d_plan
     int4  :: N_points    !! number of points on the circle
     real8, dimension(:,:), pointer :: ep_weight_gyro  !! the electrostatic potential after gyro
!  contains
!     procedure :: initialize_gyroaverage_cartesian_2d_plan
!     procedure, pass(gyroaverage) :: compute_gyroaverage_2d=>

  end type gyroaverage_2d_plan
  

 contains

   function initialize_gyroaverage_2d_plan(N_points,mu,m_x1,m_x2,field_2d,geometry,boundary) result(gyro)
 !   subroutine initialize_gyroaverage_2d_plan(gyro,N_points,mu,m_x1,m_x2,field_2d,geometry,boundary)
     type(gyroaverage_2d_plan), pointer :: gyro
     class(field_2d_plan), pointer,intent(in) :: field_2d
     class(cartesian_mesh_1d), pointer,intent(in) :: m_x1, m_x2
     character(len=*), intent(in) :: geometry,boundary
  !   int4, intent(in),optional :: flag
     real8, intent(in) :: mu
     int4, intent(in) :: N_points
     int4 :: i,j, ierr

     allocate(gyro, stat=ierr)
     allocate(gyro%ep_weight_gyro(m_x1%nodes,m_x2%nodes),stat=ierr)
     call gp_error(ierr,"ep_weight_gyro")
     gyro%N_points=N_points

     select case (geometry)
      case ("cartesian")
         select case(boundary)
         case("double_per")
           call compute_gyroweight_2d_double_per(m_x1,m_x2,gyro,field_2d,mu)
 !   print*, "per_per,ep_weight_gyro=", gyro%ep_weight_gyro
  ! print*, "per_per,field_2d%ep_weight_init=", field_2d%ep_weight_init
   case ("nat_per")
           call  compute_gyroweight_2d_nat_per(m_x1,m_x2,gyro,field_2d,mu,geometry)
 !    print*, "nat_per,ep_weight_gyro=", gyro%ep_weight_gyro

         case default
           print*, "boundary is not right,boundary=",boundary 
           stop
         end select
       
      case("polar")

     case default
           print*, "The input geometry is not right, geometry=", geometry
           stop
    end select

   end function initialize_gyroaverage_2d_plan
    
!   subroutine compute_gyro1st_spline_mesh_weight_double_per(field_in,field_weight,m_x1,m_x2)
!     real8, dimension(:,:), intent(in) :: field_in    
!     real8, dimension(:,:), intent(inout) :: field_weight
!     class(cartesian_mesh_1d) :: m_x1,m_x2
!!     class(gyroaverage_cartesian_2d_plan), intent(inout) :: gyro
!!     real8, dimension(m_x1%nodes*m_x2%nodes) :: weight
!     int4 :: i,j
!     comp8 :: fft_array(m_x1%nodes,m_x2%nodes)
! !    comp8, dimension(:,:,:), intent(out) :: mat
!     comp8, dimension(:,:,:), allocatable :: mat_1,mat_2
!     comp8, dimension(:,:,:), allocatable :: D_spl2D
!     comp8, dimension(:,:,:), allocatable :: D_contr
!     real8, dimension(:,:,:), allocatable :: coef
!
!     int4 :: ierr
!     allocate(D_contr(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!     allocate(coef(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!     allocate(D_spl2D(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!     allocate(mat_1(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!     allocate(mat_2(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!
!     call precompute_gyroaverage_double_per_coeff_splines(     &
!          m_x1%num_cells, &
!          m_x2%num_cells, &
!          rho, &
!          points, &
!          N_points, &
!          coef)
!
!     call compute_D_spl2d_per_per( &
!          m_x1%num_cells, &
!          m_x2%num_cells, &
!          D_spl2d)
!
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_3d(:,1,i)
!!     end do       
!
!     call compute_D_contr( &
!          m_x1%nodes, &
!          m_x2%nodes, &
!          coef, &
!          D_contr)
!
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_com(:,1,i)
!!     end do
!     
!     call inverse_cubic_splines_matrix( &
!          coef, &
!          m_x1%nodes-1, &
!          m_x2%nodes)
!     
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_com(:,1,i)
!!     end do
!     
!     call matrix_solve_polar_splines( &
!          field_in, &
!          field_weight, &
!          coef, &
!          m_x1%nodes-1, &
!          m_x2%nodes)
!
!     deallocate(buf_3d)
!     deallocate(coef)
!  
!   end subroutine compute_gyro1st_spline_mesh_weight_double_per


!   subroutine compute_gyro1st_spline_mesh_weight_nat_per(field_in,field_weight,m_x1,m_x2)
!     real8, dimension(:,:), intent(in) :: field_in    
!     real8, dimension(:,:), intent(inout) :: field_weight
!     class(cartesian_mesh_1d) :: m_x1,m_x2
!     
!     int4 :: i,j
!     comp8 :: fft_array(m_x1%nodes,m_x2%nodes)
!     real8, dimension(:), allocatable :: buf_fft
!     real8, dimension(:,:,:), allocatable ::  buf_3d
!     comp8, dimension(:,:,:), allocatable ::  buf_com
!     int4 :: ierr
!
!   !  allocate(buf_fft(4*m_x2%nodes+100), stat=ierr)
!     allocate(buf_3d(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)
!     allocate(buf_com(m_x2%nodes,m_x1%nodes,m_x1%nodes), stat=ierr)    
!!     call zffti(m_x2%nodes,buf_fft)
!!     
!!     do i=1,m_x1%nodes
!!        fft_array(i,:)=field_in(i,:)*(1.0_f64,0.0_f64)
!!        call zfftf(m_x2%nodes,fft_array(i,:),buf_fft)
!!     end do
!     
!     call compute_D_spl2d_nat_per( &
!          m_x1%num_cells, &
!          m_x2%num_cells, &
!          buf_3d)
!
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_3d(:,1,i)
!!     end do
!       
!     call compute_D_contr( &
!          m_x1%nodes, &
!          m_x2%nodes, &
!          buf_3d, &
!          buf_com)
!
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_com(:,1,i)
!!     end do
!     
!     call inverse_cubic_splines_matrix( &
!          buf_com, &
!          m_x1%nodes-1, &
!          m_x2%nodes)
!     
!!     do i=1, m_x1%nodes
!!        print*, i
!!        print*, buf_com(:,1,i)
!!     end do
!     
!     call matrix_solve_polar_splines( &
!          field_in, &
!          field_weight, &
!          buf_com, &
!          m_x1%nodes-1, &
!          m_x2%nodes)
!
!     
!
!!     deallocate(buf_fft)
!     deallocate(buf_3d)
!     deallocate(buf_com)
!   
!
!   end subroutine compute_gyro1st_spline_mesh_weight_nat_per


  function compute_gyroaverage_2d_point_double_per(x1,m_x1,m_x2,gyro,field_2d,mu,gyroaverage) result(pafter)
       type(gyroaverage_2d_plan), pointer, intent(in) ::  gyro
       class(field_2d_plan),pointer,intent(in) :: field_2d
       class(cartesian_mesh_1d), pointer,intent(in) :: m_x1, m_x2
       real8, dimension(:,:), pointer :: points
       real8, intent(in) :: mu
       real8, dimension(2),intent(in) :: x1
       character(len=*),intent(in), optional :: gyroaverage
!       character(len=*),intent(in) :: boundary, geometry
       real8 :: pafter
       real8 :: rho
       int4 :: i,j, k
!       real8 :: points(3,gyro%N_points) 
       real8 :: x(2),eta_min(2),eta_max(2),eta_star(2)
       int4 :: ii(2), Nc(2), flag, ind(2), ell_1, ell_2
       real8 :: val(-1:2,-1:2)

       allocate(points(3,gyro%N_points))
       call s_compute_shape_circle(points,gyro%N_points)
    
       Nc(1) = m_x1%nodes
       Nc(2) = m_x2%nodes
       rho = sqrt(2.0d0*mu)
       eta_min(1)=m_x1%eta_min
       eta_min(2)=m_x2%eta_min
       eta_max(1)=m_x1%eta_max
       eta_max(2)=m_x2%eta_max
   
    pafter=0._f64 
    if(present(gyroaverage)) then
          if(abs(x1(1)-m_x1%eta_min).lt.rho.or. abs(x1(1)-m_x1%eta_max).lt.rho.or.&
            abs(x1(2)-m_x2%eta_min).lt.rho.or. abs(x1(2)-m_x2%eta_max).lt.rho)    then
            pafter=compute_spl2d_field_point_per_per(x1, &
                       gyro%ep_weight_gyro,  &
                       (/m_x1%eta_min,m_x2%eta_min/),    &
                       (/m_x1%eta_max,m_x2%eta_max/),    &
                       (/m_x1%delta_eta,m_x2%delta_eta/),  &
                       (/m_x1%nodes,m_x2%nodes/))
         else                            
       do k=1,gyro%N_points
           flag=1
        x(1) = x1(1)+rho*points(1,k)
        x(2) = x1(2)+rho*points(2,k)
          call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,Nc,flag)
    
          if(flag==1) then
             call s_contribution_spl(eta_star,val)
             val=val*points(3,k)  
             do ell_2=-1,2
                ind(2)=modulo(ii(2)+ell_2,Nc(2))
                  do ell_1=-1,2           
                    ind(1)=modulo(ii(1)+ell_1,Nc(1))   ! this is change to
                    pafter=pafter+val(ell_1,ell_2)*(gyro%ep_weight_gyro(ind(1)+1,ind(2)+1))
                  enddo
             enddo
          end if           
        enddo

       end if

   else 
     
         if(abs(x1(1)-m_x1%eta_min).lt.rho.or. abs(x1(1)-m_x1%eta_max).lt.rho.or.&
!! drift kineitc is implemented at the  boundary 
            abs(x1(2)-m_x2%eta_min).lt.rho.or. abs(x1(2)-m_x2%eta_max).lt.rho)    then
           pafter=compute_spl2d_field_point_per_per(x1, &
                       field_2d%ep_weight_init,  &
                       (/m_x1%eta_min,m_x2%eta_min/),    &
                       (/m_x1%eta_max,m_x2%eta_max/),    &
                       (/m_x1%delta_eta,m_x2%delta_eta/),  &
                       (/m_x1%nodes,m_x2%nodes/))
         else

     do k=1,gyro%N_points
           flag=1
        x(1) = x1(1)+rho*points(1,k)
        x(2) = x1(2)+rho*points(2,k)
          call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,Nc,flag)
    
          if(flag==1) then
             call s_contribution_spl(eta_star,val)
             val=val*points(3,k)  
             do ell_2=-1,2
                ind(2)=modulo(ii(2)+ell_2,Nc(2))
                  do ell_1=-1,2           
                    ind(1)=modulo(ii(1)+ell_1,Nc(1))   ! this is change to
                    pafter=pafter+val(ell_1,ell_2)*(field_2d%ep_weight_init(ind(1)+1,ind(2)+1))
                  enddo
             enddo
          end if
      enddo

      end if

    end if
    
 !   print*, "pafter=",pafter, val(0,0:1)
  end function compute_gyroaverage_2d_point_double_per

  subroutine compute_gyroweight_2d_double_per(m_x1,m_x2,gyro,field_2d,mu)
    type(gyroaverage_2d_plan), pointer,intent(inout) ::  gyro
    class(field_2d_plan),pointer,intent(in) :: field_2d
    class(cartesian_mesh_1d), pointer,intent(in) :: m_x1, m_x2
    real8, intent(in) :: mu
    character(25) :: gyroaverage
    real8, dimension(:,:), pointer :: gyroep
    int4 :: i,j,ierr,xnodes,ynodes
    real8 :: x(2), rho
    allocate(gyroep(m_x1%nodes,m_x2%nodes),stat=ierr)
    xnodes=m_x1%nodes
    ynodes=m_x2%nodes
    do i=1, xnodes
      do j=1, ynodes
         x(1)=m_x1%eta_min+real(i-1,8)*m_x1%delta_eta
         x(2)=m_x2%eta_min+real(j-1,8)*m_x2%delta_eta
            gyroep(i,j)=compute_gyroaverage_2d_point_double_per(x,m_x1,m_x2,gyro,field_2d,mu)
      end do
    end do
!    do i=1,m_x1%nodes
!      print*,"i=",i,gyroep(:,i)
!    enddo
    call compute_spl2D_double_per_weight(gyroep,gyro%ep_weight_gyro,xnodes,ynodes)
  !   print*, "ep_weight_gyro=",gyro%ep_weight_gyro
    deallocate(gyroep)
  end subroutine   

  function compute_gyroaverage_2d_point_nat_per(x1,m_x1,m_x2,gyro,field_2d,mu,geometry,gyroaverage) result(pafter)
       type(gyroaverage_2d_plan),pointer,intent(inout) ::  gyro
       class(cartesian_mesh_1d), pointer,intent(in) :: m_x1, m_x2
       class(field_2d_plan), pointer,intent(in) :: field_2d
       real8, dimension(:),intent(in) :: x1
       real8, intent(in) :: mu
    !   real8, dimension(2),intent(in) :: polar
       character(len=*),intent(in) :: geometry
       character(len=*),intent(in),optional :: gyroaverage
       real8 :: pafter
       real8 :: rho
       int4 :: i,j, k
       real8, dimension(:,:),pointer :: points
!       real8 :: points(3,gyro%N_points)
       real8 :: polar1(2),eta_min(2),eta_max(2),eta_star(2)
       real8 :: x(2)
       int4 :: ii(2), Nc(2), flag, ind(2), ell_1, ell_2
       real8 :: val(-1:2,-1:2)
 
       allocate(points(3,gyro%N_points))     
       pafter=0._f64
       call s_compute_shape_circle(points,gyro%N_points)
    
       Nc(1) = m_x1%nodes
       Nc(2) = m_x2%nodes
       rho = sqrt(2.0d0*mu)
       eta_min(1)=m_x1%eta_min
       eta_min(2)=m_x2%eta_min
       eta_max(1)=m_x1%eta_max
       eta_max(2)=m_x2%eta_max
       
    if(present(gyroaverage)) then 
         if(abs(x1(1)-m_x1%eta_min).lt.rho.or. abs(x1(1)-m_x1%eta_max).lt.rho.or.&
            abs(x1(2)-m_x2%eta_min).lt.rho.or. abs(x1(2)-m_x2%eta_max).lt.rho)    then
            pafter=compute_spl2d_field_point_nat_per(x1, &
                       gyro%ep_weight_gyro,  &
                       (/m_x1%eta_min,m_x2%eta_min/),    &
                       (/m_x1%eta_max,m_x2%eta_max/),    &
                       (/m_x1%delta_eta,m_x2%delta_eta/),  &
                       (/m_x1%nodes,m_x2%nodes/))
         else

       do k=1,gyro%N_points
           flag=1
          polar1(1) =rho*points(1,k)
          polar1(2) =rho*points(2,k)
          if(geometry=="polar") then
            call polar_coordinate_gyro(x,x1,polar1)
          else
            x(1)=x1(1)+polar1(1)
            x(2)=x1(2)+polar1(2)
          end if
          call s_localize_nat_per_new(x,eta_min,eta_max,ii,eta_star,Nc,flag)
    
          if(flag==1) then
           call s_contribution_spl(eta_star,val)
           val=val*points(3,k)  
           do ell_2=-1,2
             ind(2)=modulo(ii(2)+ell_2,Nc(2))
               do ell_1=-1,2           
                 ind(1)=ii(1)+ell_1
                
                 if(ind(1)<0) then
                   ind(1)=0
                 else if(ind(1).ge.(Nc(1)-1)) then
                 ind(1)=Nc(1)-1
                  end if
                 pafter=pafter+val(ell_1,ell_2)*(gyro%ep_weight_gyro(ind(1)+1,ind(2)+1))
               enddo
            enddo
            end if
         enddo

        end if

     else
         if(abs(x1(1)-m_x1%eta_min).lt.rho.or. abs(x1(1)-m_x1%eta_max).lt.rho.or.&
            abs(x1(2)-m_x2%eta_min).lt.rho.or. abs(x(2)-m_x2%eta_max).lt.rho)    then
           pafter=compute_spl2d_field_point_nat_per(x1, &
                       field_2d%ep_weight_init,  &
                       (/m_x1%eta_min,m_x2%eta_min/),    &
                       (/m_x1%eta_max,m_x2%eta_max/),    &
                       (/m_x1%delta_eta,m_x2%delta_eta/),  &
                       (/m_x1%nodes,m_x2%nodes/))
         else

        do k=1,gyro%N_points
          flag=1
          polar1(1) =rho*points(1,k)
          polar1(2) =rho*points(2,k)
          if(geometry=="polar") then
            call polar_coordinate_gyro(x,x1,polar1)
          else
            x(1)=x1(1)+polar1(1)
            x(2)=x1(2)+polar1(2)
          end if
          call s_localize_nat_per_new(x,eta_min,eta_max,ii,eta_star,Nc,flag)
    
          if(flag==1) then
           call s_contribution_spl(eta_star,val)
           val=val*points(3,k)  
           do ell_2=-1,2
             ind(2)=modulo(ii(2)+ell_2,Nc(2))
               do ell_1=-1,2           
                 ind(1)=ii(1)+ell_1
                
                 if(ind(1)<0) then
                   ind(1)=0
                 else if(ind(1).ge.(Nc(1)-1)) then
                 ind(1)=Nc(1)-1
                  end if
                 pafter=pafter+val(ell_1,ell_2)*(field_2d%ep_weight_init(ind(1)+1,ind(2)+1))
               enddo
            enddo
            end if
         enddo
         endif
      end if   

 ! print*, "#coef", coef(:,:,1:3)
  end function compute_gyroaverage_2d_point_nat_per
    
  subroutine compute_gyroweight_2d_nat_per(m_x1,m_x2,gyro,field_2d,mu,geometry)
    type(gyroaverage_2d_plan),pointer,intent(inout) ::  gyro
    class(field_2d_plan),pointer,intent(in) :: field_2d
    class(cartesian_mesh_1d), pointer,intent(in) :: m_x1, m_x2
    real8, intent(in) :: mu
    character(len=*),intent(in) :: geometry
    character(25) :: gyroaverage
    real8, dimension(:,:), pointer :: gyroep
    int4 :: i,j,ierr,xnodes,ynodes
    real8 :: x(2),rho

!    real8 :: eta_min(2),eta_max(2),delta_eta(2)
!    int4  ::  Nc(2)
    allocate(gyroep(m_x1%nodes,m_x2%nodes),stat=ierr)
    xnodes=m_x1%nodes
    ynodes=m_x2%nodes
    rho=sqrt(2._f64*mu)

    do i=1, xnodes
      do j=1, ynodes
         x(1)=m_x1%eta_min+real(i-1,8)*m_x1%delta_eta
         x(2)=m_x2%eta_min+real(j-1,8)*m_x2%delta_eta
          gyroep(i,j)=compute_gyroaverage_2d_point_nat_per(x,m_x1,m_x2,gyro,field_2d,mu,geometry)
      end do
    end do
!  print*, "gyroep=",gyroep
    call compute_spl2D_nat_per_weight(gyroep,gyro%ep_weight_gyro,xnodes-1,ynodes)
!   print*, "nat_per,ep_weight_gyro=",gyro%ep_weight_gyro
    deallocate(gyroep)
  end subroutine   

   subroutine compute_D_spl2D_per_per( &
    num_cells_x, &
    num_cells_y, &
    mat_spl2D)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:,:), intent(inout) :: mat_spl2D
!    real8, dimension(:,:,:),allocatable :: buf_3d
    real8, dimension(:),allocatable,target :: dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
    real8, dimension(:),pointer :: pointer_dper_x,pointer_lper_x,pointer_mper_x
    real8, dimension(:),pointer :: pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8, dimension(:,:),allocatable,target :: mat_per_x,mat_per_y
    real8, dimension(:,:),pointer :: pointer_mat_per_x,pointer_mat_per_y
    int4 :: j
    int4 :: k,i
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
!   do i=0, Nx
!      print*,"mat_per_x", i
!      print*, mat_per_x(:,i)
!   end do
!
    do j=0,Ny-1
       do i=0,Nx-1
          do k=0,Nx-1
             mat_spl2d(j+1,i+1,k+1)=mat_per_y(0,j)*mat_per_x(i,k)
          end do
       end do
!       mat_spl2d(j+1,:,:)=buf_3d(:,:,j+1)
    enddo
    
!     do i=0, Nx
!        print*,"per_per,mat_spl2D", i
!        print*, mat_spl2D(:,i,1)
!     end do

    deallocate(dper_x)
    deallocate(lper_x)
    deallocate(mper_x)
    deallocate(dper_y)
    deallocate(mper_y)
    deallocate(lper_y)
    deallocate(mat_per_x)
    deallocate(mat_per_y)
    
  end subroutine compute_D_spl2D_per_per

!  subroutine compute_D_spl2D_per_per_noblock( &
!    num_cells_x, &
!    num_cells_y, &
!    mat_spl2D)
!    int4, intent(in) :: num_cells_x
!    int4, intent(in) :: num_cells_y
!    int4 :: ierr
!    int4 :: Nx
!    int4 :: Ny
!    real8, dimension(:,:), intent(inout) :: mat_spl2D
!!    real8, dimension(:,:,:),allocatable :: buf_3d
!    real8, dimension(:),allocatable,target :: dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
!    real8, dimension(:),pointer :: pointer_dper_x,pointer_lper_x,pointer_mper_x
!    real8, dimension(:),pointer :: pointer_dper_y,pointer_lper_y,pointer_mper_y
!    real8, dimension(:,:),allocatable,target :: mat_per_x,mat_per_y
!    real8, dimension(:,:),pointer :: pointer_mat_per_x,pointer_mat_per_y
!    int4 :: k,i,h,j
!    real8 :: testm(num_cells_x,num_cells_y),test(4,4)
!        
!    Nx = num_cells_x
!    Ny = num_cells_y
!
!
!    ALLOCATE(dper_x(0:Nx-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(lper_x(0:Nx-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mper_x(0:Nx-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")    
!    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mat_per_x(0:Nx-1,0:Nx-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(buf_3d(1:size(mat_spl2d,3),1:size(mat_spl2d,1),1:size(mat_spl2d,2)),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!
!    dper_x=0.0d0
!    lper_x=0.0d0
!    mper_x=0.0d0
!    dper_y=0.0d0
!    lper_y=0.0d0
!    mper_y=0.0d0
!    mat_per_x=0.0d0
!    mat_per_y=0.0d0
!
!    pointer_dper_y => dper_y
!    pointer_lper_y => lper_y
!    pointer_mper_y => mper_y    
!    pointer_mat_per_y => mat_per_y
! !   pointer_mat_spl2D_circ => mat_spl2D_circ
!
!    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)
!    
!    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
!    
!
!    pointer_dper_x => dper_x
!    pointer_lper_x => lper_x
!    pointer_mper_x => mper_x    
!    pointer_mat_per_x => mat_per_x
!    
!    call s_splcoefper1d0old(dper_x,lper_x,mper_x,Nx)
!    
!    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_x,pointer_dper_x,pointer_lper_x,pointer_mper_x,Nx)
!
! do h=0,Ny-1
!    do j=0,Ny-1
!       do i=0,Nx-1
!          do k=0,Nx-1
!             mat_spl2d(h*Nx+i,j*Nx+k)=mat_per_y(h,j)*mat_per_x(i,k)
!          end do
!       end do
!    enddo
! end do   
!
!    deallocate(dper_x)
!    deallocate(lper_x)
!    deallocate(mper_x)
!    deallocate(dper_y)
!    deallocate(mper_y)
!    deallocate(lper_y)
!    deallocate(mat_per_x)
!    deallocate(mat_per_y)
!    
!  end subroutine compute_D_spl2D_per_per_noblock


  subroutine compute_D_spl2D_nat_per( &
    num_cells_x, &
    num_cells_y, &
    mat_spl2D)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    real8, dimension(:,:,:), intent(inout) :: mat_spl2D
    real8,dimension(:),allocatable,target::dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
    real8,dimension(:),pointer::pointer_dper_x,pointer_lper_x,pointer_mper_x
    real8,dimension(:),pointer::pointer_dper_y,pointer_lper_y,pointer_mper_y
    real8,dimension(:,:),allocatable,target:: mat_nat_x,mat_per_y
    real8,dimension(:,:),pointer::pointer_mat_nat_x,pointer_mat_per_y
!    sll_real64,dimension(:,:,:),allocatable,targemat_spl2D
!    sll_real64,dimension(:,:,:),pointer::pointer_mat_spl2D_circ
    int4 :: k,i,j    
    
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
    pointer_dper_y => dper_y
    pointer_lper_y => lper_y
    pointer_mper_y => mper_y    
    pointer_mat_per_y => mat_per_y

    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)    
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
 
    pointer_mat_nat_x=> mat_nat_x
    
    call s_compute_splines_coefs_matrix_nat_1d_new(pointer_mat_nat_x,Nx+1)
!     do i=0, Nx
!        print*,"mat_nat_x", i
!        print*, mat_nat_x(:,i)
!     end do
    do j=0,Ny-1
       do i=0,Nx
          do k=0,Nx
             mat_spl2D(i+1,k+1,j+1)=mat_per_y(0,j)*mat_nat_x(i,k)
          end do
       end do
    enddo
!     do i=0, Nx
!        print*,"nat_per,mat_spl2D", i
!        print*, mat_spl2D(:,i,1)
!     end do

    deallocate(mat_nat_x)
    deallocate(mat_per_y)
    deallocate(dper_y)
    deallocate(mper_y)
    deallocate(lper_y)
    

  end subroutine compute_D_spl2D_nat_per


!  subroutine compute_D_spl2D_nat_per_noblock( &
!    num_cells_x, &
!    num_cells_y, &
!    mat_spl2D)
!    int4, intent(in) :: num_cells_x
!    int4, intent(in) :: num_cells_y
!    int4 :: ierr
!    int4 :: Nx
!    int4 :: Ny
!    real8, dimension(:,:), intent(inout) :: mat_spl2D
!    real8,dimension(:),allocatable,target::dper_x,lper_x,mper_x,dper_y,lper_y,mper_y
!    real8,dimension(:),pointer::pointer_dper_x,pointer_lper_x,pointer_mper_x
!    real8,dimension(:),pointer::pointer_dper_y,pointer_lper_y,pointer_mper_y
!    real8,dimension(:,:),allocatable,target:: mat_nat_x,mat_per_y
!    real8,dimension(:,:),pointer::pointer_mat_nat_x,pointer_mat_per_y
!    int4 :: k,i,h,j     
!    
!    Nx = num_cells_x
!    Ny = num_cells_y
!   
!    ALLOCATE(dper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(lper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mper_y(0:Ny-1),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mat_nat_x(0:Nx,0:Nx),stat=ierr)
!call gp_error(ierr,"gyroaverage_2d")
!    ALLOCATE(mat_per_y(0:Ny-1,0:Ny-1),stat=ierr)    
!call gp_error(ierr,"gyroaverage_2d")
!    pointer_dper_y => dper_y
!    pointer_lper_y => lper_y
!    pointer_mper_y => mper_y    
!    pointer_mat_per_y => mat_per_y
!
!    call s_splcoefper1d0old(dper_y,lper_y,mper_y,Ny)    
!    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per_y,pointer_dper_y,pointer_lper_y,pointer_mper_y,Ny)
! 
!    pointer_mat_nat_x=> mat_nat_x
!    
!    call s_compute_splines_coefs_matrix_nat_1d_new(pointer_mat_nat_x,Nx+1)
!
! do h=0, Ny-1 
!    do j=0,Ny-1
!       do i=0,Nx
!          do k=0,Nx
!             mat_spl2D(h*(Nx+1)+i,j*(Nx+1)+k)=mat_per_y(h,j)*mat_nat_x(i,k)
!          end do
!       end do
!    enddo
! end do
!
!    deallocate(mat_nat_x)
!    deallocate(mat_per_y)
!    deallocate(dper_y)
!    deallocate(mper_y)
!    deallocate(lper_y)
!    
!
!  end subroutine compute_D_spl2D_nat_per_noblock


  subroutine compute_D_contr( &
    nodes_x, &
    nodes_y, &
    coef, &
    D_contr)
    int4, intent(in) :: nodes_x
    int4, intent(in) :: nodes_y
    real8, dimension(:,:,:), intent(in) :: coef
    int4 :: ierr
    int4 :: Nx
    int4 :: Ny
    comp8, dimension(:,:,:), intent(out) :: D_contr
    int4 :: j, i
    comp8, dimension(:), allocatable :: fft_array
    real8, dimension(:), allocatable :: buf_fft
    int4 :: k

    
    Nx = nodes_x
    Ny= nodes_y

    ALLOCATE(buf_fft(4*Ny+100),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(fft_array(Ny),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")


    call zffti(Ny,buf_fft)

    do k=0,Nx-1
      do i=0,Nx-1
        fft_array(1:Ny)=coef(1:Ny,i+1,k+1)*(1.0_f64,0.0_f64)
        call zfftf(Ny,fft_array(1:Ny),buf_fft)
        D_contr(1:Ny,i+1,k+1) = fft_array(1:Ny)
      enddo  
    enddo   

  end subroutine compute_D_contr  
   

 subroutine inverse_cubic_splines_matrix(&
    mat, &    
    num_cells_r, &
    num_cells_theta)
   
    comp8, dimension(:,:,:), intent(inout) :: mat
    int4, intent(in) :: num_cells_r
    int4, intent(in) :: num_cells_theta
    int4 :: ierr
    int4, dimension(:), allocatable :: IPIV
    comp8, dimension(:), allocatable :: WORK
    int4 :: INFO
    int4 :: i,m, j


  
    ALLOCATE(IPIV(num_cells_r+1),stat=ierr)
call gp_error(ierr,"gyroaverage_2d")
    ALLOCATE(WORK((num_cells_r+1)**2),stat=ierr) 
call gp_error(ierr,"gyroaverage_2d")


    do m=1,num_cells_theta
    call ZGETRF( &            
    num_cells_r+1, &          
    num_cells_r+1, &
    mat(m,:,:), &
    num_cells_r+1, &
    IPIV, &
    INFO)

   print*, "##", 3, m
    
    call ZGETRI( &
    num_cells_r+1, &
    mat(m,:,:), &
    num_cells_r+1, &
    IPIV, &
    WORK, &
    (num_cells_r+1)**2, &
    INFO)

   print*, "##", 4, m
    
    enddo

!    print *,'#here is INFO'
!    print *,INFO

    deallocate(IPIV)
    deallocate(work)   

  end subroutine inverse_cubic_splines_matrix



  subroutine matrix_solve_polar_splines( &    ! this subroutine can be used to solve the quasi-neutral as the last step
    phi_in, &
    phi_out, &   
    mat, &   
    num_cells_r, &
    num_cells_theta)

    comp8, dimension(:,:,:), intent(in) :: mat
    real8, dimension(:,:), intent(in) :: phi_in
    real8, dimension(:,:), intent(out) :: phi_out
    comp8,dimension(:,:),allocatable :: phi_comp
    comp8,dimension(:,:),allocatable :: phi_old
    real8,dimension(:),allocatable::buf_fft
    int4, intent(in) :: num_cells_r
    int4, intent(in) :: num_cells_theta

    int4 :: Nr
    int4 :: Ntheta
    int4 :: m
    int4 :: i
    int4 :: j
    int4 :: ierr
    comp8 :: result


!!!!!!!!!!!!!!!!!!!!
  
    Nr = num_cells_r
    Ntheta = num_cells_theta
    
    ALLOCATE(phi_comp(1:Nr+1,1:Ntheta),stat=ierr)
call gp_error(ierr,"phi_comp")
    ALLOCATE(phi_old(1:Nr+1,1:Ntheta),stat=ierr)
call gp_error(ierr,"phi_old")
    ALLOCATE(buf_fft(1:4*Ntheta+100),stat=ierr)
call gp_error(ierr,"buf_fft") 
 ! FFT(PHI)
    phi_comp=phi_in*(1.0d0,0.0d0)
    call zffti(Ntheta,buf_fft)
    do i=1,Nr+1
    call zfftf(Ntheta,phi_comp(i,:),buf_fft)
    enddo   
print*, 11
    phi_old=phi_comp
    do m = 1,Ntheta
      do i = 1,Nr+1
        result = (0.0d0,0.0d0)
        do j = 1,Nr+1
          result = result + mat(m,i,j)*phi_old(j,m)   ! notice the indexes i and j;it's the inverse.
        enddo
        phi_comp(i,m) = result
      enddo 
    enddo
 print*, 12   
 ! FFT^-1
    do i=1,Nr+1
    call zfftb(Ntheta,phi_comp(i,:),buf_fft)
    enddo
print*, 13
    phi_out=real(phi_comp,8)/real(Ntheta,8)
    
 ! Sorties     
 !   print *,"phi_min : ",minval(phi)
 !   print *,"phi_max : ",maxval(phi)
    
 ! Deallocate    
    DEALLOCATE(phi_comp,stat=ierr)
    DEALLOCATE(phi_old,stat=ierr)
    DEALLOCATE(buf_fft,stat=ierr)
        
  end subroutine matrix_solve_polar_splines


!  subroutine polar_coordinate_gyro(x,x1,x2)
!    real8, intent(out) :: x(2)
!    real8, intent(in ) :: x1(2),x2(2)
!  !  int4, optional  :: gyroorder
!    real8 :: y1,y2   
!    
!    y1=x1(1)*cos(x1(2))+x2(1)
!    y2=x1(1)*sin(x1(2))+x2(2)
!    x(1)=sqrt(y1**2+y2**2)
!
!    x(2)=atan2(y2,y1)
!
!end subroutine polar_coordinate_gyro  
 

 subroutine precompute_gyroaverage_double_per_coeff_splines( &
    num_cells_x, &
    num_cells_y, &
    rho, &
    points, &
    N_points, &
    coef)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    real8, dimension(:),intent(in) :: rho
    real8, dimension(:,:), intent(in) :: points
    int4, intent(in) :: N_points
    real8, dimension(:,:,:),intent(out) :: coef
    int4 ::i,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    real8::eta_star(2),eta(2),delta_eta(2),x(2),x1(2)
    int4 ::error,max_nb
!    sll_int32, intent(in) :: rank
    real8 :: eta_min(2), eta_max(2)
    int4 :: Nc(2)
    int4 :: gyroorder=1
    int4 :: flag
    real8 :: val(-1:2,-1:2)
 
    delta_eta(1)=2._f64*pi_/real(num_cells_x,f64)
    delta_eta(2)=2._f64*pi_/real(num_cells_y,f64)

    Nc(1) = num_cells_x
    Nc(2) = num_cells_y

    eta_min(1) = 0._f64
    eta_min(2) = 0._f64
    eta_max(1) = 2._f64*pi_
    eta_max(2) = 2._f64*pi_

    eta(2)=0._f64
    max_nb=0
    nb=0
    do i=1,Nc(1)
       x1(1)=eta_min(1)+real(i-1,f64)*delta_eta(1)
       x1(2)=0.0_f64
       do k=1,N_points
           flag=1
        x(1) = x1(1)+rho(i)*points(1,k)
        x(2) = x1(2)+rho(i)*points(2,k)
         call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,Nc,flag)

        if(flag==1) then
           call s_contribution_spl(eta_star,val)
           val=val*points(3,k)
           do ell_2=-1,2
             ind(2)=modulo(ii(2)+ell_2,Nc(2))
               do ell_1=-1,2
                 ind(1)=modulo(ii(1)+ell_1,Nc(1))   ! this is change to
                 coef(i,ind(1)+1,ind(2)+1)=coef(i,ind(1)+1,ind(2)+1)+val(ell_1,ell_2)
               enddo
           enddo
        end if
       end do
    enddo    

 end subroutine precompute_gyroaverage_double_per_coeff_splines

subroutine precompute_gyroaverage_nat_per_cartesian_splines( &
    num_cells_x, &
    num_cells_y, &
    rho, &
    points, &
    N_points, &
    coef)
    int4, intent(in) :: num_cells_x
    int4, intent(in) :: num_cells_y
    real8, dimension(:),intent(in) :: rho
    real8, dimension(:,:), intent(in) :: points
    int4, intent(in) :: N_points
    real8, dimension(:,:,:),intent(inout) :: coef
    int4 ::i,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    real8::eta_star(2),eta(2),delta_eta(2),x(2),x1(2)
    int4 ::error,max_nb
    real8 :: eta_min(2), eta_max(2)
    int4 :: Nc(2)
    int4 :: gyroorder=1
    int4 :: flag
    real8 :: val(-1:2,-1:2)
 
    delta_eta(1)=2._f64*pi_/real(num_cells_x,f64)
    delta_eta(2)=2._f64*pi_/real(num_cells_y,f64)

    Nc(1) = num_cells_x
    Nc(2) = num_cells_y

    eta_min(1) = 0._f64
    eta_min(2) = 0._f64
    eta_max(1) = 2._f64*pi_
    eta_max(2) = 2._f64*pi_

    eta(2)=0._f64
    max_nb=0
    nb=0
    do i=1,Nc(1)+1
       x1(1)=eta_min(1)+real(i-1,f64)*delta_eta(1)
       x1(2)=0.0_f64

       do k=1,N_points
           flag=1
        x(1) = x1(1)+rho(i)*points(1,k)
        x(2) = x1(2)+rho(i)*points(2,k)
        call s_localize_nat_per_new(x,eta_min,eta_max,ii,eta_star,Nc,flag)
       if(flag==1) then
           call s_contribution_spl(eta_star,val)
          val=val*points(3,k)
        do ell_2=-1,2
           ind(2)=modulo(ii(2)+ell_2,Nc(2))
             do ell_1=-1,2
                ind(1)=ii(1)+ell_1   ! this is change to
             if(ind(1)<0) then
                ind(1)=0
             else if(ind(1).ge.num_cells_x) then
                ind(1)=num_cells_x
             end if
               coef(i,ind(1)+1,ind(2)+1)=coef(i,ind(1)+1,ind(2)+1)+val(ell_1,ell_2)

             enddo
         enddo
      end if

    enddo

 enddo
 ! print*, "#coef", coef(:,:,1:3)

end subroutine precompute_gyroaverage_nat_per_cartesian_splines         


  subroutine polar_coordinate_gyro(x,x1,x2,gyroorder)
    real8, intent(out) :: x(2)
    real8, intent(in ) :: x1(2),x2(2)
    int4, optional  :: gyroorder
    real8 :: y1,y2

    if(.not.present(gyroorder)) then
       print*, 'gyroorder is not input'
    else if(gyroorder==1) then

    y1=x1(1)*cos(x1(2))+x2(1)
    y2=x1(1)*sin(x1(2))+x2(2)
    x(1)=sqrt(y1**2+y2**2)

    x(2)=atan2(y2,y1)

  else if(gyroorder==2) then

    y1=x1(1)*cos(x1(2))-x2(1)
    y2=x1(1)*sin(x1(2))-x2(2)
    x(1)=sqrt(y1**2+y2**2)

    x(2)=atan2(y2,y1)

   end if
  end subroutine polar_coordinate_gyro

  subroutine compute_spl2d_driftsquare_weight_per_per(driftsquare_weight,bamp,ep_weight_gyro,eta_min,eta_max,eta_delta,NC) 
    real8, intent(in) :: eta_min(2),eta_max(2),eta_delta(2)
    real8, dimension(:,:),pointer, intent(inout) :: driftsquare_weight
    real8, dimension(:,:), intent(in) :: bamp
    real8, dimension(:,:),pointer, intent(in) :: ep_weight_gyro
    int4, intent(in) :: NC(2)
    real8, dimension(:,:),pointer :: driftsquare
    real8 :: deri_firstorder(2)
    real8 :: y(2)
    int4 :: i,j,ierr
    
    allocate(driftsquare(0:NC(1)-1,0:NC(2)-1),stat=ierr)
    call gp_error(ierr,"driftsquare")    
    do i=0,Nc(1)-1
       do j=0,Nc(2)-1
          y(1)=eta_delta(1)*real(i,8)
          y(2)=eta_delta(2)*real(j,8)
          deri_firstorder=0.0_f64
          call compute_spl2d_firstorder_derivative_point_per_per(y, &
              ep_weight_gyro, &
              eta_min,   &
              eta_max,   &
              eta_delta, &
              NC,        &
              deri_firstorder)
         driftsquare(i,j)=abs(deri_firstorder(1)/bamp(i,j))**2+abs(deri_firstorder(2)/bamp(i,j))**2          
       end do
    end do

!  open(10,file="/PARA/blsc950/electrostatic_exp/run/test.txt",status="replace")
!    do, j=Nc(2)-19,Nc(2)-1
!      do i=0, Nc(1)-1 
!       write(10,"(1f16.10,2I8)"), ep_weight_gyro(i,j),j,i
!      end do    
!    end do
!  close(10)
!
    call compute_spl2d_double_per_weight(&
         driftsquare,driftsquare_weight,NC(1),NC(2))

    end subroutine compute_spl2d_driftsquare_weight_per_per
 
  subroutine compute_spl2d_driftsquare_weight_nat_per(driftsquare_weight,bamp,ep_weight_gyro,&
                                                      eta_min,eta_max,eta_delta,NC, geometry) 
    real8, intent(in) :: eta_min(2),eta_max(2),eta_delta(2)
    real8, dimension(:,:),pointer, intent(inout) :: driftsquare_weight
    real8, dimension(:,:), intent(in) :: bamp
    real8, dimension(:,:),pointer, intent(in) :: ep_weight_gyro
    int4, intent(in) :: NC(2)
    character(len=*),intent(in) :: geometry
    real8, dimension(:,:),pointer :: driftsquare
    real8 :: deri_firstorder(2)
    real8 :: y(2)
    int4 :: i,j,ierr
    
    allocate(driftsquare(0:NC(1)-1,0:NC(2)-1),stat=ierr)
    call gp_error(ierr,"driftsquare")    
    do i=0,Nc(1)-1
       do j=0,Nc(2)-1
          y(1)=eta_delta(1)*real(i,8)
          y(2)=eta_delta(2)*real(j,8)
          deri_firstorder=0.0_f64
          select case(geometry)
          case("cartesian")
          call compute_spl2d_firstorder_derivative_point_nat_per(y, &
              ep_weight_gyro, &
              eta_min,   &
              eta_max,   &
              eta_delta, &
              NC,        &
              deri_firstorder)
          case ("polar")
          call compute_spl2d_firstorder_derivative_point_polar_nat_per(y, &
              ep_weight_gyro, &
              eta_min,   &
              eta_max,   &
              eta_delta, &
              NC,        &
              deri_firstorder)

          case default
            stop
          end select
         driftsquare(i,j)=abs(deri_firstorder(1)/bamp(i,j))**2+abs(deri_firstorder(2)/bamp(i,j))**2          
       end do
    end do

!  open(10,file="/PARA/blsc950/electrostatic_exp/run/test.txt",status="replace")
!    do, j=Nc(2)-19,Nc(2)-1
!      do i=0, Nc(1)-1 
!       write(10,"(1f16.10,2I8)"), ep_weight_gyro(i,j),j,i
!      end do    
!    end do
!  close(10)
!
    call compute_spl2d_nat_per_weight(&
         driftsquare,driftsquare_weight,NC(1),NC(2))
  end subroutine compute_spl2d_driftsquare_weight_nat_per
   

end module gyroaverage_2d
