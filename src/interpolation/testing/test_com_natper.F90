program test_com_natper
#include "work_precision.h"
use spline_module, only: &
     s_compute_splines_coefs_matrix_nat_1d_new, &
     s_splcoefper1d0old, &
     s_compute_splines_coefs_matrix_per_1d, &
     compute_spl2d_field_point_nat_1d, &
     compute_spl2d_field_point_per_1d 
use cartesian_mesh, only: cartesian_mesh_1d

implicit none

      integer :: N=10,pi=3.1415926
      real8 :: eta_min,eta_max,delta_eta
      real8, dimension(:),pointer :: dper,lper,mper
      real8, dimension(:,:),pointer :: matnat, matper
      real8, dimension(:,:), pointer :: weightnat,weightper
      real8, dimension(:),allocatable :: epo
      real8, dimension(:,:),pointer :: epoone
      real8, dimension(:,:), pointer :: pos_epo
      class(cartesian_mesh_1d),pointer :: m_x 
      real8 :: x
      integer :: i,j,m
      real8 :: errper=0._f64,errnat=0._f64

      allocate(m_x)
      allocate(dper(N),lper(N),mper(N),matnat(N,N),matper(N,N),weightnat(N,1),weightper(N,1),epo(N))
      allocate(epoone(N,1))
     ! allocate(pos_epo(N,2))
      allocate(pos_epo(N,3))
      eta_min=0.0
      eta_max=2.0*pi
      delta_eta=(eta_max-eta_min)/real(N,8)

      m_x%eta_min=0.0
      m_x%eta_max=eta_max
      m_x%delta_eta=delta_eta

      !!!!=++++initialize the position of a series of points 
      do i=1, N
         pos_epo(i,1)=m_x%eta_min+(real(i-1,8)+0.3)*m_x%delta_eta
      end do
      call field_initialize(epo,N,m_x)
      epoone(:,1)=epo(:)
      call s_compute_splines_coefs_matrix_nat_1d_new(matnat,N)
      
      do i=1,N
         print*, "i=",i,"matnat(i,:)=",matnat(i,:)
      end do

      weightnat=matmul(matnat,epoone)

      do i=1, N
         x=pos_epo(i,1)         
         pos_epo(i,2)=compute_spl2d_field_point_nat_1d(x,weightnat(:,1),eta_min,eta_max,delta_eta,N)          
      end do
      
      call s_splcoefper1d0old(dper,lper,mper,N)
      call s_compute_splines_coefs_matrix_per_1d(matper,dper,lper,mper,N) 
      
      weightper=matmul(matper, epoone)
      do i=1,N
        x=pos_epo(i,1)
        pos_epo(i,3)=compute_spl2d_field_point_per_1d(x,weightper(:,1),eta_min,eta_max,delta_eta,N)
      end do

!      do i=1, N
!         print*, i,weightper(i,1),weightnat(i,1)
!      end do
!
!      print*,  "++++++++++++++++"
!
!      do i=1, N
!         print*,i, weightper(i,1)-weightnat(i,1)
!      end do

      open(50,file="/Users/zhang/program/debug/electrostatic_exp/run/com_natper.txt",status="replace")
       do i=1,N
          write(50, "(I3,3F19.10)") i,pos_epo(i,2),pos_epo(i,3),epo(i)
       end do
      close(50)
     
      do i=1, N
        errper=errper+(pos_epo(i,3)-epo(i))**2*delta_eta
        errnat=errnat+(pos_epo(i,2)-epo(i))**2*delta_eta
      end do

      print*, "errper", errper
      print*, "errnat",errnat
contains
      subroutine field_initialize(epo,N,m_x)
        real(8),dimension(:),intent(inout) :: epo
        integer, intent(in) :: N
        class(cartesian_mesh_1d), intent(in) :: m_x
        real(8) :: x
        integer :: i,j 
   
        do i=1,N
           x=m_x%eta_min+real(i-1,8)*m_x%delta_eta
           epo(i)=sin(x)
        end do         
      end subroutine
end program


