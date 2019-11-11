program test_interpotwo
#include "work_precision.h"
  use  spline_module, only : &
       s_splcoefper1d0old, &
       s_compute_splines_coefs_matrix_nat_1d_new, &
       s_compute_splines_coefs_matrix_per_1d


implicit none

real8 :: r_min=1.0_f64, r_max=3.0_f64
int4  :: Nr=4, Ntheta=4, &
     N_points=4,gyroorder=2
real8,allocatable :: rho(:)
int4, allocatable  :: pre_compute_N(:)
real8, allocatable :: points(:,:)
int4  :: size_pre_compute
int4  :: i,j,m,ierr
real8  :: mu=0.01_f64
comp8, dimension(:,:,:), allocatable :: mat,mat_real
real8, dimension(:), allocatable :: buf_fft
comp8, dimension(:), allocatable :: fft_array
real8, dimension(:), allocatable :: unit_1,unit_2
int4, dimension(:,:), allocatable :: pre_compute_index
real8, dimension(:), allocatable  :: pre_compute_coeff_spl
real8,dimension(:,:,:),allocatable,target::mat_spl2D_circ
real8,dimension(:,:,:),pointer::pointer_mat_spl2D_circ


allocate(rho(Nr+1))
allocate(pre_compute_N(Nr+1))
allocate(points(3,N_points))
allocate(mat(Ntheta,Nr+1,Nr+1))
allocate(mat_real(Ntheta,Nr+1,Nr+1))
allocate(buf_fft(4*Ntheta+15))
allocate(fft_array(Ntheta))
allocate(unit_1(Ntheta*(Nr+1)))
allocate(unit_2(Nr+1))
ALLOCATE(mat_spl2D_circ(0:Ntheta-1,0:Nr,0:Nr))

call  compute_D_spl2D
    num_cells_r, &
    Nr, &
    Ntheta, &
    D_spl2D)

end program test_interpotwo

  subroutine compute_D_spl2D
    num_cells_r, &
    num_cells_theta, &
    D_spl2D)
    int4, intent(in) :: num_cells_r
    int4, intent(in) :: num_cells_theta
    int4 :: ierr
    int4 :: Nr
    int4 :: Ntheta
    comp8, dimension(:,:,:), intent(out) :: D_spl2D
    real8,,dimension(:),allocatable,target::dnat,lnat,dper,lper,mper
    real8,dimension(:),pointer::pointer_dnat,pointer_lnat
    real8,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    real8,dimension(:,:),allocatable,target::mat_nat,mat_per
    real8,dimension(:,:),pointer::pointer_mat_nat,pointer_mat_per
    real8,dimension(:,:,:),allocatable,target::mat_spl2D_circ
   pointer_mat_per => mat_per
    pointer_mat_spl2D_circ => mat_spl2D_circ

    call s_splcoefper1d0old(dper,lper,mper,Ntheta)

    call s_compute_splines_coefs_matrix_nat_1d_new(pointer_mat_nat,Nr+1)
!    print*, 'cont2'
    call s_compute_splines_coefs_matrix_per_1d(pointer_mat_per,pointer_dper,pointer_lper,pointer_mper,Ntheta)
!    print*, 'cont3'
    do j=0,Ntheta-1
      pointer_mat_spl2D_circ(j,:,:)=pointer_mat_per(0,j)*pointer_mat_nat
   enddo

    do k=0,Nr
      do j=0,Nr
        D_spl2D(1:Ntheta,j+1,k+1) = pointer_mat_spl2D_circ(0:Ntheta-1,j,k)*(1._f64,0._f64)

    enddo
 enddo




  end subroutine compute_D_spl2D
