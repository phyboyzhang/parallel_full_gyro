program test_alltoallv_per_per
#include "work_precision.h"
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base
use paradata_layout, only: initialize_pic_para_2d_base  
use utilities_module, only: f_is_power_of_two
use m_mpilayout, only : initialize_layout_with_distributed_2d_array
use m_parautilities, only: mpi2d_alltoallv_box_per_per
implicit none
include "mpif.h"

print*, 0
call initialize_pic2d_base()


contains

 subroutine initialize_pic2d_base()      
    class(pic_para_2d_base),pointer :: pic2d
    int4 :: rank,size,global_sz(2)
    real8 :: delta(2)
    int4  :: num1,num2,row
    int4  :: ierr, boxindex(4)
    real8 :: amp=1.0,wave_one=1.0,wave_two=1.0
    character(90) :: geometry="cartesian"
    int4  :: i,j
    print*, 1

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    pic2d => initialize_pic_para_2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/3.0,3.0/)
    pic2d%para2d%N_points=100
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/3,3/) 
    row=pic2d%para2d%row
    do i=1,2
       delta(i)=1._f64/real(pic2d%para2d%cell_per_unit(i),8)
    end do
    print*, 2

    !!! initialize layout2d      
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
!    if(.not.f_is_power_of_two(int(size,8))) then
!      print*, "size is not a numer of power of 2"
!      stop
!    end if
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    do i=1,2
       global_sz(i)=(pic2d%para2d%cell_per_unit(i)+1)*(pic2d%para2d%gxmax(i)-pic2d%para2d%gxmin(i))
    end do
print*, 4
    call initialize_layout_with_distributed_2d_array( &
      global_sz(1), &
      global_sz(2), &
      pic2d%para2d%numproc(1), &
      pic2d%para2d%numproc(2), &
      pic2d%layout2d, &
      pic2d%para2d%boundary)
  
!    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))
    do i=0,size-1
    if(i==0) then
      pic2d%para2d%gboxmin(0,:)=pic2d%para2d%gxmin(:)
      pic2d%para2d%gboxmax(i,1)=pic2d%para2d%gboxmin(i,1)+delta(1)*real(pic2d%layout2d%boxes(i)%i_max &
        -pic2d%layout2d%boxes(i)%i_min,8)
      pic2d%para2d%gboxmax(i,2)=pic2d%para2d%gboxmin(i,2)+delta(2)*real(pic2d%layout2d%boxes(i)%j_max &
        -pic2d%layout2d%boxes(i)%j_min,8) 
    else
      pic2d%para2d%gboxmin(i,:)=pic2d%para2d%gboxmax(i-1,:)+(/delta(1),delta(2)/)
      pic2d%para2d%gboxmax(i,1)=pic2d%para2d%gboxmin(i,1)+delta(1)*real(pic2d%layout2d%boxes(i)%i_max &
        -pic2d%layout2d%boxes(i)%i_min,8)
      pic2d%para2d%gboxmax(i,2)=pic2d%para2d%gboxmin(i,2)+delta(2)*real(pic2d%layout2d%boxes(i)%j_max &
        -pic2d%layout2d%boxes(i)%j_min,8)
    end if
    end do
print*, 4
    pic2d%para2d%m_x1=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1,&
      pic2d%para2d%gboxmin(rank,1),pic2d%para2d%gboxmax(rank,1),delta(1))
    pic2d%para2d%m_x2=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1,&
    pic2d%para2d%gboxmin(rank,2),pic2d%para2d%gboxmax(rank,2),delta(2)) 
    num1=pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1
    num2=pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1                    
    allocate(pic2d%field2d%ep(num1,num2),pic2d%field2d%ep_w(num1,row),&
             pic2d%field2d%ep_e(num1,row),pic2d%field2d%ep_n(row,num2),pic2d%field2d%ep_s(row,num2), &
             pic2d%field2d%ep_sw(row,row),pic2d%field2d%ep_se(row,row),pic2d%field2d%ep_nw(row,row), &
             pic2d%field2d%ep_ne(row,row))
print*, 5
    !!!! initialize the value of epotential
    call compute_field_2d_mesh(pic2d%field2d%ep,pic2d%para2d%m_x1,pic2d%para2d%m_x2,amp,wave_one,wave_two,geometry) 
   
    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max
print*, 7
    call mpi2d_alltoallv_box_per_per(row,mpi_comm_world,rank,pic2d%para2d%numproc, &
                                     pic2d%field2d%ep,pic2d%field2d%ep_w,pic2d%field2d%ep_e,pic2d%field2d%ep_n,&
                                     pic2d%field2d%ep_s,pic2d%field2d%ep_sw,pic2d%field2d%ep_se,pic2d%field2d%ep_nw,&
                                     pic2d%field2d%ep_ne,boxindex)  
print*, 6
    if(rank==0) then
       do i=1,num2
          print*,"i=",i, pic2d%field2d%ep_e(i,:)
       end do
    end if                                  
   
   call MPI_FINALIZE(ierr) 
   end subroutine initialize_pic2d_base



  subroutine compute_field_2d_mesh( &
      field_2d, &
      m_x1,  &
      m_x2,  &
       amp,&
       wave_one, &
       wave_two, &
      geometry)
      real8,dimension(:,:), pointer,intent(inout) :: field_2d
      class(cartesian_mesh_1d), pointer,intent(in) :: m_x1,m_x2
      real8, intent(in) :: amp,wave_one,wave_two
!      real8, intent(in) :: contr
 !     character(len=*), intent(in) :: boundary
      character(len=*), intent(in) :: geometry
      real8 :: x(2),x1(2)
      int4 :: i, j

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
!                wave_two, x)

           field_2d(i+1,j+1) = x(2)

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
               field_2d(i+1,j+1) = epotential_2d_analytical(amp,wave_one, &
                    wave_two, x1)
            end do
         end do

    case default
       print*, "input the correct geometry, geometry=", geometry
       stop
    end select

  end subroutine compute_field_2d_mesh


   function epotential_2d_analytical(waveamplitude, &
                                    wavenumber_one,  &
                                    wavenumber_two,  &
                                    x) result(epo)   ! the fist kind of electric potential. k is the wave vector needing to be scanned
  real8, intent(in) :: x(2)
  real8 :: epo
  real8, intent(in) :: waveamplitude, wavenumber_one, wavenumber_two
!  call wavearray(input)
 !    epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
 !    efi(2)=input(2,2)*sin(input(1,2)*x(1))*cos(input(1,2)*x(2))
!     epo=waveamplitude*sin(wavenumber_one*x(1))*cos(wavenumber_two*x(2))
! epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))
!  epo=waveamplitude*(sin(wavenumber_one*x(1))+cos(wavenumber_two*x(2)))+waveamplitude*(x(2)+x(1))/4.0
       epo=waveamplitude*x(2)
  !    epo=waveamplitude
     return
   end function epotential_2d_analytical

end program test_alltoallv_per_per
