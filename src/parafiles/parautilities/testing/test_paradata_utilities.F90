Program test_paradata_utilities
#include "work_precision.h"
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base
use paradata_layout, only: initialize_pic_para_2d_base
                              
use utilities_module, only: f_is_power_of_two
use m_mpilayout, only : initialize_layout_with_distributed_2d_array
use piclayout, only :   root_precompute_data, &
                        initialize_rootdata_structure
use m_parautilities, only: mpi2d_alltoallv_box_per_per, &
                           gather_field_to_rootprocess_per_per, &
                           scatter_field_from_rootprocess_per_per

use paradata_utilities, only: compute_rank_from_globalind_2d, &
                              compute_process_of_point_per_per, &
                              startind_of_process, &
                              globalind_from_localind_2d, &
                              dimsize_of_rank_per_per, &
                              copy_boundary_value_per_per
implicit none
include "mpif.h"

    class(pic_para_2d_base),pointer :: pic2d
    class(root_precompute_data), pointer :: rootdata
    int4 :: rank,size,global_sz(2)
    real8 :: delta(2)
    int4  :: num1,num2,row
    int4  :: ierr, boxindex(4)
    real8 :: amp=1.0,wave_one=1.0,wave_two=1.0
    character(90) :: geometry="cartesian"
    int4  :: i,j,size1
    int4  :: rankone,startind(2),globalind(2)
    int4, dimension(:),pointer :: num_p
    real8 :: x(2)
    int4 :: dimsize(2)

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(num_p(0:size-1))
    num_p=0
    pic2d => initialize_pic_para_2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/2.0,2.0/)
    pic2d%para2d%N_points=100
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=2
    pic2d%para2d%cell_per_unit=(/3,3/) 
    row=pic2d%para2d%row
    do i=1,2
       delta(i)=1._f64/real(pic2d%para2d%cell_per_unit(i),8)
    end do

    !!! initialize layout2d      
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    if(.not.f_is_power_of_two(int(size,8))) then
      print*, "size is not a numer of power of 2"
      stop
    end if
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    do i=1,2
       global_sz(i)=pic2d%para2d%cell_per_unit(i)*(pic2d%para2d%gxmax(i)-pic2d%para2d%gxmin(i))+1
    end do

    call initialize_layout_with_distributed_2d_array( &
      global_sz(1), &
      global_sz(2), &
      pic2d%para2d%numproc(1), &
      pic2d%para2d%numproc(2), &
      pic2d%layout2d, &
      pic2d%para2d%boundary)
!    if(rank==0) then
!       do i=0, size-1
!         print*, "boxes(i),i=",i,pic2d%layout2d%boxes(i)
!       enddo  
!    end if 
!    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))
    do i=0,pic2d%para2d%numproc(2)-1      
    do j=0,pic2d%para2d%numproc(1)-1
    size1=i*pic2d%para2d%numproc(1)+j
    if(size1==0) then
    pic2d%para2d%gboxmin(0,:)=pic2d%para2d%gxmin(:)
    else
    startind=startind_of_process(size1,pic2d%para2d%numproc,pic2d%layout2d)  
    pic2d%para2d%gboxmin(size1,:)=pic2d%para2d%gxmin(:)+(startind(:)-(/1,1/))*delta(:)
    end if

    pic2d%para2d%gboxmax(size1,1)=pic2d%para2d%gboxmin(size1,1)+delta(1)*real(pic2d%layout2d%boxes(size1)%i_max &
    -pic2d%layout2d%boxes(size1)%i_min,8)
    pic2d%para2d%gboxmax(size1,2)=pic2d%para2d%gboxmin(size1,2)+delta(2)*real(pic2d%layout2d%boxes(size1)%j_max &
    -pic2d%layout2d%boxes(size1)%j_min,8)

    end do
  end do

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
    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

  !!!periodic boundary condition
 
  dimsize=dimsize_of_rank_per_per(rank,pic2d%para2d%numproc,pic2d%layout2d)
  do i=1,dimsize(1)
    do j=1, dimsize(2)
       globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
       pic2d%field2d%ep(i,j)=real(globalind(1)+globalind(2), 8)
    end do
  end do
     if(rank==2) then
        print*, "ep=", pic2d%field2d%ep
    end if

call copy_boundary_value_per_per(pic2d%field2d%ep,rank,pic2d%para2d%numproc,pic2d%layout2d)

     if(rank==2) then
        print*, "ep=", pic2d%field2d%ep
    end if

    rankone=compute_rank_from_globalind_2d(5,3,pic2d%para2d%numproc,pic2d%layout2d)
!    print*, "rank=",rank,"rankone=",rankone
  
!    print*, "rank=",rank, boxindex(:)
!    if(rank==0) then
!      do i=0,size-1
!          print*,"rank=",i,pic2d%para2d%gboxmin(i,:),pic2d%para2d%gboxmax(i,:)
!      end do
!    end if

     x=(/1.2,1.3/)      
    rankone=compute_process_of_point_per_per(x,pic2d%para2d%numproc,pic2d%para2d%gxmin,pic2d%para2d%gxmax, &
         pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)


    rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)
    call gather_field_to_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,rank,size,boxindex,&
          pic2d%para2d%numproc,pic2d%layout2d) 

    if(rank==0) then
       print*, 'rootfield=', rootdata%field
    end if

    pic2d%field2d%ep=0._f64
    call scatter_field_from_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,size, &
         pic2d%para2d%numproc,(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/),pic2d%layout2d)
     
    if(rank==0) then
       print*, 'ep=', pic2d%field2d%ep
    endif
 
    call MPI_FINALIZE(IERR) 
end program test_paradata_utilities
