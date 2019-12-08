Program test_random_sample
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                              pic_para_total2d_base
use paradata_layout, only:    initialize_pic_para_total2d_base, &
                              allocate_parameters_array_2d

use m_mpilayout, only : initialize_layout_with_distributed_2d_array, &
                        get_layout_2d_box_index
use piclayout, only :   root_precompute_data, &
                        parameters_array_2d, &
                        initialize_rootdata_structure, &
                        ful2d_node,ful2dsend_node,gy2d_node, &
                        gy2dsend_node, gy2dmu_node
implicit none
include "mpif.h"

    class(pic_para_total2d_base),pointer :: pic2d
    class(parameters_array_2d), pointer :: pamearray
    int4 :: rank,size,global_sz(2)
    int4  :: num1,num2,row
    character(90) :: geometry="cartesian"
    int4  :: i,j,size1,k
!    int4, dimension(:),pointer :: num_p
    int4 :: dimsize(2), comm, numproc(2),boxindex(4)
    int4 :: ierr, cell_per_unit(2)
    real8 :: mu,mumax

    int4 :: mu_num


    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

      
!    pic2d=> initialize_pic_para_total2d_base(size)
!    pic2d%para2d%gxmin=(/0.0,0.0/)
!    pic2d%para2d%gxmax=(/2.0*pi_,2.0*pi_/)
!    pic2d%para2d%mu=1.0
!    pic2d%para2d%row=3
!    pic2d%para2d%cell_per_unit=(/10,10/) 
!    pic2d%para2d%mumin=0.0_f64
!    pic2d%para2d%mumax=20._F64
!    pic2d%para2d%mu_num=10
!    pic2d%para2d%boundary = "double_per"
!    cell_per_unit=pic2d%para2d%cell_per_unit
!
!    !!! initialize layout2d      
!    pic2d%layout2d%collective%rank=rank
!    pic2d%layout2d%collective%size=size
!    pic2d%layout2d%collective%comm=mpi_comm_world
!    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
!    numproc=pic2d%para2d%numproc 
!    comm=pic2d%layout2d%collective%comm
!    mu_num=pic2d%para2d%mu_num
!
!    do i=1,2
!       global_sz(i)=pic2d%para2d%cell_per_unit(i)*numproc(i)+1
!    end do

!    call initialize_layout_with_distributed_2d_array( &
!      global_sz(1), &
!      global_sz(2), &
!      pic2d%para2d%numproc(1), &
!      pic2d%para2d%numproc(2), &
!      pic2d%layout2d, &
!      pic2d%para2d%boundary)
!    if(rank==0) then
!       do i=0, size-1
!         print*, "size,boxes(i),",i,pic2d%layout2d%boxes(i)
!       enddo  
!    end if
     mu_num=10
     global_sz(1)= 2
     global_sz(2)= 2
!    global_sz(1)=pic2d%layout2d%global_sz1
!    global_sz(2)=pic2d%layout2d%global_sz2
!print*, mu_num,global_sz
    pamearray=> allocate_parameters_array_2d(mu_num,global_sz)
deallocate(pamearray)
call mpi_barrier(mpi_comm_world) 
print*, rank


!    deallocate(pic2d,pamearray,rootdata)
     call MPI_FINALIZE(IERR) 
  end program  test_random_sample

