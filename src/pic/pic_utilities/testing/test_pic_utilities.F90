Program test_pic_utilities
#include "work_precision.h"
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                         pic_para_total2d_base
use paradata_layout, only:    initialize_pic_para_2d_base, &
                              initialize_pic_para_total2d_base, &
                              allocate_memory_to_field_2d_ful, &
                              allocate_memory_to_field_2d_gy, &
                              allocate_memory_to_magfield_2d
                             
use utilities_module, only: f_is_power_of_two
use m_mpilayout, only : initialize_layout_with_distributed_2d_array, &
                        get_layout_2d_box_index
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

use spline_module, only: compute_D_spl2D_per_per_noblock
use m_para_spline, only: para_compute_spl2D_weight, &
                         para_compute_spl2d_point_per_per_weight, &
                         para_compute_spl2d_field_point_per_per, &
                         para_spl2d_firstorder_derivatve_point_per_per

use para_gyroaverage_2d_one,only: sort_quadraturepoint_among_process, &
                                  para_compute_gyroaverage_stiff_matrix, &
                                  para_compute_gyroaverage_mesh_field
use gyroaverage_2d_base, only: gyropoint_node

use m_picutilities, only: singlepart_to_mesh, & 
                          grid_global_2dcoord, &
                          particle_to_2dgrid,  &
                          low_localindex_of_particle

implicit none
include "mpif.h"

    class(pic_para_total2d_base),pointer :: pic2d
 !   class(pic_para_total2d_base)
    class(root_precompute_data), pointer :: rootdata
    int4 :: rank,size,global_sz(2)
    real8 :: delta(2)
    int4  :: num1,num2,row
    real8 :: amp=1.0,wave_one=1.0,wave_two=1.0
    character(90) :: geometry="cartesian"
    int4  :: i,j,size1,k
    int4  :: rankone,startind(2),globalind(2)
    int4, dimension(:),pointer :: num_p
    real8 :: x(2),x1(2), fieldvalue
    int4 :: dimsize(2), comm, numproc(2),boxindex(4)
    int4 :: ierr
    real8, dimension(:,:), pointer :: weight
    real8 :: deri_firstorder(2)


!!!!! gyroaverage

    real8 :: mu
    class(gyropoint_node), dimension(:),pointer :: pointhead


!!!  PIC part
    real8 :: xmin(2)
    int4  :: index(2),low(2)
    real8 :: trap_coord(4,2),ptc_ratio(4)   


    allocate(weight(-1:2,-1:2))

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(num_p(0:size-1))
    allocate(pointhead(0:size-1))
    do i=0, size-1
      allocate(pointhead(i)%ptr)
    end do

    num_p=0

    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/3.0,3.0/)
    pic2d%para2d%N_points=3
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/2,2/) 
    row=pic2d%para2d%row

    do i=1,2
       delta(i)=1._f64/real(pic2d%para2d%cell_per_unit(i),8)
    end do

    !!! initialize layout2d      
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    numproc=pic2d%para2d%numproc 
    comm=pic2d%layout2d%collective%comm
 
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
    if(rank==0) then
       do i=0, size-1
         print*, "size,boxes(i),",i,pic2d%layout2d%boxes(i)
       enddo  
    end if 
!    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))
  do i=0,pic2d%para2d%numproc(2)-1      
    do j=0,pic2d%para2d%numproc(1)-1
    size1=i*pic2d%para2d%numproc(1)+j
    if(size1==0) then
       pic2d%para2d%gboxmin(0,:)=pic2d%para2d%gxmin(:)
    else
       startind=startind_of_process(size1,pic2d%para2d%numproc,pic2d%layout2d) 
       if(rank==0) then
         print*, "size1=",size1,"startind=",startind
       end if
       do k=1,2 
          pic2d%para2d%gboxmin(size1,k)=pic2d%para2d%gxmin(k)+(startind(k)-1)*delta(k)
       end do
    end if

    pic2d%para2d%gboxmax(size1,1)=pic2d%para2d%gboxmin(size1,1)+delta(1)*real(pic2d%layout2d%boxes(size1)%i_max &
    -pic2d%layout2d%boxes(size1)%i_min,8)
    pic2d%para2d%gboxmax(size1,2)=pic2d%para2d%gboxmin(size1,2)+delta(2)*real(pic2d%layout2d%boxes(size1)%j_max &
    -pic2d%layout2d%boxes(size1)%j_min,8)

    end do
  end do

if(rank==0) then
   print*, "gboxmin(:,1)",pic2d%para2d%gboxmin(:,1)

   print*, "gboxmin(:,2)",pic2d%para2d%gboxmin(:,2)
end if

    pic2d%para2d%m_x1=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1,&
  pic2d%para2d%gboxmin(rank,1),pic2d%para2d%gboxmax(rank,1),delta(1))
    pic2d%para2d%m_x2=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1,&
    pic2d%para2d%gboxmin(rank,2),pic2d%para2d%gboxmax(rank,2),delta(2)) 
    num1=pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1
    num2=pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1                    

call allocate_memory_to_field_2d_ful(pic2d%field2d,num1,num2,row)
call allocate_memory_to_field_2d_gy(pic2d%field2d,num1,num2,row,1)
call allocate_memory_to_magfield_2D(pic2d%field2d,num1,num2,row)
    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

  !!!periodic boundary condition
 
  dimsize=dimsize_of_rank_per_per(rank,pic2d%para2d%numproc,pic2d%layout2d)
  do i=1,dimsize(1)
    do j=1, dimsize(2)
       globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
       pic2d%field2d%ep(i,j)=1.0    ! real(globalind(1)+globalind(2), 8)
    end do
  end do
 
  rankone=3
  if(rank==rankone) then
    x(1)=pic2d%para2d%gboxmin(rankone,1)+0.3
    x(2)=pic2d%para2d%gboxmin(rankone,2)+0.4
    xmin(:)=pic2d%para2d%gboxmin(rank,:)
    delta(:)=(/pic2d%para2d%m_x1%delta_eta,pic2d%para2d%m_x2%delta_eta/)
    call get_layout_2d_box_index( pic2d%layout2d, rank, boxindex )    
    index(1)=boxindex(1) 
    index(2)=boxindex(3)
 
    call grid_global_2dcoord(trap_coord,x,delta,xmin)
    do i=1,4
    print*, trap_coord(i,:)
    end do
    call particle_to_2dgrid(ptc_ratio, trap_coord, x)
 !   print*, "ptc_ratio=",ptc_ratio
    call low_localindex_of_particle(low,x,xmin,delta)
 !   print*, "low=",low

    call singlepart_to_mesh(pic2d%field2d%denf,x,delta,xmin) 
    print*, pic2d%field2d%denf 


  end if




    call MPI_FINALIZE(IERR) 
  end program test_pic_utilities
