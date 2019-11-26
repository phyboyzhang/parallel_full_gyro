Program test_para_partgrid_coupl_one
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
                        initialize_rootdata_structure, &
                        ful2d_node, &
                        ful2dsend_node
use m_parautilities, only: mpi2d_alltoallv_box_per_per, &
                           gather_field_to_rootprocess_per_per, &
                           scatter_field_from_rootprocess_per_per

use paradata_utilities, only: compute_rank_from_globalind_2d, &
                              compute_process_of_point_per_per, &
                              startind_of_process, &
                              globalind_from_localind_2d, &
                              dimsize_of_rank_per_per, &
                              copy_boundary_value_per_per, &
                              coordinates_pointoutbound_per_per
                            

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
                          low_localindex_of_particle, &
                          partition_density_to_grid_ful, &
                          mpi2d_alltoallv_send_particle_2d

use para_random_sample, only: para_accept_reject_gaussian1d_ful2d_per_per, &
                              para_accprej_gaus2d2v_fulgyro_unifield_per_per, &
                              para_accprej_gaus1d2v_fulgyro_unifield_per_per

use m_fieldsolver, only: solve_weight_of_field_among_processes
use m_precompute,  only: precompute_ASPL
use m_moveparticles, only: new_position_per_per_ful, &
                           push_particle_among_box_ful2d_per_per, &
                           new_list_ful2d_per_per
use para_write_file, only: open_file,&
                           close_file, &
                           para_write_field_file_2d
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
    real8 :: x(2),v(2),x1(2), fieldvalue
    int4 :: dimsize(2), comm, numproc(2),boxindex(4)
    int4 :: ierr
    real8, dimension(:,:), pointer :: weight
    real8 :: deri_firstorder(2)


!!!!! gyroaverage

    real8 :: mu,mumax
    class(gyropoint_node), dimension(:),pointer :: pointhead


!!!  PIC part
    real8 :: xmin(2)
    int4  :: index(2),low(2)
    real8 :: trap_coord(4,2),ptc_ratio(4)   

    class(ful2d_node),pointer :: curful
    class(ful2dsend_node), dimension(:),pointer ::currk
    real8 :: coords(4)   
    int4 :: rank1 

!!!!  outputdata 
    int4 :: fileitem1
    character(100) :: filepath1

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
    pic2d%para2d%gxmax=(/5.0,5.0/)
    pic2d%para2d%N_points=3
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/10,10/) 

    !!! particle in cell part
    pic2d%para2d%numparticle=550000
    pic2d%para2d%numcircle = 10
    pic2d%para2d%sigma = 1.0
    pic2d%para2d%tempt = 1.0
    pic2d%para2d%mumin=0.0_f64
    pic2d%para2d%mumax=20._F64
    row=pic2d%para2d%row

    do i=1,2
       delta(i)=1._f64/real(pic2d%para2d%cell_per_unit(i),8)
    end do

    !!! initialize layout2d      
    pic2d%layout2d%collective%thread_level_provided=1
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

    global_sz(1)=pic2d%layout2d%global_sz1
    global_sz(2)=pic2d%layout2d%global_sz2 
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
   rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)

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
num_p=0
   curful=>pic2d%ful2d_head 
   fileitem1=10
   filepath1="/PARA/blsc950/electrostatic_exp/run/density.txt"
!  call para_accept_reject_gaussian1d_ful2d_per_per(pic2d)
!  call para_accprej_gaus1d2v_fulgyro_unifield_per_per(curful,pic2d)  
!  call sort_particle_among_box_ful2d_per_per(pic2d,num_p)
  call partition_density_to_grid_ful(curful,pic2d) 
!
   call open_file(fileitem1,filepath1,rank)  
   
   call para_write_field_file_2d(pic2d%field2d%denf,fileitem1,rootdata,pic2d)

   call close_file(fileitem1,rank)

    call MPI_FINALIZE(IERR) 
  end program  test_para_partgrid_coupl_one
