Program test_gyorbit_allmu
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                              pic_para_total2d_base
use paradata_layout, only:    initialize_pic_para_total2d_base, &
                              allocate_parameters_array_2d, &
                              allocate_memory_to_field_2d_ful, &
                              allocate_memory_to_field_2d_gy,  &
                              allocate_memory_to_magfield_2d, &
                              initialize_parameters_2d,  &
                              initialize_parameters_array_2d, &
                              computing_mu_number

use utilities_module, only: f_is_power_of_two, &
                            muarray_euler_maclaurin_choice
use m_mpilayout, only : initialize_layout_with_distributed_2d_array, &
                        get_layout_2d_box_index
use piclayout, only :   root_precompute_data, &
                        parameters_array_2d, &
                        initialize_rootdata_structure, &
                        ful2d_node,ful2dsend_node,gy2d_node,gy2dsend_node
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
                          mpi2d_alltoallv_send_particle_2d,&
                          tp_mpi2d_alltoallv_send_particle_2d, &
                          tp_mpi2d_alltoallv_send_particle_2d_gy, &
                          partition_density_to_grid_gy_allmu

use para_random_sample, only: para_accept_reject_gaussian1d_ful2d_per_per, &
                              para_accprej_gaus2d2v_fulgyro_unifield_per_per, &
                              para_accprej_gaus1d2v_fulgyro_unifield_per_per

use m_moveparticles, only: push_particle_among_box_ful2d_per_per
use m_fieldsolver, only: solve_weight_of_field_among_processes, &
                         solve_field_quasi_neutral, &
                         solve_field_ful, &
                         solve_gyfieldweight_from_fulfield, &
                         compute_equdensity_for_ful, &
                         compute_equdensity_for_gy, &
                         compute_gyrodensity_perturbation

use m_precompute,  only: precompute_ASPL
use m_moveparticles, only: new_position_per_per_ful, &
                           push_particle_among_box_ful2d_per_per, &
                           new_list_ful2d_per_per
use para_write_file, only: open_file,&
                           close_file, &
                           para_write_field_file_2d, &
                           para_write_orbit_file_2d, &
                           para_write_orbit_file_2d_gy
use m_para_orbit, only: borissolve, fulrk4solve, gyrork4solveallmu, &
                        para_obtain_interpolation_elefield_per_per_ful, &
                        sortposition_by_process_ful2d, &
                        compute_f_of_points_out_orbit_ful2d, &
                        fulrkfunc_f_per_per, &
                        compute_f_of_points_in_orbit_ful2d, &
                        gyrork4solveallmu
use orbit_data_base, only: rk4ful2dnode, pointin_node,tp_ful2d_node, &
                           tp_ful2dsend_node, &
                           rk4gy2dnode, tp_gy2d_node,tp_gy2dsend_node
use field_initialize,only: para_initialize_field_2d_mesh
use m_tp_para_orbit, only: tp_push_ful_orbit, &
                           tp_push_gy_orbit_allmu, &
                           tp_ful_solve, &
                           tp_ful2dlist_to_ful2dlist, &
                           tp_sort_particles_among_ranks
use m_precompute, only: precompute_ASPL, &
                        precompute_doublegyroaverage_matrix
use piclayout, only: ful2d_node,gy2d_node,gy2dmu_node

implicit none
include "mpif.h"

    class(pic_para_total2d_base),pointer :: pic2d
    class(root_precompute_data), pointer :: rootdata
    class(parameters_array_2d),  pointer :: pamearray
    int4 :: rank,size,global_sz(2)
    real8 :: delta(2)
    int4  :: num1,num2,row
    real8 :: amp,wave_one,wave_two,amp_eq
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

!!!  PIC part
    real8 :: xmin(2)
    int4  :: index(2),low(2)
    real8 :: trap_coord(4,2),ptc_ratio(4)   

    real8 :: coords(4)   
    int4 :: rank1 
!    real8 :: elef(3)

!!!!  outputdata 
    int4 :: fileitem_boris,fileitem_rk4,fileitem_gy
    character(100) :: filepath_boris,filepath_gy,filepath_rk4,filepath1

!!!! Test para_orbit
    class(rk4ful2dnode), pointer :: partlist,partmp
    class(pointin_node), pointer :: inlist, intmp
    real8 :: f1(6), f2(6),f3(6),f4(6)
    real8 :: elef(3), magf(3)
!   int4 :: row,prank,size,numgr,rank
    int4, dimension(:), pointer :: num,recnum
    real8 :: vec(6),dt
    int4 :: order,h,num_time

!!! Test test_particles
    class(tp_ful2d_node), pointer :: tpful2d_head,tpful2dtmp,tprk4ful2d_head,tprk4ful2dtmp
    class(tp_gy2d_node),  pointer :: tpgy2d_head,tpgy2dtmp
    class(tp_ful2dsend_node), dimension(:), pointer :: tpful2dsend_head,tpful2dsendtmp, &
                             tprk4ful2dsend_head, tprk4ful2dsendtmp
    class(tp_gy2dsend_node), dimension(:), pointer :: tpgy2dsend_head,tpgy2dsendtmp
    int4 :: numleft_rk4,numleft_boris,numleft_gy,numcircle,numgr,numgr_gy
    character(25) :: pushkind
    int4, dimension(:), pointer :: num_rk4,num_gy
    real8 :: rho, theta,x2(2)
    int4 :: rk4order,cell_per_unit(2)
    int4 :: orbit_field=1

!!!! integrated simulation
    class(ful2d_node), pointer :: ful2d_head,ful2dtmp
!   class(gy2d_node),  pointer :: gy2d_head, gy2dtmp
!   class(ful2dsend_node), dimension(:),pointer ::currk
!   class(gy2dsend_node), pointer :: gy2dsend_head,gy2dsendtmp
    class(gy2dmu_node), dimension(:), pointer :: gy2dmu_head, gy2dmutmp 
    int4 :: mu_num,mutest
    real8, dimension(:), pointer :: mu_nodes,mu_weights
    int4, dimension(:), pointer :: munum_partition
    real8 :: integ
    int4 :: sum

!!!! store density
   character(100) :: filepathgy,filepathful
   int4 :: fileitemful,fileitemgy

    allocate(weight(-1:2,-1:2))
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
 
    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/2.0*pi_,2.0*pi_/)
    pic2d%para2d%N_points=50
    pic2d%para2d%iter_number=10
    pic2d%para2d%numcircle=8
    pic2d%para2d%numparticle=300000
    pic2d%para2d%dtgy=1.0
    pic2d%para2d%num_time=15
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
 !   pic2d%para2d%mu=1.0
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/10,10/) 
    pic2d%para2d%dtful=pic2d%para2d%dtgy/real(pic2d%para2d%num_time,8)
    pic2d%para2d%mu_scheme = 1
    !!! particle in cell part
    pic2d%para2d%sigma = 25.0
    pic2d%para2d%tempt = 1.0
    pic2d%para2d%mumin=0.0_f64
    pic2d%para2d%mumax=20._F64
    pic2d%para2d%mu_tail=1000
    pic2d%para2d%mulast = 40
!    pic2d%para2d%mu_num=39
    pic2d%para2d%gyroorder=1
    row=pic2d%para2d%row
    amp=0.02
    amp_eq=0.0
    wave_one=20.0
    wave_two=20.0
    num_time=pic2d%para2d%num_time
    cell_per_unit=pic2d%para2d%cell_per_unit

    !!! initialize layout2d      
    pic2d%layout2d%collective%thread_level_provided=1
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    numproc=pic2d%para2d%numproc 
    comm=pic2d%layout2d%collective%comm
!    mu_num=pic2d%para2d%mu_num


!    allocate(num_p(0:size-1),recnum(0:size-1))
!    allocate(partlist,inlist) 
!
!    allocate(tpful2d_head,tprk4ful2d_head,tpgy2d_head)
!    allocate(tpful2dsend_head(0:size-1),tpful2dsendtmp(0:size-1))
!    allocate(tprk4ful2dsend_head(0:size-1), tprk4ful2dsendtmp(0:size-1))
!    allocate(tpgy2dsend_head(0:size-1), tpgy2dsendtmp(0:size-1))
!    allocate(num_rk4(0:size-1),num_gy(0:size-1))


    
    allocate(mu_nodes(100),mu_weights(100),munum_partition(100))
    munum_partition=0
    call computing_mu_number(mu_nodes,mu_weights,munum_partition,mu_num, pic2d)
    pic2d%para2d%mu_num=mu_num
!   munum_partition=0
    allocate(ful2d_head)
    allocate(gy2dmu_head(1:mu_num),gy2dmutmp(1:mu_num)) 
 
    do i=1, mu_num
      allocate(gy2dmu_head(i)%ptr, gy2dmutmp(i)%ptr)
    end do
!    num_p=0

    do i=1,2
       delta(i)=pic2d%para2d%gxmax(i)/real(cell_per_unit(i)*numproc(i),8)
    end do

    do i=1,2
       global_sz(i)=pic2d%para2d%cell_per_unit(i)*numproc(i)+1
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

!print*, "rank1=",rank
 
    global_sz(1)=pic2d%layout2d%global_sz1
    global_sz(2)=pic2d%layout2d%global_sz2

     pamearray=>allocate_parameters_array_2d(mu_num,global_sz(1:2))

   call initialize_parameters_2d(pic2d,pamearray) 
   call initialize_parameters_array_2d(pic2d%para2d%mumax,mu_num,pic2d%para2d%mu_scheme,pamearray, &
        mu_nodes,mu_weights,munum_partition)  

   sum=0
   do i=1,mu_num
     sum=pamearray%munum_partition(i)+sum
   end do

if(rank==0) then
!print*, "mu_num=",mu_num, "mutest=",mutest
print*, "sum=",sum*size
print*, "munum_partition=",pamearray%munum_partition
print*, "mu_num=",pic2d%para2d%mu_num
print*, "mu_nodes=",pamearray%mu_nodes
print*, "mu_weights=",pamearray%mu_weights
end if

if(rank==0) then
   print*, "gboxmin(:,1)",pic2d%para2d%gboxmin(:,1)

   print*, "gboxmin(:,2)",pic2d%para2d%gboxmin(:,2)
end if

    num1=pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1
    num2=pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1                    

    call allocate_memory_to_field_2d_ful(pic2d%field2d,num1,num2,row)

    call allocate_memory_to_field_2d_gy(pic2d%field2d,num1,num2,row,mu_num)
    call allocate_memory_to_magfield_2D(pic2d%field2d,num1,num2,row)

    rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)

    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

 !!! prepare the initial field
  dimsize=dimsize_of_rank_per_per(rank,pic2d%para2d%numproc,pic2d%layout2d)

  if(orbit_field==0) then
  do i=1,dimsize(1)
    do j=1, dimsize(2)
       globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
       pic2d%field2d%ep(i,j)=real(globalind(2),8)*0.1    !real(globalind(1)+globalind(2), 8)
       pic2d%field2d%Bf03(i,j)=1.0
    end do
  end do
  else
    call para_initialize_field_2d_mesh(amp,amp_eq,wave_one,wave_two, pic2d)
  endif
    
  !!! prepare the inital distribution of particles
  call para_accprej_gaus1d2v_fulgyro_unifield_per_per(ful2d_head,gy2dmu_head,pic2d,pamearray)  

if(rank==0) then
print*, "#samping particles is finished."
end if


 !!! precomputing
  call precompute_ASPL(rank,global_sz,rootdata%ASPL)
if(rank==0) then
print*, "precomputing ASPL is finished."
endif

  call precompute_doublegyroaverage_matrix(rootdata,pic2d,pamearray)
if(rank==0) then
print*, "precomputing is finished."
endif
  call solve_weight_of_field_among_processes(pic2d%field2d%Bf03,rootdata,pic2d, &
       pic2d%field2d%bf03wg,pic2d%field2d%BF03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
       pic2d%field2d%bf03wg_s, pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
       pic2d%field2d%bf03wg_nw,pic2d%field2d%bf03wg_ne)

    call solve_gyfieldweight_from_fulfield(rootdata,pic2d,pamearray)

    call gyrork4solveallmu(gy2dmu_head,pic2d,pic2d%para2d%iter_number)


!  do i=1,pic2d%para2d%iter_number
!    if(rank==0) then
!      print*, "#iter_number=", i
!    endif
!    !!! and store the equilibrium distirbution on the mesh
!    call solve_gyfieldweight_from_fulfield(rootdata,pic2d,pamearray)
!print*, 1
!    call gyrork4solveallmu(gy2dmu_head,pic2d,pic2d%para2d%iter_number)
!print*, 2
!    call partition_density_to_grid_gy_allmu(gy2dmu_head,pic2d)
!!    pic2d%field2d%deng=pic2d%field2d%deng-pic2d%field2d%dengeq
!    call compute_gyrodensity_perturbation(rootdata,pic2d,pamearray)
!print*,3
!    call solve_field_quasi_neutral(rank,rootdata,pic2d,pamearray)  !!! solve the electrostatic potential
!
!!    call diagnoize
!     do j=1,pic2d%para2d%num_time
!       if(rank==0) then
!         print*, "#j=",j
!       endif
!       call solve_weight_of_field_among_processes(pic2d%field2d%ep,rootdata%ASPL,rootdata,pic2d, &
!       pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e,pic2d%field2d%epwg_n,&
!       pic2d%field2d%epwg_s, pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se, &
!       pic2d%field2d%epwg_nw,pic2d%field2d%epwg_ne)
!
!       call borissolve(ful2d_head,pic2d,(i-1)*pic2d%para2d%num_time+j)
!
!       call partition_density_to_grid_ful(ful2d_head,pic2d)
!
!       pic2d%field2d%denf=pic2d%field2d%denf-pic2d%field2d%denfeq
!       call solve_field_ful(rootdata,pic2d,pamearray)
!
!
!!       call diagnoize ????
!    end do
!
!  end do
!print*, rank
    deallocate(pic2d,pamearray,rootdata)
    deallocate(ful2d_head,gy2dmu_head,gy2dmutmp)
 !!! and store the equilibrium distirbution on the mesh
    call MPI_FINALIZE(IERR) 
  end program test_gyorbit_allmu

