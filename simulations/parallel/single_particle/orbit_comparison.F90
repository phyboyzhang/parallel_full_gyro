Program orbit_comparison
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                              pic_para_total2d_base
use paradata_layout, only:    initialize_pic_para_total2d_base, &
                              initialize_parameters_array_2d, &
                              allocate_parameters_array_2d, &
                              allocate_memory_to_field_2d_ful, &
                              allocate_memory_to_field_2d_gy,  &
                              allocate_memory_to_magfield_2d
                              
use utilities_module, only: f_is_power_of_two
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
                          tp_mpi2d_alltoallv_send_particle_2d_gy

use para_random_sample, only: para_accept_reject_gaussian1d_ful2d_per_per, &
                              para_accprej_gaus2d2v_fulgyro_unifield_per_per, &
                              para_accprej_gaus1d2v_fulgyro_unifield_per_per

use m_moveparticles, only: push_particle_among_box_ful2d_per_per
use m_fieldsolver, only: solve_weight_of_field_among_processes, &
                         solve_gyfieldweight_from_fulfield
use m_precompute,  only: precompute_ASPL
use m_moveparticles, only: new_position_per_per_ful, &
                           push_particle_among_box_ful2d_per_per, &
                           new_list_ful2d_per_per
use para_write_file, only: open_file,&
                           close_file, &
                           para_write_field_file_2d, &
                           para_write_orbit_file_2d, &
                           para_write_orbit_file_2d_gy, &
                           para_write_orbit_file_2d_gy_allmu
use m_para_orbit, only: para_obtain_interpolation_elefield_per_per_ful, &
                        sortposition_by_process_ful2d, &
                        compute_f_of_points_out_orbit_ful2d, &
                        fulrkfunc_f_per_per, &
                        compute_f_of_points_in_orbit_ful2d
use orbit_data_base, only: rk4ful2dnode, pointin_node,tp_ful2d_node, &
                           tp_ful2dsend_node, &
                           rk4gy2dnode, tp_gy2d_node,tp_gy2dsend_node, &
                           tp_gy2dallmu_node
use field_initialize,only: para_initialize_field_2d_mesh
use m_tp_para_orbit, only: tp_push_ful_orbit, &
                           tp_push_gy_orbit_allmu, &
                           tp_ful_solve, &
                           tp_ful2dlist_to_ful2dlist, &
                           tp_sort_particles_among_ranks
implicit none
include "mpif.h"

    class(pic_para_total2d_base),pointer :: pic2d
 !   class(pic_para_total2d_base)
    class(root_precompute_data), pointer :: rootdata 
    class(parameters_array_2d),  pointer :: pamearray
    int4 :: rank,size,global_sz(2)
    real8 :: delta(2)
    int4  :: num1,num2,row
    real8 :: amp,wave_one,wave_two,amp_eq
    character(90) :: geometry="cartesian"
    int4  :: i,j,size1,k,h
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

    class(ful2d_node),pointer :: curful,ful2d_head
    class(ful2dsend_node), dimension(:),pointer ::currk
    class(gy2d_node), pointer :: gy2d_head, gy2dtmp
    class(gy2dsend_node), pointer :: gy2dsend_head,gy2dsendtmp
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
   int4 :: order,num_time

!!! Test test_particles
   class(tp_ful2d_node), pointer :: tpful2d_head,tpful2dtmp,tprk4ful2d_head,tprk4ful2dtmp
   class(tp_gy2dallmu_node),dimension(:), pointer :: tpgy2dmu_head,tpgy2dmutmp
   class(tp_ful2dsend_node), dimension(:), pointer :: tpful2dsend_head,tpful2dsendtmp, &
                             tprk4ful2dsend_head, tprk4ful2dsendtmp
   class(tp_gy2dsend_node), dimension(:), pointer :: tpgy2dsend_head,tpgy2dsendtmp
   int4 :: numleft_rk4,numleft_boris,numcircle,numgr,numgr_gy
   int4, dimension(:), pointer :: numleft_gy
   character(25) :: pushkind
   int4, dimension(:), pointer :: num_rk4
   int4, dimension(:), pointer :: num_gy 
   real8 :: rho, theta,x2(2)
   int4 :: rk4order,cell_per_unit(2)
   int4 :: mu_num
   int4 :: orbit_field=1
   character(250),dimension(:),pointer :: mugyfileitems
   character(100) :: muth  

    allocate(weight(-1:2,-1:2))

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/2.0*pi_,2.0*pi_/)
    pic2d%para2d%N_points=50
    pic2d%para2d%iter_number=100
    pic2d%para2d%numcircle=3
    pic2d%para2d%dtgy=1.0
    pic2d%para2d%num_time=15
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=1.0
    pic2d%para2d%mu_num=1.0
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/30,30/) 
    pic2d%para2d%dtful=pic2d%para2d%dtgy/real(pic2d%para2d%num_time,8)
    !!! particle in cell part
    pic2d%para2d%sigma = 1.0
    pic2d%para2d%tempt = 1.0
    pic2d%para2d%mumin=0.0_f64
    pic2d%para2d%mumax=20._F64
    pic2d%para2d%gyroorder=1
    row=pic2d%para2d%row
    amp=0.001  !0.02
    amp_eq=0.0
    wave_one=1.0
    wave_two=1.0
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
    mu_num=pic2d%para2d%mu_num

    allocate(num_p(0:size-1),recnum(0:size-1))
    allocate(partlist,inlist) 

    allocate(tpful2d_head,tprk4ful2d_head)
    allocate(tpful2dsend_head(0:size-1),tpful2dsendtmp(0:size-1))
    allocate(tprk4ful2dsend_head(0:size-1), tprk4ful2dsendtmp(0:size-1))
    allocate(tpgy2dsend_head(0:size-1), tpgy2dsendtmp(0:size-1))
    allocate(num_rk4(0:size-1),num_gy(0:size-1))
    allocate(ful2d_head)

    allocate(tpgy2dmu_head(1),tpgy2dmutmp(1))
    allocate(tpgy2dmu_head(1)%ptr) 
    allocate(numleft_gy(mu_num))
    allocate(mugyfileitems(mu_num))

    allocate(pointhead(0:size-1)) 
    do i=0, size-1
      allocate(pointhead(i)%ptr)
    end do

    num_p=0
    do i=1,2
       delta(i)=pic2d%para2d%gxmax(i)/real(cell_per_unit(i)*numproc(i),8)
    end do

    do i=1,2
       global_sz(i)=pic2d%para2d%cell_per_unit(i)*numproc(i)+1
    end do 


    call initialize_layout_with_distributed_2d_array( &
      global_sz(1), &
      global_sz(2), &
      numproc(1), &
      numproc(2), &
      pic2d%layout2d, &
      pic2d%para2d%boundary)
    if(rank==0) then
       do i=0, size-1
         print*, "size,boxes(i),",i,pic2d%layout2d%boxes(i)
       enddo  
    end if
    
    global_sz(1)=pic2d%layout2d%global_sz1
    global_sz(2)=pic2d%layout2d%global_sz2 
    pamearray=>allocate_parameters_array_2d(mu_num,global_sz)

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
call allocate_memory_to_field_2d_gy(pic2d%field2d,num1,num2,row,mu_num)
call allocate_memory_to_magfield_2D(pic2d%field2d,num1,num2,row)

   rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)

    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

  !!!periodic boundary condition
!   call para_initialize_field_2d_mesh(amp,wave_one,wave_two, pic2d)
  dimsize=dimsize_of_rank_per_per(rank,pic2d%para2d%numproc,pic2d%layout2d)

  if(orbit_field==0) then
  do i=1,dimsize(1)
    do j=1, dimsize(2)
       globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
       pic2d%field2d%gep(i,j)=real(globalind(2),8)*0.1    ! real(globalind(1)+globalind(2), 8)
       pic2d%field2d%ep(i,j)=real(globalind(2),8)*0.1    ! real(globalind(1)+globalind(2), 8) 
       pic2d%field2d%Bf03(i,j)=1.0
    end do
  end do
  else
    call para_initialize_field_2d_mesh(amp,amp_eq,wave_one,wave_two, pic2d) 
  endif

  !  call para_compute_gyroaverage_mesh_field(pic2d%para2d%mu,1,pic2d)

  if(rank==0) then
    call compute_D_spl2D_per_per_noblock( &
         pic2d%layout2d%global_sz1, &
         pic2d%layout2d%global_sz2, &
         rootdata%ASPL)
   end if
  call solve_weight_of_field_among_processes(pic2d%field2d%ep,rootdata,pic2d, &
       pic2d%field2d%ep_weight, pic2d%field2d%epwg_w,pic2d%field2d%epwg_e,pic2d%field2d%epwg_n, &
       pic2d%field2d%epwg_s, pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se, &
       pic2d%field2d%epwg_nw,pic2d%field2d%epwg_ne)
 
  call solve_weight_of_field_among_processes(pic2d%field2d%gep,rootdata,pic2d, &
       pic2d%field2d%gep_weight, pic2d%field2d%gepwg_w,pic2d%field2d%gepwg_e,pic2d%field2d%gepwg_n, &
       pic2d%field2d%gepwg_s, pic2d%field2d%gepwg_sw,pic2d%field2d%gepwg_se, &
       pic2d%field2d%gepwg_nw,pic2d%field2d%gepwg_ne)
  
  call para_compute_gyroaverage_mesh_field(pic2d%para2d%mu,1,pic2d)
!print*, rank
  
  call solve_gyfieldweight_from_fulfield(rootdata,pic2d,pamearray)

  call solve_weight_of_field_among_processes(pic2d%field2d%Bf03,rootdata,pic2d, &
       pic2d%field2d%bf03wg, pic2d%field2d%BF03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
       pic2d%field2d%bf03wg_s, pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
       pic2d%field2d%bf03wg_nw,pic2d%field2d%bf03wg_ne)

!if(rank==3)  then
!print*, pic2d%field2d%ep_weight
!print*, pic2d%field2d%bf03wg
!end if

   do i=0, size-1
     allocate(tpful2dsend_head(i)%ptr)
     allocate(tprk4ful2dsend_head(i)%ptr)
     allocate(tpgy2dsend_head(i)%ptr)
     tpful2dsendtmp(i)%ptr=>tpful2dsend_head(i)%ptr
     tprk4ful2dsendtmp(i)%ptr=>tprk4ful2dsend_head(i)%ptr
     tpgy2dsendtmp(i)%ptr=>tpgy2dsend_head(i)%ptr
   end do

   numcircle=pic2d%para2d%numcircle 
   tpful2dtmp=>tpful2d_head
   tprk4ful2dtmp=>tprk4ful2d_head
!   do i=1,mu_num
      tpgy2dmutmp(1)%ptr=>tpgy2dmu_head(1)%ptr
!   end do   

   numgr=5
   numgr_gy=4
   num_p=0
   num_rk4=0
   num_gy=0
 
  if(rank==0) then
   x1(1)=(pic2d%para2d%gxmin(1)+pic2d%para2d%gxmax(1))/2.0+0.5
   x1(2)=(pic2d%para2d%gxmin(2)+pic2d%para2d%gxmax(2))/2.0
   rho=sqrt(2._f64*pic2d%para2d%mu)

   rank1=compute_process_of_point_per_per(x1(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
         pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)   
!   print*, rank1,"x1=",x1
       if(rank1==rank) then
         tpgy2dmutmp(1)%ptr%coords(1:2)=x1(1:2)
         tpgy2dmutmp(1)%ptr%coords(3)  =pic2d%para2d%mu
         tpgy2dmutmp(1)%ptr%tp=1
         num_gy(rank1)=num_gy(rank1)+1
         allocate(tpgy2dmutmp(1)%ptr%next)
         tpgy2dmutmp(1)%ptr=>tpgy2dmutmp(1)%ptr%next
       else
         tpgy2dsendtmp(rank1)%ptr%coords(1:2)=x1(1:2)
         tpgy2dsendtmp(rank1)%ptr%coords(3)  =pic2d%para2d%mu
         tpgy2dsendtmp(rank1)%ptr%tp=1
         num_gy(rank1)=num_gy(rank1)+1
         allocate(tpgy2dsendtmp(rank1)%ptr%next)
         tpgy2dsendtmp(rank1)%ptr=>tpgy2dsendtmp(rank1)%ptr%next
       endif

   do j=0,numcircle-1
       theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)
       coords(1)=x1(1)+rho*cos(theta)
       coords(2)=x1(2)+rho*sin(theta)
       coords(3)=-rho*cos(theta+pi_/2.0_f64)
       coords(4)=-rho*sin(theta+pi_/2.0_f64)
!print*, "coords=",coords(1:4)
!       call coordinates_pointoutbound_per_per(coords(1:2),pic2d)
       rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
             pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
       if(rank1==rank) then
         tpful2dtmp%coords(1:4)=coords(1:4)
         tpful2dtmp%tp=j+1
         num_p(rank1)=num_p(rank1)+1
         allocate(tpful2dtmp%next)
         tpful2dtmp=>tpful2dtmp%next

         tprk4ful2dtmp%coords(1:4)=coords(1:4)
         tprk4ful2dtmp%tp=j+1
         num_rk4(rank1)=num_rk4(rank1)+1
         allocate(tprk4ful2dtmp%next)
         tprk4ful2dtmp=>tprk4ful2dtmp%next
 
       else
         tpful2dsendtmp(rank1)%ptr%coords(1:4)=coords(1:4)
         tpful2dsendtmp(rank1)%ptr%tp=j+1
         num_p(rank1)=num_p(rank1)+1
         allocate(tpful2dsendtmp(rank1)%ptr%next)
         tpful2dsendtmp(rank1)%ptr=>tpful2dsendtmp(rank1)%ptr%next

         tprk4ful2dsendtmp(rank1)%ptr%coords(1:4)=coords(1:4)
         tprk4ful2dsendtmp(rank1)%ptr%tp=j+1
         num_rk4(rank1)=num_rk4(rank1)+1
         allocate(tprk4ful2dsendtmp(rank1)%ptr%next)
         tprk4ful2dsendtmp(rank1)%ptr=>tprk4ful2dsendtmp(rank1)%ptr%next
 
       endif
       
     end do       
   end if
    call tp_mpi2d_alltoallv_send_particle_2d(tpful2d_head,numleft_boris,tpful2dsend_head,num_p,numgr,pic2d)
   call tp_mpi2d_alltoallv_send_particle_2d(tprk4ful2d_head,numleft_rk4,tprk4ful2dsend_head,num_rk4,numgr,pic2d)
       call tp_mpi2d_alltoallv_send_particle_2d_gy(tpgy2dmu_head,numleft_gy(1),tpgy2dsend_head, &
                                                num_gy,numgr_gy,1,pic2d)


!tpgy2dtmp=>tpgy2d_head
!do while(associated(tpgy2dtmp)) 
!   if(.not.associated(tpgy2dtmp%next)) then
!      exit
!   else
!      print*, "rank=",rank,"coords=",tpgy2dtmp%coords(1:3)
!      tpgy2dtmp=>tpgy2dtmp%next
!   end if
!end do
!print*, numleft_gy,"rank=",rank,"num_gy=",num_gy

!print*, "rank=",rank,"num_P", num_p
!         call tp_push_ful_orbit(tpful2d_head,numleft,numgr,pic2d,pushkind) 

! print*, "rank=",rank,"num_p",num_p   
!    deallocate(tpful2dsend_head)
!    allocate(tpful2dsend_head(0:size-1))
!    do i=0,size-1
!      allocate(tpful2dsend_head(i)%ptr)
!    end do
!    num_p=0
!    if(numleft==0) then
!      goto 100
!    else
!      print*, "rank=",rank
!      call tp_ful_solve(tpful2d_head,pic2d,pushkind)
!      print*, "rank1=",rank
!    end if

!call tp_sort_particles_among_ranks(tpful2d_head,tpful2dsend_head,pic2d,num_p)
!print*, "rank1=",rank, "num=",num_p 


    pushkind="rk4"
    if(rank==0) then 
      print*, "#pushkind=", pushkind
    end if
   fileitem_boris=100
   fileitem_rk4  =200
   filepath1="/home/qmlu/zsx163/parallel_full_gyro/run/orbit/orbit_"

   filepath_boris=trim(filepath1)//"boris"//".txt"
   filepath_rk4  =trim(filepath1)//"rk4"//".txt"
     call open_file(fileitem_boris,filepath_boris,rank)
     call open_file(fileitem_rk4,  filepath_rk4,  rank)
     i=1
       write(unit=muth,fmt=*) i
       filepath_gy=trim(filepath1)//"gy_"//trim(adjustl(muth))//".txt"        
       call open_file(i,   filepath_gy,   rank)
    
     
!tpful2dtmp=>tpful2d_head
!do while(associated(tpful2dtmp)) 
!   if(.not.associated(tpful2dtmp%next)) then
!      exit
!   else
!      print*, "rank=",rank,"coords=",tpful2dtmp%coords(1:2)
!      tpful2dtmp=>tpful2dtmp%next
!   end if
!end do
!print*, "rank2=",rank,num_p


rk4order=4

     do i=1, 10000
        if(rank==0) then
           print*, "#i=",i
        end if
 
        do j=1,pic2d%para2d%num_time 
           call tp_push_ful_orbit(tpful2d_head,numleft_boris,numgr,pic2d,"boris",rk4order,(i-1)*num_time+j)

           call tp_push_ful_orbit(tprk4ful2d_head,numleft_rk4,numgr,pic2d,"rk4", rk4order,(i-1)*num_time+j)

           call para_write_orbit_file_2d(tpful2d_head,numleft_boris,numgr,fileitem_boris,pic2d,(i-1)*num_time+j) 
           call para_write_orbit_file_2d(tprk4ful2d_head,numleft_rk4,numgr,fileitem_rk4,pic2d, (i-1)*num_time+j)
        end do

        call tp_push_gy_orbit_allmu(tpgy2dmu_head,numleft_gy,numgr_gy,pic2d,rk4order,i) 

        call para_write_orbit_file_2d_gy(tpgy2dmu_head,numleft_gy(1),numgr_gy,1,1,pic2d,i)
   end do

     call close_file(fileitem_boris,rank)
     call close_file(fileitem_rk4,rank)
!     do i=1,mu_num
       call close_file(1,rank)
!     end do

!      tpgy2dmutmp(1)%ptr=>tpgy2dmu_head(1)%ptr
!      do while(associated(tpgy2dmutmp(1)%ptr))
!        if(.not.associated(tpgy2dmutmp(1)%ptr%next)) then
!           exit
!        else
!           print*, "i=",i,"rank=",rank, tpgy2dmutmp(1)%ptr%coords(1:2)
!          tpgy2dmutmp(1)%ptr=>tpgy2dmutmp(1)%ptr%next
!        end if
!      end do




100  continue
       call mpi_barrier(comm)
     call MPI_FINALIZE(IERR) 
  end program  orbit_comparison

