Program integ_simulation
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_structure, only: pic_para_2d_base, &
                              pic_para_total2d_base, &
                              initialize_pic_para_2d_base, &
                              initialize_pic_para_total2d_base, &
                              allocate_memory_to_field_2d
                              
use utilities_module, only: f_is_power_of_two
use m_mpilayout, only : initialize_layout_with_distributed_2d_array, &
                        get_layout_2d_box_index
use piclayout, only :   root_precompute_data, &
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
use m_fieldsolver, only: solve_weight_of_field_among_processes
use m_precompute,  only: precompute_ASPL
use m_moveparticles, only: new_position_per_per_ful, &
                           push_particle_among_box_ful2d_per_per, &
                           new_list_ful2d_per_per
use para_write_file, only: open_file,&
                           close_file, &
                           para_write_field_file_2d, &
                           para_write_orbit_file_2d, &
                           para_write_orbit_file_2d_gy
use m_para_orbit, only: borissolve, fulrk4solve, gyrork4solve, &
                        para_obtain_interpolation_elefield_per_per_ful, &
                        sortposition_by_process_ful2d, &
                        compute_f_of_points_out_orbit_ful2d, &
                        fulrkfunc_f_per_per, &
                        compute_f_of_points_in_orbit_ful2d
use orbit_data_base, only: rk4ful2dnode, pointin_node,tp_ful2d_node, &
                           tp_ful2dsend_node, &
                           rk4gy2dnode, tp_gy2d_node,tp_gy2dsend_node
use field_initialize,only: para_initialize_field_2d_mesh
use m_tp_para_orbit, only: tp_push_ful_orbit, &
                           tp_push_gy_orbit, &
                           tp_ful_solve, &
                           tp_ful2dlist_to_ful2dlist, &
                           tp_sort_particles_among_ranks
implicit none
include "mpif.h"

    class(pic_para_total2d_base),pointer :: pic2d
 !   class(pic_para_total2d_base)
    class(root_precompute_data), pointer :: rootdata
    real8 :: gxmin(2),gxmax(2)
    int4 :: cell_per_unit(2)
    real8 :: mumin,mumax
    int4 :: gyroorder
    int4 :: numcircle
    real8 :: dtgy,dtful
    int4 :: num_time,iter_num
    int4 :: numparticle
    character(20) :: boundary,geometry
    int4 :: row


    allocate(weight(-1:2,-1:2))

    call mpi_initialization(rank,size)
    pic2d%layout2d%collective%thread_level_provided=1
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
 
    namelist /mesh/  &
       gxmin,   &
       gxmax,   &
       cell_per_unit
    namelist /gyro/ &
       mumin,   &
       mumax,   &
       N_points,  &
       gyroorder
    namelist  /test_partiles/ &
       numcircle
    namelist /time_step/  &
       dtgy,  &
       dtful,  &
       num_time, &
       iter_num
    namelist  /pic/  &
       numparticle
    namelist  /profile/
    namelist  /conditions/ &
       boundary, &
       geometry
    namelist  /parallel/ &
       row

  call get_command_argument(1,filename)

  open(unit = input_file, file=trim(filename),IOstat=IO_stat)
  if(IO_stat /= 0) then
     print*, "#it's failed to open the initial data input file"
     stop
  end if

  read(input_file, mesh)
  read(input_file, gyro)
  read(input_file, test_particles)
  read(input_file, time_step)
  read(input_file, pic)
  read(input_file, profile)
  read(input_file, conditions)
  read(input_file, parallel)

  close(input_file) 






    allocate(num_p(0:size-1),recnum(0:size-1))
    allocate(pointhead(0:size-1))
    allocate(partlist,inlist) 

    allocate(tpful2d_head,tprk4ful2d_head,tpgy2d_head)
    allocate(tpful2dsend_head(0:size-1),tpful2dsendtmp(0:size-1))
    allocate(tprk4ful2dsend_head(0:size-1), tprk4ful2dsendtmp(0:size-1))
    allocate(tpgy2dsend_head(0:size-1), tpgy2dsendtmp(0:size-1))
    allocate(num_rk4(0:size-1),num_gy(0:size-1))
    allocate(ful2d_head, gy2d_head)
  

    do i=0, size-1
      allocate(pointhead(i)%ptr)
    end do

    num_p=0

    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=gxmin
    pic2d%para2d%gxmax=gxmax
    pic2d%para2d%N_points=N_points
    pic2d%para2d%iter_number=iter_number
    pic2d%para2d%numcircle=numcircle
    pic2d%para2d%dtgy=dtgy
    pic2d%para2d%num_time=num_time
    pic2d%para2d%boundary=boundary
    pic2d%para2d%geometry=geometry
    pic2d%para2d%row=row
    pic2d%para2d%cell_per_unit=cell_per_unit 
    pic2d%para2d%dtful=pic2d%para2d%dtgy/real(pic2d%para2d%num_time,8)
    !!! particle in cell part
    pic2d%para2d%sigma =
    pic2d%para2d%tempt = 
    pic2d%para2d%mumin=mumin
    pic2d%para2d%mumax=mumax
    pic2d%para2d%gyroorder=gyroorder

    call initialize_parameters2d(pic2d)

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

    num1=pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1
    num2=pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1                    

    call allocate_memory_to_field_2d(pic2d%field2d,num1,num2,row)
    rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)

    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

  !!!periodic boundary condition
!   call para_initialize_field_2d_mesh(amp,wave_one,wave_two, pic2d)
    dimsize=dimsize_of_rank_per_per(rank,pic2d%para2d%numproc,pic2d%layout2d)
    do i=1,dimsize(1)
      do j=1, dimsize(2)
        globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
        pic2d%field2d%ep(i,j)=1.0    ! real(globalind(1)+globalind(2), 8)
        pic2d%field2d%Bf03(i,j)=1.0
      end do
    end do

    if(rank==0) then
      call compute_D_spl2D_per_per_noblock( &
         pic2d%layout2d%global_sz1, &
         pic2d%layout2d%global_sz2, &
         rootdata%ASPL)
    end if
    call solve_weight_of_field_among_processes(pic2d%field2d%ep,rootdata%ASPL,rootdata,pic2d, &
       pic2d%field2d%ep_weight, pic2d%field2d%epwg_w,pic2d%field2d%epwg_e,pic2d%field2d%epwg_n, &
       pic2d%field2d%epwg_s, pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se, &
       pic2d%field2d%epwg_nw,pic2d%field2d%epwg_ne)

    call solve_weight_of_field_among_processes(pic2d%field2d%Bf03,rootdata%ASPL,rootdata,pic2d, &
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
   tpgy2dtmp=>tpgy2d_head
   tprk4ful2dtmp=>tprk4ful2d_head
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
         tpgy2dtmp%coords(1:2)=x1(1:2)
         tpgy2dtmp%coords(3)  =pic2d%para2d%mu
         tpgy2dtmp%tp=1
         num_gy(rank1)=num_gy(rank1)+1
         allocate(tpgy2dtmp%next)
         tpgy2dtmp=>tpgy2dtmp%next
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
       coords(2)=x2(2)+rho*sin(theta)
       coords(3)=rho*cos(theta+pi_/2.0_f64)
       coords(4)=rho*sin(theta+pi_/2.0_f64)
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
    call tp_mpi2d_alltoallv_send_particle_2d_gy(tpgy2d_head,numleft_gy,tpgy2dsend_head, &
                                                num_gy,numgr_gy,pic2d)

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
   fileitem_boris=10
   fileitem_rk4  =20
   fileitem_gy   =30
   filepath1="/home/qmlu/zsx163/electrostatic_exp/run/orbit_"

   filepath_boris=trim(filepath1)//"boris"//".txt"
   filepath_rk4  =trim(filepath1)//"rk4"//".txt"
   filepath_gy   =trim(filepath1)//"gy"//".txt"

     call open_file(fileitem_boris,filepath_boris,rank)
     call open_file(fileitem_rk4,  filepath_rk4,  rank)
     call open_file(fileitem_gy,   filepath_gy,   rank)
     
     
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

     do i=1, 1000
if(rank==0) then
print*, "#i=",i
end if
!         call tp_push_ful_orbit(tpful2d_head,numleft_boris,numgr,pic2d,"boris",rk4order,i)
!         call tp_push_ful_orbit(tprk4ful2d_head,numleft_rk4,numgr,pic2d,"rk4", rk4order,i)
         call tp_push_gy_orbit(tpgy2d_head,numleft_gy,numgr_gy,pic2d,rk4order,i) 

!         call para_write_orbit_file_2d(tpful2d_head,numleft_boris,numgr,fileitem_boris,pic2d,i) 
!         call para_write_orbit_file_2d(tprk4ful2d_head,numleft_rk4,numgr,fileitem_rk4,pic2d,i)
!         call para_write_orbit_file_2d_gy(tpgy2d_head,numleft_gy,numgr_gy,fileitem_gy,pic2d,i)
call mpi_barrier(pic2d%layout2d%collective%comm)
     end do

     call close_file(fileitem_boris,rank)
     call close_file(fileitem_rk4,rank)
     call close_file(fileitem_gy,rank)



   !  print*, rank
100  continue
     call MPI_FINALIZE(IERR) 
  end program  integ_simulation

