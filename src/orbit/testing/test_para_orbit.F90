Program test_para_orbit
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                         pic_para_total2d_base
use paradata_layout,only:     initialize_pic_para_2d_base, &
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
                        ful2dsend_node,gy2d_node,gy2dsend_node, gy2dmu_node
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
                          mpi2d_alltoallv_send_particle_2d, &
                          sort_particles_among_ranks, &
                          mpi2d_alltoallv_send_particle_2d_gy, &
                          sort_particles_among_ranks_gy

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
                           para_write_field_file_2d
use m_para_orbit, only: borissolve, fulrk4solve, &
                        para_obtain_interpolation_elefield_per_per_ful, &
                        sortposition_by_process_ful2d, &
                        compute_f_of_points_out_orbit_ful2d, &
                        fulrkfunc_f_per_per, &
                        compute_f_of_points_in_orbit_ful2d
use orbit_data_base, only: rk4ful2dnode, pointin_node, rk4gy2dnode, &
                           tp_gy2d_node,tp_gy2dsend_node
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
    class(gy2d_node), pointer :: gy2d_head,gy2dtmp
    class(gy2dsend_node), dimension(:), pointer :: gy2dsend_head, gy2dsendtmp
    real8 :: coords(4)   
    int4 :: rank1
    int4, dimension(:), pointer :: num_gy 
!    real8 :: elef(3)

!!!!  outputdata 
    int4 :: fileitem1
    character(100) :: filepath1

!!!! Test para_orbit
   class(rk4ful2dnode), pointer :: partlist,partmp
   class(pointin_node), pointer :: inlist, intmp
   class(ful2d_node), pointer :: ful2dtmp,ful2d_head
   class(gy2dmu_node), dimension(:),  pointer :: gy2dmu_head,gy2dmutmp
   real8 :: f1(6), f2(6),f3(6),f4(6)
   real8 :: elef(3), magf(3)
   real8 :: rho,theta,x2(2)
!   int4 :: row,prank,size,numgr,rank
   int4, dimension(:), pointer :: num,recnum
   real8 :: vec(6),dt
   int4 :: h,numgr,numleft,numcircle,order,rk4order

    allocate(weight(-1:2,-1:2))

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(num_p(0:size-1),recnum(0:size-1),num_gy(0:size-1))
    allocate(pointhead(0:size-1))
    allocate(partlist,inlist)
    allocate(ful2d_head) 
    allocate(gy2d_head)  

    allocate(gy2dmu_head(1),gy2dmutmp(1))
    allocate(gy2dmu_head(1)%ptr)

!    allocate(numleft_gy(mu_num))

    do i=0, size-1
      allocate(pointhead(i)%ptr)
    end do

    num_p=0

    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/12.0,12.0/)
    pic2d%para2d%N_points=3
    pic2d%para2d%iter_number=500000
    pic2d%para2d%dtgy=1.0
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/10,10/) 
    pic2d%para2d%dtful=pic2d%para2d%dtgy/real(pic2d%para2d%num_time,8)
    !!! particle in cell part
    pic2d%para2d%numparticle=5000
    pic2d%para2d%numcircle = 3
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
       pic2d%field2d%bf03(i,j)=1.0
    end do
  end do

  call copy_boundary_value_per_per(pic2d%field2d%ep,rank,pic2d%para2d%numproc,pic2d%layout2d)
  call copy_boundary_value_per_per(pic2d%field2d%bf03,rank,pic2d%para2d%numproc,pic2d%layout2d)
  num_p=0

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

!if(rank==0) then
!   print*,"bf03wg=", pic2d%field2d%bf03wg
!end if

!  call mpi2d_alltoallv_box_per_per(row,comm,rank,pic2d%para2d%numproc, &
!       pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n,&
!       pic2d%field2d%bf03wg_s, pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
!       pic2d%field2d%bf03wg_nw,pic2d%field2d%bf03wg_ne,boxindex)
 
   fileitem1=10
   filepath1="/PARA/blsc950/electrostatic_exp/run/density.txt"
   num_p=0
   rho=sqrt(2.0*pic2d%para2d%mu)
   numcircle=pic2d%para2d%numcircle

   ful2dtmp=>ful2d_head
   gy2dmutmp(1)%ptr=>gy2dmu_head(1)%ptr
    allocate(currk(0:size-1))
    allocate(gy2dsendtmp(0:size-1),gy2dsend_head(0:size-1))
    do i=0,size-1
       currk(i)%ptr=>pic2d%ful2dsend_head(i)%ptr
       allocate(gy2dsend_head(i)%ptr)
       gy2dsendtmp(i)%ptr=>gy2dsend_head(i)%ptr
    end do
!   if(rank==1.or.rank==2.or.rank==4.or.rank==6.or.rank==5.or.rank==7) then
!     goto 10 
!   else

   gy2dmutmp(1)%ptr=>gy2dmu_head(1)%ptr
   num_gy=0

   do j=1,pic2d%para2d%m_x1%nodes-1
      do i=1,pic2d%para2d%m_x2%nodes-1
          do k=1,2
              if(k==1) then
                 x(1)=pic2d%para2d%m_x1%eta_min+(real(i-1,8)+0.5)*pic2d%para2d%m_x1%delta_eta
                 x(2)=pic2d%para2d%m_x2%eta_min+(real(j-1,8)+0.5)*pic2d%para2d%m_x2%delta_eta 
              else
                 x(1)=pic2d%para2d%m_x1%eta_min+(real(i-1+pic2d%para2d%m_x1%nodes-1,8)+0.5) &
                      *pic2d%para2d%m_x1%delta_eta
                 x(2)=pic2d%para2d%m_x2%eta_min+(real(j-1+pic2d%para2d%m_x2%nodes-1,8)+0.5) &
                      *pic2d%para2d%m_x2%delta_eta 
              end if

              coords(1)=x(1)
              coords(2)=x(2)
              coords(3)=0.0
              coords(4)=0.0
!              call coordinates_pointoutbound_per_per(coords(1:2),pic2d)
              rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
               pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
!if(rank==0) then
!print*, "rank1=",rank1
!endif
            if(rank1==rank) then 
              ful2dtmp%coords(1:4)=coords(1:4)
              num_p(rank1)=num_p(rank1)+1
              allocate(ful2dtmp%next)
              ful2dtmp=>ful2dtmp%next

              gy2dmutmp(1)%ptr%coords(1:2)=coords(1:2)
              gy2dmutmp(1)%ptr%coords(3)=pic2d%para2d%mu 
              num_gy(rank1)=num_gy(rank1)+1
              allocate(gy2dmutmp(1)%ptr%next)
              gy2dmutmp(1)%ptr=>gy2dmutmp(1)%ptr%next
            else
              currk(rank1)%ptr%coords(1)=coords(1)
              currk(rank1)%ptr%coords(2)=coords(2)
              currk(rank1)%ptr%coords(3)=coords(3)
              currk(rank1)%ptr%coords(4)=coords(4)
              num_p(rank1)=num_p(rank1)+1
              allocate(currk(rank1)%ptr%next)
              currk(rank1)%ptr=>currk(rank1)%ptr%next 


              gy2dsendtmp(rank1)%ptr%coords(1)=coords(1)
              gy2dsendtmp(rank1)%ptr%coords(2)=coords(2)
              gy2dsendtmp(rank1)%ptr%coords(3)=pic2d%para2d%mu
              num_gy(rank1)=num_gy(rank1)+1
              allocate(gy2dsendtmp(rank1)%ptr%next)
              gy2dsendtmp(rank1)%ptr=>gy2dsendtmp(rank1)%ptr%next
!print*,"rank2=",rank
            end if
          end do        
      end do
   end do
!10  end if
print*, rank
call mpi_barrier(comm)
   do i=0,size-1
     currk(i)%ptr=>pic2d%ful2dsend_head(i)%ptr
     gy2dsendtmp(i)%ptr=>gy2dsend_head(i)%ptr
   end do
   ful2dtmp=>ful2d_head
   gy2dtmp=>gy2d_head
   call mpi2d_alltoallv_send_particle_2d(ful2dtmp,currk,num_p,pic2d)

   call mpi2d_alltoallv_send_particle_2d_gy(gy2dmu_head,gy2dsend_head,num_gy,1,pic2d)



!   ful2dtmp=>ful2d_head
!    allocate(currk(0:size-1))
!    do i=0,size-1
!       currk(i)%ptr=>pic2d%ful2dsend_head(i)%ptr
!    end do
!    num_p=0
! 
!  if(rank==0) then
!   x1(1)=(pic2d%para2d%gxmin(1)+pic2d%para2d%gxmax(1))/2.0+0.5
!   x1(2)=(pic2d%para2d%gxmin(2)+pic2d%para2d%gxmax(2))/2.0
!
!   do j=0,numcircle-1
!       theta=real(j,8)*2.0_f64*pi_/real(numcircle,8)
!       coords(1)=x1(1)+rho*cos(theta)
!       coords(2)=x2(2)+rho*sin(theta)
!       coords(3)=rho*cos(theta+pi_/2.0_f64)
!       coords(4)=rho*sin(theta+pi_/2.0_f64)
!       call coordinates_pointoutbound_per_per(coords(1:2),pic2d)
!       rank1=compute_process_of_point_per_per(coords(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
!             pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
!            if(rank1==rank) then
!              ful2dtmp%coords(1:4)=coords(1:4)
!              num_p(rank)=num_p(rank)+1
!              allocate(ful2dtmp%next)
!              ful2dtmp=>ful2dtmp%next
!            else
!              currk(rank1)%ptr%coords(1)=coords(1)
!              currk(rank1)%ptr%coords(2)=coords(2)
!              currk(rank1)%ptr%coords(3)=coords(3)
!              currk(rank1)%ptr%coords(4)=coords(4)
!              num_p(rank1)=num_p(rank1)+1
!              allocate(currk(rank1)%ptr%next)
!              currk(rank1)%ptr=>currk(rank1)%ptr%next
!            end if
!     end do
!   end if
!
!   do i=0,size-1
!     currk(i)%ptr=>pic2d%ful2dsend_head(i)%ptr
!   end do
!   ful2dtmp=>ful2d_head
!   call mpi2d_alltoallv_send_particle_2d(ful2dtmp,currk,num_p,pic2d)

!ful2dtmp=>ful2d_head
!do while(associated(ful2dtmp)) 
!   if(.not.associated(ful2dtmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",ful2dtmp%coords(1:2)
!      ful2dtmp=>ful2dtmp%next
!   end if
!end do



  
!    dt=pic2d%para2d%dtful
!    numgr=4
!    recnum=0
!    partmp=>partlist
!    curful=>ful2d_head 
!    do while(associated(curful))
!      if(.not.associated(curful%next)) then
!        exit
!      else
!         partmp%coords(1:4)=curful%coords(1:4)
!         partmp%vec(1:2)=partmp%coords(1:2)
!         partmp%vec(4:5)=partmp%coords(3:4)
!         partmp%vec(3)=0.0
!         partmp%vec(6)=0.0
!         partmp%at=0
!         partmp%vec0=partmp%vec
!
!         curful=>curful%next
!         allocate(partmp%next)
!         partmp=>partmp%next
!      end if
!    end do
!        
!    num_p=0
!    recnum=0 
!    order=1
!
!  call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)


!   call sortposition_by_process_ful2d(num_p,recnum,partlist,inlist,pic2d,dt,order)
!
!partmp=>partlist
!do while(associated(partmp)) 
!   if(.not.associated(partmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"coords1=",partmp%coords(1:2)
!      partmp=>partmp%next
!   end if
!end do

!   deallocate(inlist)
!   nullify(inlist)
!   allocate(inlist)
!
!   num_p=0    
!   recnum=0
!   order=2 
!  call  sortposition_by_process_ful2d(num_p,recnum,partlist,inlist,pic2d,dt,order)
!  call  compute_f_of_points_out_orbit_ful2d(num_p,numgr,recnum,inlist,partlist,pic2d,order)
!  call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)
!
!
!   deallocate(inlist)
!   nullify(inlist)
!   allocate(inlist)
!   recnum=0
!   num_p=0   
!   order=3 
!  call  sortposition_by_process_ful2d(num_p,recnum,partlist,inlist,pic2d,dt,order)
!  call  compute_f_of_points_out_orbit_ful2d(num_p,numgr,recnum,inlist,partlist,pic2d,order)
!  call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)
!
!    deallocate(inlist)
!   nullify(inlist)
!   allocate(inlist)
!   recnum=0
!   num_p=0   
!   order=4 
!   call sortposition_by_process_ful2d(num_p,recnum,partlist,inlist,pic2d,dt,order)
!  call compute_f_of_points_out_orbit_ful2d(num_p,numgr,recnum,inlist,partlist,pic2d,order)
!  call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)

  rk4order=4

  do i=1, 100
!   if(rank==0) then
!   print*, "#iter_numeber=", i
!   end if
!   ful2dtmp=>ful2d_head 
   call fulrk4solve(ful2dtmp,pic2d,rk4order,i)
!   ful2dtmp=>ful2d_head
!   do j=0,size-1
!      nullify(pic2d%ful2dsend_head(j)%ptr%next) 
!      currk(j)%ptr=>pic2d%ful2dsend_head(j)%ptr
!   end do   
!   num_p=0
!   call sort_particles_among_ranks(ful2dtmp,currk,pic2d,num_p)
!   call mpi2d_alltoallv_send_particle_2d(ful2dtmp,currk,num_p,pic2d)

!   call gyrork4solve(gy2dtmp,pic2d)
!   do j=0, size-1
!      nullify(gy2dsend_head(j)%ptr%next)
!      gy2dsendtmp(j)%ptr=>gy2dsend_head(j)%ptr
!   end do 
!   num_gy=0
!   call sort_particles_among_ranks_gy(gy2d_head,gy2dsend_head,pic2d,num_gy)
!   call mpi2d_alltoallv_send_particle_2d_gy(gy2d_head,gy2dsend_head,num_gy,pic2d)
 end do
!print*, rank
 
!  call close_file(fileitem1,rank)

    call MPI_FINALIZE(IERR) 
  end program  test_para_orbit



