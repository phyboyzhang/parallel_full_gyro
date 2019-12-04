Program test_precomputing
#include "work_precision.h"
use constants, only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base, &
                              pic_para_total2d_base
use paradata_layout, only:    initialize_pic_para_2d_base, &
                              initialize_pic_para_total2d_base, &
                              allocate_memory_to_field_2d_ful, &
                              allocate_memory_to_field_2d_gy, &
                              allocate_memory_to_magfield_2d
                              
use utilities_module, only: f_is_power_of_two, &
                            gp_error
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
                                  store_data_on_rootprocess
use gyroaverage_2d_base, only: gyropoint_node
use m_fieldsolver, only: solve_weight_of_field_among_processes
use m_precompute, only: precompute_doublegyroaverage_matrix 

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
    int4 :: cell_per_unit(2)
 
    allocate(weight(-1:2,-1:2))

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(num_p(0:size-1))
    allocate(pointhead(0:size-1),stat=ierr)
    call gp_error(ierr,"pointhead")
    do i=0, size-1
      allocate(pointhead(i)%ptr,stat=ierr)
      call gp_error(ierr,"pointhead")
    end do

    num_p=0

    pic2d=> initialize_pic_para_total2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/2.0*pi_,2.0*pi_/)
    pic2d%para2d%N_points=3
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/10,10/) 
    row=pic2d%para2d%row

    !!! initialize layout2d      
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    numproc=pic2d%para2d%numproc 
    comm=pic2d%layout2d%collective%comm
    cell_per_unit=pic2d%para2d%cell_per_unit   
 
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
 
!     if(rank==2) then
!        print*, "ep=", pic2d%field2d%ep
!    end if

 call copy_boundary_value_per_per(pic2d%field2d%ep,rank,pic2d%para2d%numproc,pic2d%layout2d)

!     if(rank==2) then
!        print*, "ep=", pic2d%field2d%ep
!    end if

!    rankone=compute_rank_from_globalind_2d(5,3,pic2d%para2d%numproc,pic2d%layout2d)
!    print*, "rank=",rank,"rankone=",rankone
  
!    print*, "rank=",rank, boxindex(:)
!    if(rank==0) then
!      do i=0,size-1
!          print*,"rank=",i,pic2d%para2d%gboxmin(i,:),pic2d%para2d%gboxmax(i,:)
!      end do
!    end if

     x=(/2.2,1.3/) 
     x1=x
!    rankone=compute_process_of_point(x,pic2d%para2d%numproc,pic2d%para2d%gxmin,pic2d%para2d%gxmax, &
!         pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)

     rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)

    call gather_field_to_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,rank,size,boxindex,&
          pic2d%para2d%numproc,pic2d%layout2d) 

!    if(rank==0) then
!       print*, 'rootfield=', rootdata%field
!    end if

!    pic2d%field2d%ep=0._f64
!    call scatter_field_from_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,size, &
!         pic2d%para2d%numproc,(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/),pic2d%layout2d)
    
!    if(rank==0) then
!       print*, 'ep=', pic2d%field2d%ep
!    endif


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

!  call para_compute_spl2D_weight(rootdata%ASPL,rootdata%field,pic2d%field2d%ep_weight, &
!       pic2d%para2d%numproc,pic2d%layout2d,pic2d%para2d%boundary)
!
!
!
!!   print*, "rank=",rank,pic2d%field2d%ep_weight
! 
!   call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
!    call mpi2d_alltoallv_box_per_per(row,comm,rank,numproc,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w, &
!       pic2d%field2d%epwg_e,pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw, &
!       pic2d%field2d%epwg_se,pic2d%field2d%epwg_nw,pic2d%field2d%epwg_ne,boxindex)       

! if(rank==rankone) then
!!     print*, "rankone=",rank
!!     print*, pic2d%para2d%m_x1%eta_min,pic2d%para2d%m_x1%eta_max
!!     x=x1
!    call para_compute_spl2d_point_per_per_weight(weight,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
!       x,row,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
!       pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
!       pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw)
!
!  fieldvalue=para_compute_spl2d_field_point_per_per(x,pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
!             pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
!             pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
!             pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw)
!
!  print*, "fieldvalue=",fieldvalue
!
!
!  call para_spl2d_firstorder_derivatve_point_per_per(deri_firstorder,x,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
!        row, pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
!        pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
!        pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw)
!
! print*, "deri_firstorder=",deri_firstorder
!
! end if
!
! call mpi_barrier(comm)

!  num_p=0
!! if(rank==2) then
!  call sort_quadraturepoint_among_process(pic2d%para2d%mu,rank,pointhead,num_p,pic2d%para2d%N_points, & 
!          pic2d)
!if(rank==0) then
!  do i=0, size-1
!     print*, "rank=",rank,"i=",i,"num_p(i)=",num_p(i)
!  end do
! end if
!
!!num_p=0
!!if(rank==0) then
!  call para_compute_gyroaverage_stiff_matrix(pointhead,num_p,pic2d%para2d%mu,1,pic2d%para2d%N_points,pic2d,rootdata)
!end if
!if(rank==0) then
!
!print*, "acontri=",rootdata%acontri
!end if

!  call store_data_on_rootprocess(pic2d%para2d%mu,1,rank,rootdata,pic2d)

   call precompute_doublegyroaverage_matrix(rootdata,pic2d,pamearray)  

    call MPI_FINALIZE(IERR) 
  end program test_precomputing
