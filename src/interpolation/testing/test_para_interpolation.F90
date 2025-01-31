Program test_para_interpolation
#include "work_precision.h"
use constants,only: pi_
use cartesian_mesh,   only: cartesian_mesh_1d, &
                            init_para_cartesian_mesh_1d
use paradata_type, only: pic_para_2d_base
use paradata_layout, only: initialize_pic_para_2d_base
                              
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

use spline_module, only: compute_D_spl2D_per_per_noblock, &
                         s_contribution_spl, &
                         s_localize_new
use m_para_spline, only: para_compute_spl2D_weight, &
                         para_compute_spl2d_point_per_per_weight, &
                         para_compute_spl2d_field_point_per_per, &
                         para_spl2d_firstorder_derivatve_point_per_per

implicit none
include "mpif.h"

    class(pic_para_2d_base),pointer :: pic2d
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
    real8, dimension(:,:), pointer :: weight,val
    real8 :: deri_firstorder(2)
    int4 :: l,m,cell_per_unit(2), II(2),NC(2),flag
    real8 :: eta_star(2),eta_min(2),eta_max(2)
    

    allocate(weight(-1:2,-1:2),val(-1:2,-1:2))

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(num_p(0:size-1))
    num_p=0
    pic2d => initialize_pic_para_2d_base(size)
!!! initialize parameter_2d_sets
    pic2d%para2d%gxmin=(/0.0,0.0/)
    pic2d%para2d%gxmax=(/2.0*pi_,2.0*pi_/)
    pic2d%para2d%N_points=4
    pic2d%para2d%iter_number=20000
    pic2d%para2d%dtgy=0.5
    pic2d%para2d%num_time=20
    pic2d%para2d%boundary="double_per"
    pic2d%para2d%geometry="cartesian"
    pic2d%para2d%mu=0.2
    pic2d%para2d%row=3
    pic2d%para2d%cell_per_unit=(/20,20/) 
    row=pic2d%para2d%row
    cell_per_unit=pic2d%para2d%cell_per_unit

 
    !!! initialize layout2d      
    pic2d%layout2d%collective%rank=rank
    pic2d%layout2d%collective%size=size
    pic2d%layout2d%collective%comm=mpi_comm_world
    pic2d%para2d%numproc=NINT(sqrt(real(size,8)))
    numproc=pic2d%para2d%numproc 
    comm=pic2d%layout2d%collective%comm

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

   print*, "gboxmin(:,1)",pic2d%para2d%gboxmax(:,1)
end if

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


     allocate(pic2d%field2d%ep_weight(num1,num2),pic2d%field2d%epwg_w(num1,row),&
  pic2d%field2d%epwg_e(num1,row),pic2d%field2d%epwg_n(row,num2),pic2d%field2d%epwg_s(row,num2), &
  pic2d%field2d%epwg_sw(row,row),pic2d%field2d%epwg_se(row,row),pic2d%field2d%epwg_nw(row,row), &
  pic2d%field2d%epwg_ne(row,row))

    boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
    boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
    boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
    boxindex(4)=pic2d%layout2d%boxes(rank)%j_max 

  !!!periodic boundary condition
 
  do i=1,pic2d%para2d%m_x1%nodes
    do j=1, pic2d%para2d%m_x2%nodes
       globalind=globalind_from_localind_2d((/i,j/),pic2d%para2d%numproc,rank,pic2d%layout2d,pic2d%para2d%boundary)
       pic2d%field2d%ep(i,j)=cos(real(globalind(2)-1,8)*delta(2))    ! real(globalind(1)+globalind(2), 8)
    end do
  end do
!     if(rank==2) then
!        print*, "ep=", pic2d%field2d%ep
!    end if

call copy_boundary_value_per_per(pic2d%field2d%ep,rank,pic2d%para2d%numproc,pic2d%layout2d)

!     if(rank==2) then
!        print*, "ep=", pic2d%field2d%ep
!    end if

    rankone=compute_rank_from_globalind_2d(5,3,pic2d%para2d%numproc,pic2d%layout2d)
!    print*, "rank=",rank,"rankone=",rankone
  
!    print*, "rank=",rank, boxindex(:)
!    if(rank==0) then
!      do i=0,size-1
!          print*,"rank=",i,pic2d%para2d%gboxmin(i,:),pic2d%para2d%gboxmax(i,:)
!      end do
!    end if

   rootdata=>initialize_rootdata_structure(pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2)
    call gather_field_to_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,rank,size,boxindex,&
          pic2d%para2d%numproc,pic2d%layout2d) 

  if(rank==0) then
    call compute_D_spl2D_per_per_noblock( &
         pic2d%layout2d%global_sz1, &
         pic2d%layout2d%global_sz2, &
         rootdata%ASPL)
  end if


 !    if(rank==0) then
!       print*, 'rootfield=', rootdata%field
!    end if
!call mpi_barrier(comm)
!print*, rank

!    pic2d%field2d%ep=0._f64
!    call scatter_field_from_rootprocess_per_per(rootdata%field,pic2d%field2d%ep,size, &
!         pic2d%para2d%numproc,(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/),pic2d%layout2d)
!     
!    if(rank==0) then
!       print*, 'ep=', pic2d%field2d%ep
!    endif
!call mpi_barrier(comm)


!call mpi_barrier(comm)
!print*, rank

  call para_compute_spl2D_weight(rootdata%ASPL,rootdata%field,pic2d%field2d%ep_weight, &
       pic2d%para2d%numproc,pic2d%layout2d,pic2d%para2d%boundary)


!   print*, "rank=",rank,pic2d%field2d%ep_weight
!   pic2d%field2d%ep_weight=1.0 
!   call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
    call mpi2d_alltoallv_box_per_per(row,comm,rank,numproc,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w, &
       pic2d%field2d%epwg_e,pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw, &
       pic2d%field2d%epwg_se,pic2d%field2d%epwg_nw,pic2d%field2d%epwg_ne,boxindex)       

! print*, "rank=",rank,"weight=", pic2d%field2d%ep_weight
! print*, "epwg_w=", pic2d%field2d%epwg_w
! print*, "epwg_e=", pic2d%field2d%epwg_e
! print*, "epwg_n=", pic2d%field2d%epwg_n
! print*, "epwg_s=",  pic2d%field2d%epwg_s
! print*, "epwg_sw=", pic2d%field2d%epwg_sw
! print*, "epwg_se=", pic2d%field2d%epwg_se
! print*, "epwg_nw=", pic2d%field2d%epwg_nw
! print*, "epwg_ne=", pic2d%field2d%epwg_ne
! if(rank==rankone) then
!     print*, "rankone=",rank
!     print*, pic2d%para2d%m_x1%eta_min,pic2d%para2d%m_x1%eta_max
!     x=x1

     x=(/2.1324555367458369,3.0000000000000000/) 
!   x=(/1.0,2.0/)  
   x1=x
!    call compute_process_of_point_per_per(x,pic2d%para2d%numproc,pic2d%para2d%gxmin,pic2d%para2d%gxmax, &
!                                          pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
    rankone=compute_process_of_point_per_per(x,pic2d%para2d%numproc,pic2d%para2d%gxmin,pic2d%para2d%gxmax, &
         pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)

!  if(rank==rankone) then  
!  print*, "rankone=",rankone
!  endif

  if(rank==4) then
    x1(1)=pic2d%para2d%gboxmin(rank,1)
    x1(2)=pic2d%para2d%gboxmin(rank,2) 

    fieldvalue= para_compute_spl2d_field_point_per_per(x1,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
        row,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
       pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
       pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw)

print*, "field=",fieldvalue, pic2d%field2d%ep(1,1)

   call para_compute_spl2d_point_per_per_weight(weight,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
       x1,row,pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
       pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
       pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw)
!   print*, weight
    eta_min=(/pic2d%para2d%m_x1%eta_min,pic2d%para2d%m_x2%eta_min/)
    eta_max=(/pic2d%para2d%m_x1%eta_max,pic2d%para2d%m_x2%eta_max/) 
    NC=(/pic2d%para2d%m_x1%nodes,pic2d%para2d%m_x1%nodes/)
    call s_localize_new(x1,eta_min,eta_max,ii,eta_star,NC-(/1,1/),flag)
    call s_contribution_spl(eta_star,val) 
!            val=val/real(pic2d%para2d%N_points,8)
            fieldvalue=0._f64
            do m=-1,2
              do l=-1,2
                fieldvalue=fieldvalue+val(m,l)*weight(m,l)
              end do
            end do
    print*, "fieldvalue=",fieldvalue
    print*, "ep=", pic2d%field2d%ep(1,1)
  end if
!  c2d%para2d%gboxmin(rank,1)fieldvalue=para_compute_spl2d_field_point_per_per(x,pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
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

! end if


 
    call MPI_FINALIZE(IERR) 
end program test_para_interpolation
