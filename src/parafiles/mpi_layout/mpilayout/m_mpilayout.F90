!For periodic boundary condition, the maxiumal line in both dimension records the same values of the first line in both 
!dimensions.

module m_mpilayout
#include "work_precision.h"
#include "gp_assert.h"
  
  use m_collective, only: &
       f_collectives_are_same, &
       o_collective_allgather,  &
       o_collective_allreduce,  &
       o_collective_alltoall,   &
       o_collective_alltoallv,  &
       mpi_collective_t,        &
       f_get_collective_rank,   &
       f_get_collective_size,   &
       f_new_collective
 !      f_get_layout_3D_collective

  use utilities_module, only: &
       s_int2string, &
       f_is_even, &
       f_is_power_of_two, &
       gp_error

  implicit none
include "mpif.h"

  public :: t_layout_2d, t_layout_3d, &
            get_layout_2d_box_index,  &
            initialize_layout_with_distributed_2d_array, &
            split_domain_into_box1d
  private

  type :: box_2D
     int4  ::  i_min, i_max
     int4  ::  j_min, j_max
!!$     real8, dimension(:), pointer :: send_il
!!$     real8, dimension(:), pointer :: send_ir
!!$     real8, dimension(:), pointer :: send_jl
!!$     real8, dimension(:), pointer :: send_jr
     
!!$     real8, dimension(:), pointer :: recv_il
!!$     real8, dimension(:), pointer :: recv_ir
!!$     real8, dimension(:), pointer :: recv_jl
!!$     real8, dimension(:), pointer :: recv_jr
  end type box_2D

  type :: box_3D
     int4  ::  i_min, i_max
     int4  ::  j_min, j_max
     int4  ::  k_min, k_max
!!$     real8, dimension(:), pointer :: send_il
!!$     real8, dimension(:), pointer :: send_ir
!!$     real8, dimension(:), pointer :: send_jl
!!$     real8, dimension(:), pointer :: send_jr
!!$     real8, dimension(:), pointer :: send_kl
!!$     real8, dimension(:), pointer :: send_kr
!!$     
!!$     real8, dimension(:), pointer :: recv_il
!!$     real8, dimension(:), pointer :: recv_ir
!!$     real8, dimension(:), pointer :: recv_jl
!!$     real8, dimension(:), pointer :: recv_jr
!!$     real8, dimension(:), pointer :: recv_jl
!!$     real8, dimension(:), pointer :: recv_jr     
  end type box_3D

  type :: box_4d
     int4  ::  i_min, i_max
     int4  ::  j_min, j_max
     int4  ::  k_min, k_max
     int4  ::  l_min, l_max
  end type box_4d
  

  
  
!!$  type :: box_4D
!!$     int4, private  ::  i_min, i_max
!!$     int4, private  ::  j_min, j_max
!!$     int4, private  ::  k_min, k_max
!!$     int4, private  ::  l_min, l_max
!!$  end type box_4D  

  type :: t_layout_2d
     type(mpi_collective_t), pointer  :: collective
     int4                           :: global_sz1
     int4                           :: global_sz2
     type(box_2d), dimension(:), pointer :: boxes
  end type t_layout_2d

  type :: t_layout_3d
     type(mpi_collective_t), pointer :: collective
     int4                           :: global_sz1
     int4                           :: global_sz2
     int4                           :: global_sz3
     type(box_3d), dimension(:), pointer :: boxes
  end type t_layout_3d

   type :: t_layout_4d
     type(mpi_collective_t), pointer :: collective
     int4                          :: global_sz1
     int4                           :: global_sz2
     int4                           :: global_sz3
     int4                           :: global_sz4
     type(box_4d), dimension(:), pointer, private :: boxes
  end type t_layout_4d

  
  type(t_layout_2d), pointer :: t_layout_2d_ptr



contains

  
  function f_new_layout_2d( col )                       
    intrinsic :: associated                               
    type(t_layout_2d), pointer  :: f_new_layout_2d               
    type(mpi_collective_t), pointer :: col                 
    int4                            :: n_nodes             
    int4                            :: ierr                
    if(.not.associated(col) ) then                         
       write(*,'(a)') 'error, uninitialized collective'    
       stop 'new_layout_function'                          
    end if                                                 
    allocate( f_new_layout_2d, stat=ierr)                             
    f_new_layout_2d%collective => col                            
    n_nodes              = f_get_collective_size(col)      
    allocate( f_new_layout_2d%boxes(0:(n_nodes-1)), stat=ierr)         
  end function f_new_layout_2d

  
  subroutine delete_layout_2d( layout ) 
    type(t_layout_2d), pointer   :: layout
    int4                         :: ierr     
    nullify( layout%collective)            
    deallocate( layout%boxes, stat=ierr)   
    deallocate( layout, stat=ierr)         
  end subroutine delete_layout_2d
  

  function f_get_layout_2d_box( layout, rank )
    type(box_2d)             :: f_get_layout_2d_box
    type(t_layout_2d), pointer :: layout
    int4, intent(in)      :: rank
    ASSERT((rank.ge.0).and.(rank.le.(get_layout_num_nodes(layout)-1)))
    f_get_layout_2d_box = layout%boxes(rank)
  end function f_get_layout_2d_box

  
  subroutine get_layout_2d_box_index( layout, rank, index )
    type(t_layout_2d), pointer    :: layout                         
    int4, intent(in)      :: rank
    int4, dimension(4), intent(inout)   :: index
    
    index(1) = layout%boxes(rank)%i_min
    index(2) = layout%boxes(rank)%i_max
    index(3) = layout%boxes(rank)%j_min
    index(4) = layout%boxes(rank)%j_max    
  end subroutine get_layout_2d_box_index

  function get_box_2d_i_min(box)
    type(box_2d) :: box
    int4 get_box_2d_i_min
    get_box_2d_i_min=box%i_min
  end function get_box_2d_i_min

  function get_box_2d_i_max(box)
    type(box_2d) :: box
    int4 get_box_2d_i_max
    get_box_2d_i_max=box%i_max
  end function  

  function get_box_2d_j_min(box)
    type(box_2d) :: box
    int4 get_box_2d_j_min
    get_box_2d_j_min=box%j_min
  end function get_box_2d_j_min


  function get_box_2d_j_max(box)
    type(box_2d) :: box
    int4 get_box_2d_j_max
    get_box_2d_j_max=box%j_max
  end function  
  
  subroutine set_layout_2d_box( layout, rank, val )                          
    type(t_layout_2d), pointer    :: layout                         
    int4, intent(in)      :: rank                           
    int4,dimension(4), intent(in)      :: val                            
    layout%boxes(rank)%i_min = val(1)
    layout%boxes(rank)%i_max = val(2)
    layout%boxes(rank)%j_min = val(3)
    layout%boxes(rank)%j_max = val(4)      
  end subroutine set_layout_2d_box
  
  
  function get_layout_3d_global_size( layout ) result(res);                                  
    int4, dimension(2)                :: res                                  
    type(t_layout_2d), pointer :: layout                              
    res(1) = layout%global_sz1
    res(2) = layout%global_sz2
  end function get_layout_3d_global_size
  

  function get_layout_2d_collective( layout )
    intrinsic                       :: associated
    type(mpi_collective_t), pointer :: get_layout_2d_collective
    type(t_layout_2d), pointer      :: layout
    if( .not. associated(layout) ) then
       stop 'ERROR: uninitialized argument, get_layout_XD_collective()'
    end if
    get_layout_2d_collective => layout%collective
  end function get_layout_2d_collective
  

  function get_layout_2d_size( layout )
    intrinsic                  :: associated
    int4                  :: get_layout_2d_size
    type(t_layout_2d), pointer :: layout
    if( .not. associated(layout) ) then
       STOP 'ERROR: not associated argument passed to get_layout_size().'
    end if
    get_layout_2d_size = f_get_collective_size( layout%collective )
  end function get_layout_2d_size
  

  function linear_index_2D(npx1, i, j)
    int4, intent(in) :: npx1
    int4, intent(in) :: i
    int4, intent(in) :: j
    int4 :: linear_index_2D
    linear_index_2D = i + npx1*j
  end function linear_index_2D

  function linear_index_3D(npx1, npx2, i, j, k)
    int4, intent(in) :: npx1
    int4, intent(in) :: npx2
    int4, intent(in) :: i
    int4, intent(in) :: j
    int4, intent(in) :: k
    int4 :: linear_index_3D
    linear_index_3D = i + npx1*(j + npx2*k)
  end function linear_index_3D

  subroutine get_2d_box_index(box, index)
    type(box_2d) :: box
    int4, dimension(:), intent(inout) :: index

    index(1)=box%i_min
    index(2)=box%i_max
    index(1)=box%j_min
    index(2)=box%j_max
  end subroutine get_2d_box_index


  function local_to_global_2D( layout, doublet )
    int4, dimension(1:2)             :: local_to_global_2D
    type(t_layout_2d), pointer              :: layout
    int4, intent(in), dimension(1:2) :: doublet
    type(mpi_collective_t), pointer       :: col
    int4                             :: my_rank
    type(box_2D)                     :: box
    ! fixme: arg checking
    col                => get_layout_2D_collective( layout )
    my_rank            =  f_get_collective_rank( col )
    box                =  f_get_layout_2D_box( layout, my_rank )
    local_to_global_2D(1) = get_box_2D_i_min(box) + doublet(1) - 1
    local_to_global_2D(2) = get_box_2D_j_min(box) + doublet(2) - 1
  end function local_to_global_2D

!!$  function local_to_global_3D( layout, triplet )
!!$    int4, dimension(1:3)             :: local_to_global_3D
!!$    type(t_layout_3d), pointer            :: layout
!!$    int4, intent(in), dimension(1:3) :: triplet
!!$    type(mpi_collective_t), pointer       :: col
!!$    int4                             :: my_rank
!!$    type(box_3D)                          :: box
!!$    ! fixme: arg checking
!!$    col                => get_layout_3D_collective( layout )
!!$    my_rank            =  get_collective_rank( col )
!!$    box                =  get_layout_3D_box( layout, my_rank )
!!$    local_to_global_3D(1) = get_box_3D_i_min(box) + triplet(1) - 1
!!$    local_to_global_3D(2) = get_box_3D_j_min(box) + triplet(2) - 1
!!$    local_to_global_3D(3) = get_box_3D_k_min(box) + triplet(3) - 1
!!$  end function local_to_global_3D  

  subroutine initialize_layout_with_distributed_2d_array( &
       global_npx1, &
       global_npx2, &
       num_proc_x1, &
       num_proc_x2, &
!       num_send,    &
!       num_recv,    &
       layout, &
       boundary)

    int4, intent(in) :: global_npx1
    int4, intent(in) :: global_npx2
    int4, intent(in) :: num_proc_x1
    int4, intent(in) :: num_proc_x2
    character(len=*), intent(in) :: boundary
    type(t_layout_2d), pointer,intent(inout) :: layout
    int4 :: i,j
    int4 :: total_num_processors
    int4 :: node
    int4 :: collective_size
    int4 :: ierr
    int4, dimension(:,:), pointer :: box1d1
    int4, dimension(:,:), pointer :: box1d2

    int4 :: i_min
    int4 :: i_max
    int4 :: j_min
    int4 :: j_max
    
    int4 :: val(4)
    int4  :: err

! if( &
!       .not. f_is_power_of_two(int(num_proc_x1,i64)) .or. &
!       .not. f_is_power_of_two(int(num_proc_x2,i64)) ) then
!       print *, 'ERROR: distribute_2D_array() needs that the integers that',&
!            'describe the process mesh are powers of 2.'
!       STOP
!    end if

    if( &
       .not. (global_npx1 .gt. 0) .or. &
       .not. (global_npx2 .gt. 0) ) then
       print *, 'ERROR: distribute_2D_array() needs that the array dimensions',&
            'be greater than zero.'
       STOP
    end if
    
    select case(boundary)
      case("double_per")
        layout%global_sz1 = global_npx1-1
        layout%global_sz2 = global_npx2-1
       
      case("nat_per")
        layout%global_sz1 = global_npx1
        layout%global_sz2 = global_npx2-1  

      case default
        print*, "input the correct boundary condition. Boundary=", boundary
        stop
     end select
       



    ALLOCATE( box1d1(0:num_proc_x1-1,2), stat=err )
    ALLOCATE( box1d2(0:num_proc_x2-1,2), stat=err )

    ! Allocate the layout to be returned.    
    total_num_processors = num_proc_x1*num_proc_x2
!    collective_size = get_layout_2D_size(layout)
!    if( total_num_processors .ne. collective_size ) then
!       print *, 'ERROR, initialize_layout_with_distributed_2d_array(): ',&
!            'requested size of the processor mesh is inconsistent with ', &
!            'the size of the collective.', 'number of processors = ', &
!            total_num_processors, ' collective size = ', collective_size
!       STOP
!    end if

    ! Compute the arrays with the split index information along the different
    ! dimensions.
    box1d1 = split_domain_into_box1d( 1, global_npx1, num_proc_x1)
    box1d2 = split_domain_into_box1d( 1, global_npx2, num_proc_x2)

    ! Fill the layout array.
    do j=0, num_proc_x2-1
       j_min = box1d2(j,1)
       j_max = box1d2(j,2)
       do i=0, num_proc_x1-1
          node = linear_index_2D( num_proc_x1, i, j )
          i_min = box1d1(i,1)
          i_max = box1d1(i,2)
          val(1)=i_min
          val(2)=i_max
          val(3)=j_min
          val(4)=j_max
          call set_layout_2d_box( layout, node, val )
         
       end do
    end do   
    
    DEALLOCATE( box1d1, stat=err )
    DEALLOCATE( box1d2, stat=err )

  end subroutine initialize_layout_with_distributed_2d_array


  function split_domain_into_box1d( min, max, numproc) result(box1d)
    int4, intent(in)                       :: numproc
    int4, dimension(0:numproc-1,1:2) :: box1d
    int4, intent(in)                       :: min
    int4, intent(in)                       :: max
    int4                                   :: num_cells
    int4 :: num1,num2,i
    real8 :: x
    num_cells = max - min 
    if(num_cells.le.numproc) then
       print*, 'the number of the cells of the 1d domain must be equal to or large than  &
               the processors allocated for this dimension.'
    end if
    
    x=real(num_cells,8)/real(numproc,8)
    num1=floor(x)
    do while(num1*numproc.lt.num_cells) 
       num1=num1+1
    end do
    if(num1*(numproc-1)==num_cells) then
       print*, "The number of numproc should minus 1."
       stop
    end if
    num2=num_cells-num1*(numproc-1)
   
    do i=0, numproc-1
       box1d(i,1)=i*num1+min
       box1d(i,2)=(i+1)*num1+min
       if(i==numproc-1) then
          box1d(numproc-1,1)=(numproc-1)*num1+min
          box1d(numproc-1,2)=(numproc-1)*num1+min+num2
       end if
    end do

  end function split_domain_into_box1d

  function intersect_boxes_2D( b1, b2, ans )
    intrinsic                 :: min, max
    logical                   :: intersect_boxes_2D
    type(box_2D), intent(in)  :: b1, b2
    type(box_2D), intent(out) :: ans
    int4                 :: loi, hii
    int4                 :: loj, hij
    int4                 :: loib1, hiib1
    int4                 :: lojb1, hijb1
    int4                 :: loib2, hiib2
    int4                 :: lojb2, hijb2
    int4                 :: index_1(4), index_2(4)
    ! FIXME: add error checking, if boxes are null, for instance.
    call get_2d_box_index(b1,index_1)
    loib1 = index_1(1)
    hiib1 = index_1(2)
    lojb1 = index_1(3)
    hijb1 = index_1(4)

    call get_2d_box_index(b2,index_2)
    loib2 = index_2(1)
    hiib2 = index_2(2)
    lojb2 = index_2(3)
    hijb2 = index_2(4)

    ASSERT( (loib1 .le. hiib1) .and. (loib2 .le. hiib2) )
    ASSERT( (lojb1 .le. hijb1) .and. (lojb2 .le. hijb2) )

    loi = max(loib1, loib2)
    hii = min(hiib1, hiib2)
    loj = max(lojb1, lojb2)
    hij = min(hijb1, hijb2)

    if( (loi .gt. hii) .or. (loj .gt. hij) ) then 
       ans%i_min = 0
       ans%i_max = 0
       ans%j_min = 0
       ans%j_max = 0
       intersect_boxes_2D = .false.
    else
       ans%i_min = loi
       ans%i_max = hii
       ans%j_min = loj
       ans%j_max = hij
       intersect_boxes_2D = .true.
    end if
  end function intersect_boxes_2D  
  

  subroutine mpi_initialization(rank,size)
    int4, intent(in) :: rank,size
    int4 :: ierr

    call MPI_Init(ierr)
    call MPI_COMM_Size(mpi_comm_world,size,ierr)
    call MPI_COMM_Rank(mpi_comm_world,rank,ierr)

  end subroutine


end module m_mpilayout
