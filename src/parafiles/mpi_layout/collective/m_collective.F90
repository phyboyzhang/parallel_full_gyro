module m_collective
#include "work_precision.h"
#include "gp_assert.h"

use utilities_module, only:&
  gp_error
!use mpi, only: &
!    mpi_allgather, &
!    mpi_allgatherv, &
!    mpi_allreduce, &
!    mpi_alltoall, &
!    mpi_alltoallv, &
!    mpi_barrier, &
!    mpi_bcast, &
!    mpi_comm_rank, &
!    mpi_comm_size, &
!    mpi_comm_split, &
!    mpi_comm_world, &
!    mpi_complex, &
!    mpi_double_complex, &
!    mpi_double_precision, &
!    mpi_finalize, &
!    mpi_gather, &
!    mpi_gatherv, &
!    mpi_init_thread, &
!    mpi_integer, &
!    mpi_integer8, &
!    mpi_logical, &
!    mpi_real, &
!    mpi_real8, &
!    mpi_reduce, &
!    mpi_scatter, &
!    mpi_scatterv, &
!    mpi_success, &
!    mpi_sum, &
!    mpi_thread_funneled
  
  implicit none
include "mpif.h"

  public :: &
       f_collectives_are_same, &
       o_collective_allgather,  &
       o_collective_allreduce,  &
       o_collective_alltoall,   &
       o_collective_alltoallv,  &
       mpi_collective_t,          &
       f_get_collective_rank,   &
       f_get_collective_size,   &
       f_new_collective

     private 
  
  type mpi_collective_t
     type(mpi_collective_t), pointer :: parent=>null()
     int4       :: comm  ! communicator
     int4       :: rank  !
     int4       :: key   ! control of the rank assignment
     int4       :: color !
     int4       :: size  ! communicator size
     int4       :: thread_level_required
     int4       :: thread_level_provided

  end type mpi_collective_t

  type(mpi_collective_t),pointer :: mpi_world_collective

  interface o_collective_allgather
     module procedure collective_allgather_int, &
            collective_allgather_real64
  end interface o_collective_allgather

!  interface o_collective_allgatherv
!     module procedure collective_allgatherv_real64
  !        collective_allgatherv_real64
!  end interface

  interface o_collective_allreduce
     module procedure collective_allreduce_real32, &
                      collective_allreduce_real64, &
                      collective_allreduce_comp64
  end interface
  
  
  interface o_collective_alltoall
     module procedure  collective_alltoall_int, &
          collective_alltoall_double, &
          collective_alltoall_complex_double
  end interface o_collective_alltoall

  interface o_collective_alltoallv
     module procedure collective_alltoallv_int, &
                      collective_alltoallv_double, &
                      collective_alltoallV_complex_double
  end interface  
  
contains


  subroutine check_collective_ptr(ptr)
    type(mpi_collective_t), pointer :: ptr
    if(.not.associated(ptr)) then
       write(*,'(a)') "check_collective_ptr: non-associated pointer"
       stop "check_collective_ptr"
    end if
  end subroutine check_collective_ptr

  subroutine test_mpi_error( ierr, descriptor )
    int4, intent(in)        :: ierr
    character(len=*), intent(in) :: descriptor

    if( ierr .ne. MPI_SUCCESS ) then
       write (*, '(a, a)') 'MPI error code failure: ', descriptor
    end if
  end subroutine test_mpi_error
  
  subroutine mpi_boot_collective(required)
    int4, intent(in), optional :: required
    int4 :: ierr

    allocate(mpi_world_collective, stat=ierr)

    if(present(required)) then
       mpi_world_collective%thread_level_required = required
    else
       mpi_world_collective%thread_level_required = MPI_THREAD_FUNNELED
    end if

    call MPI_Init_Thread(mpi_world_collective%thread_level_required, &
                    mpi_world_collective%thread_level_required, &
                    ierr)
    mpi_world_collective%comm    = MPI_COMM_WORLD
    mpi_world_collective%color   = 0
    mpi_world_collective%key     = 0
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_world_collective%rank, ierr)
    call test_mpi_error( ierr, 'boot_collective(): MPI_COMM_RANK()' )
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_collective%size, ierr )
    call test_mpi_error( ierr, 'boot_collective(): MPI_COMM_SIZE()')    
    
  end subroutine mpi_boot_collective

  subroutine mpi_halt_collective
    int4 :: ierr
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    call test_mpi_error( ierr, &
         'mpi_halt_collective(): MPI_BARRIER()' )
    call MPI_Finalize(ierr)
  end subroutine mpi_halt_collective
  

  function f_collectives_are_same( col1, col2 )
    logical :: f_collectives_are_same
    type(mpi_collective_t), pointer :: col1
    type(mpi_collective_t), pointer :: col2
    ASSERT( associated(col1) )
    ASSERT( associated(col2) )
    if( col1%comm == col2%comm ) then
       f_collectives_are_same = .true.
    else
       f_collectives_are_same = .false.
    end if
  end function f_collectives_are_same


  function f_get_collective_rank( col )
    type(mpi_collective_t), pointer :: col
    int4                       :: f_get_collective_rank
    call check_collective_ptr( col )
    f_get_collective_rank = col%rank
  end function f_get_collective_rank

  function f_get_collective_size( col )
    type(mpi_collective_t), pointer :: col
    int4                     :: f_get_collective_size
    call check_collective_ptr( col )
    f_get_collective_size = col%size
  end function f_get_collective_size

  function get_collective_parent( col )
    type(mpi_collective_t), pointer :: col
    type(mpi_collective_t), pointer :: get_collective_parent
    call check_collective_ptr( col )
    get_collective_parent => col%parent
  end function get_collective_parent  


  function f_new_collective( parent, color, key )
    type(mpi_collective_t), pointer  :: parent
    int4, intent(in)            :: color
    int4, intent(in)            :: key
    int4                        :: ierr
    type(mpi_collective_t), pointer  :: f_new_collective

    call check_collective_ptr( parent )
    ALLOCATE( f_new_collective, stat=ierr )
    call gp_error(ierr,"f_new_collective")
    f_new_collective%parent => parent
    f_new_collective%color  = color
    f_new_collective%key    = key
    call MPI_COMM_SPLIT( parent%comm, color, key, f_new_collective%comm, &
                         ierr)
    call test_mpi_error( ierr, 'f_new_collective(): MPI_COMM_SPLIT()')
    ! fill out the rest of the fields.
    call MPI_COMM_RANK(f_new_collective%comm, f_new_collective%rank, ierr)
    call test_mpi_error(ierr, 'f_new_collective(): MPI_COMM_RANK()')
    call MPI_COMM_SIZE(f_new_collective%comm, f_new_collective%size, ierr)
    call test_mpi_error(ierr, 'f_new_collective(): MPI_COMM_RANK()')
  end function f_new_collective  
  
  subroutine collective_allgather_int( col, send_buf, send_sz, &
       recv_buf, recv_sz )
    type(mpi_collective_t), pointer        :: col
    int4, dimension(:), intent(in)    :: send_buf ! what would change...
    int4, intent(in)                  :: send_sz
    int4, dimension(:), intent(inout) :: recv_buf ! would also change
    int4, intent(in)                  :: recv_sz
    int4                              :: ierr
    ! FIXME: Argument checking
    call check_collective_ptr( col )
    call MPI_ALLGATHER( send_buf(:), send_sz, MPI_INTEGER, &
                        recv_buf(:), recv_sz, MPI_INTEGER, col%comm, ierr )
    call test_mpi_error( ierr, &
         'ollective_allgather_int(): MPI_ALLGATHER()' )

  end subroutine collective_allgather_int

  subroutine collective_allgather_real64( &
    col, &
    send_buf, &
    send_sz, &
    recv_buf, &
    recv_sz )

    type(mpi_collective_t), pointer         :: col
    real8, dimension(:), intent(in)    :: send_buf
    int4, intent(in)                   :: send_sz
    real8, dimension(:), intent(out)   :: recv_buf ! would change
    int4, intent(in)                   :: recv_sz
    int4                               :: ierr
    ! FIXME: Argument checking
    call check_collective_ptr( col )
    call MPI_ALLGATHER( send_buf(:), send_sz, MPI_REAL8, &
                        recv_buf(:), recv_sz, MPI_REAL8, &
                        col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_allgather_int(): MPI_ALLGATHER()' )

  end subroutine collective_allgather_real64

  subroutine collective_allreduce_real32( col, send_buf, count, op, &
       rec_buf )
    type(mpi_collective_t), pointer       :: col
    real4, dimension(:), intent(in)  :: send_buf ! what would change...
    int4, intent(in)                 :: count
    int4, intent(in)                :: op
    real4, dimension(:), intent(out) :: rec_buf  ! would also change
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_REAL, &
      op, &
      col%comm, &
      ierr )
    call test_mpi_error( &
      ierr, &
      'collective_allreduce_real32(): MPI_ALLREDUCE()' )
  end subroutine collective_allreduce_real32  

  subroutine collective_allreduce_real64( col, send_buf, count, op, &
       rec_buf )
    type(mpi_collective_t), pointer       :: col
    real8, dimension(:), intent(in)  :: send_buf ! what would change...
    int4, intent(in)                 :: count
    int4, intent(in)                :: op
    real8, dimension(:), intent(out) :: rec_buf  ! would also change
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_DOUBLE_PRECISION, &
      op, &
      col%comm, &
      ierr )
    call test_mpi_error( &
      ierr, &
      'collective_allreduce_real64(): MPI_ALLREDUCE()' )
  end subroutine collective_allreduce_real64


  subroutine collective_allreduce_real64_2darray( col, send_buf, count, op, &
       rec_buf )
    type(mpi_collective_t), pointer       :: col
    real8, dimension(:,:), intent(in)  :: send_buf ! what would change...
    int4, intent(in)                 :: count
    int4, intent(in)                :: op
    real8, dimension(:,:), intent(out) :: rec_buf  ! would also change
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_BARRIER( col%comm, ierr )
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_DOUBLE_PRECISION, &
      op, &
      col%comm, &
      ierr )
    call test_mpi_error( &
      ierr, &
      'collective_allreduce_real64(): MPI_ALLREDUCE()' )
  end subroutine collective_allreduce_real64_2darray


  subroutine collective_allreduce_comp64( col, send_buf, count, op, &
       rec_buf )
    type(mpi_collective_t), pointer       :: col
    comp8, dimension(:), intent(in)  :: send_buf ! what would change...
    int4, intent(in)                 :: count
    int4, intent(in)                :: op
    comp8, dimension(:), intent(out) :: rec_buf  ! would also change
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLREDUCE( &
      send_buf, &
      rec_buf, &
      count, &
      MPI_DOUBLE_COMPLEX , &
      op, &
      col%comm, &
      ierr )
    call test_mpi_error( &
      ierr, &
      'collective_allreduce_comp64(): MPI_ALLREDUCE()' )
  end subroutine collective_allreduce_comp64

  subroutine collective_alltoall_int( send_buf, send_count, &
                                          recv_count, recv_buf, col )
    int4, dimension(:), intent(in)  :: send_buf
    int4, intent(in)                :: send_count
    int4, intent(in)                :: recv_count
    int4, dimension(:), intent(out) :: recv_buf
    type(mpi_collective_t), pointer      :: col
    int4                            :: ierr
    call check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_ALLTOALL( send_buf, send_count, MPI_INTEGER, &
                       recv_buf, recv_count, MPI_INTEGER, col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoall_int(): MPI_ALLTOALLV()' )

  end subroutine collective_alltoall_int

  subroutine collective_alltoall_double( send_buf, send_count, &
                                             recv_count, recv_buf, col )
    real8, dimension(:), intent(in)  :: send_buf
    int4, intent(in)                 :: send_count
    int4, intent(in)                 :: recv_count
    real8, dimension(:), intent(out) :: recv_buf
    type(mpi_collective_t), pointer       :: col
    int4                             :: ierr
    call check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_ALLTOALL( send_buf, send_count, MPI_DOUBLE_PRECISION, &
                       recv_buf, recv_count, MPI_DOUBLE_PRECISION, &
                       col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoall_double(): MPI_ALLTOALL()' )
  end subroutine collective_alltoall_double

  subroutine collective_alltoall_complex_double( send_buf, send_count, &
                                                     recv_count, recv_buf, col )
    comp8, dimension(:), intent(in)  :: send_buf
    int4, intent(in)                 :: send_count
    int4, intent(in)                 :: recv_count
    comp8, dimension(:), intent(out) :: recv_buf
    type(mpi_collective_t), pointer       :: col
    int4                             :: ierr
    call check_collective_ptr( col )
    ! FIXME: MORE ARG CHECKING
    call MPI_ALLTOALL( send_buf, send_count, MPI_DOUBLE_COMPLEX, &
                       recv_buf, recv_count, MPI_DOUBLE_COMPLEX, &
                       col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoall_complex_double(): MPI_ALLTOALL()' )

  end subroutine collective_alltoall_complex_double


  subroutine collective_alltoallv_int( send_buf, send_cnts, &
                                           send_displs, &
                                           recv_buf, recv_cnts, &
                                           recv_displs, col )
    int4, dimension(:), intent(in) :: send_buf
    int4, dimension(:), intent(in) :: send_cnts
    int4, dimension(:), intent(in) :: send_displs
    int4, dimension(:), intent(out) :: recv_buf
    int4, dimension(:), intent(in) :: recv_cnts
    int4, dimension(:), intent(in) :: recv_displs
    type(mpi_collective_t), pointer     :: col
    int4                           :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLTOALLV( send_buf(:), send_cnts(:), send_displs(:), MPI_INTEGER,&
         recv_buf(:), recv_cnts(:), recv_displs(:), MPI_INTEGER,&
         col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoallv_int(): MPI_ALLTOALLV()' )

 end subroutine collective_alltoallv_int  


  subroutine collective_alltoallv_double( send_buf, send_cnts, &
                                              send_displs, &
                                              recv_buf, recv_cnts, &
                                              recv_displs, col )
    real8, dimension(:), intent(in)  :: send_buf
    int4,  dimension(:), intent(in)  :: send_cnts
    int4,  dimension(:), intent(in)  :: send_displs
    real8, dimension(:), intent(out) :: recv_buf
    int4,  dimension(:), intent(in)  :: recv_cnts
    int4,  dimension(:), intent(in)  :: recv_displs
    type(mpi_collective_t), pointer       :: col
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLTOALLV( send_buf, send_cnts, send_displs, &
                        MPI_DOUBLE_PRECISION, &
                        recv_buf, recv_cnts, recv_displs, &
                        MPI_DOUBLE_PRECISION, &
                        col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoallv_double(): MPI_ALLTOALLV()' )

  end subroutine collective_alltoallv_double


 subroutine collective_alltoallV_complex_double( send_buf, send_cnts, &
                                                      send_displs, &
                                                      recv_buf, recv_cnts, &
                                                      recv_displs, col )
    comp8, dimension(:), intent(in)  :: send_buf
    int4,  dimension(:), intent(in)  :: send_cnts
    int4,  dimension(:), intent(in)  :: send_displs
    comp8, dimension(:), intent(out) :: recv_buf
    int4,  dimension(:), intent(in)  :: recv_cnts
    int4,  dimension(:), intent(in)  :: recv_displs
    type(mpi_collective_t), pointer       :: col
    int4                             :: ierr
    ! FIXME: ARG CHECKING!
    call check_collective_ptr( col )
    call MPI_ALLTOALLV( send_buf, send_cnts, send_displs, &
                        MPI_DOUBLE_COMPLEX, &
                        recv_buf, recv_cnts, recv_displs, &
                        MPI_DOUBLE_COMPLEX, &
                        col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_alltoallV_complex_double(): MPI_ALLTOALLV()' )

  end subroutine collective_alltoallV_complex_double


    subroutine collective_scatterv_real64( col, send_buf, send_count,&
                                        displs,recv_count, root,rec_buf)
    type(mpi_collective_t), pointer      :: col
    real8, dimension(:), intent(in) :: send_buf ! what would change...
    int4, dimension(:), intent(in)  :: send_count
    int4, dimension(:), intent(in)  :: displs
    int4, intent(in)                :: recv_count
    real8, dimension(:), intent(out) :: rec_buf  ! would also change
    int4, intent(in)                :: root
    int4                            :: ierr
    ! FIXME: ARG CHECKING!
    ASSERT( SIZE(send_count) .eq. col%size )
    ASSERT( SIZE(displs) .eq. col%size )
    call MPI_SCATTERV( send_buf, send_count, displs, MPI_REAL, rec_buf,&
         recv_count, MPI_REAL, root, col%comm, ierr )
    call test_mpi_error( ierr, &
         'collective_scatterv_real(): MPI_SCATTERV()' )
  end subroutine collective_scatterv_real64
  

  subroutine collective_allgatherv_real64( &
    col, &
    send_buf, &
    send_cnt, &
    rec_cnt, &
    displs, &
    rec_buf )

    type(mpi_collective_t), pointer       :: col
    real8, dimension(:), intent(in)  :: send_buf ! what would change...
    int4, intent(in)                 :: send_cnt
    int4, dimension(:), intent(in)   :: rec_cnt
    int4, dimension(:), intent(in)   :: displs
    real8, dimension(:), intent(out) :: rec_buf  ! would also change
    int4                            :: ierr
    ! FIXME: argument checking
    ASSERT(col%size .eq. SIZE(displs))
    ASSERT(col%size .eq. SIZE(rec_cnt))
    call MPI_ALLGATHERV( &
         send_buf, &
         send_cnt, &
         MPI_DOUBLE_PRECISION, &
         rec_buf, &
         rec_cnt, &
         displs, &
         MPI_DOUBLE_PRECISION, &
         col%comm, &
         ierr )
    call test_mpi_error( ierr, &
         'collective_allgatherv_real64(): MPI_ALLGATHERV()' )
  end subroutine collective_allgatherv_real64
  
end module m_collective


