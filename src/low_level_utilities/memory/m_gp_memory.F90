module gp_m_memory
  implicit none

  public :: gp_s_test_error_code


contains

  subroutine gp_s_test_error_code(err_code,descriptor,file_name,line_number)
    integer :: err_code
    character(len=*), intent(in) :: descriptor
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: line_number

    if(err_code .ne. 0) then
       write(*, '(a, a, i8)') descriptor, ' Triggered in FILE '//file_name// &
            ', in LINE: ', line_number
       stop 'ERROR: sll_s_test_error_code(): exiting program'
    end if       
  
  end subroutine gp_s_test_error_code






end module gp_m_memory
