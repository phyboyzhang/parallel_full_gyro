module gp_m_assert
  implicit none

  public :: gp_s_assertion


  private

    ! Instead of using the non-standard subroutine abort() provided by the compiler,
  ! use abort() from the C standard library "stdlib.h"
  
  interface
    subroutine c_abort() bind(C, name="abort")
    end subroutine
  end interface

contains

  subroutine gp_s_assertion(msg,file,line)
    character(len=*), intent(in) :: msg
    character(len=*), intent(in) :: file
    integer,  intent(in) :: line

    write(*,'(a)') "ASSERTION FAITURE: condition ("//trim(msg)//") is not satisfied."
    write(*,'(a,i0)') 'Triggered at'//file//':', line

    call c_abort()
  end subroutine gp_s_assertion

end module gp_m_assert
