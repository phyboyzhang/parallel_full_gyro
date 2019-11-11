module m_work_precision
  implicit none
  
  public :: f32, f64, i32, i64
    
  private

  integer, parameter :: i32 = selected_int_kind(9)

  integer, parameter :: i64 = selected_int_kind(18)

  integer, parameter :: f32 = selected_real_kind(1,37)

  integer, parameter :: f64 = selected_real_kind(1,99)


end module 
