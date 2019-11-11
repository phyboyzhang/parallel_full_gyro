#ifndef _memory_h_
#define _memory_h_

use gp_m_memory

#define gpallocate(array_name_and_lims, error_var)   \
  allocate(array_name_and_lims,stat=error_var);      \
  call gp_s_test_error_code(error_var, 'Memory allocation Failure.', __FILE__, \
  __LINE__);


#endif
