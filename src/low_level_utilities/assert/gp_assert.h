#ifndef _assert_h_
#define _assert_h_

#if(defined (GFORTRAN) || defined (_PGI) || defined(MPIF90))
#define STRNG(x) "x"
#else
#define STRNG(x) #x
#endif

#ifdef DEBUG
#define ASSERT(x) if(.not.(x)) then; \
  call gp_s_assertion( STRNG(x), __FILE__, __LINE__); \
  end if;
#else
#define ASSERT(x)
#endif

! Version of GP_ASSERT that is also to run in release mode.
!#ifdef DEBUG
!#define ASSERT_ALWAY(x) if(.not.(x)) then; \
!  call gp_s_assertion( STRNG(x), __FILE__, __LINE__); \
!  end if;
!#else
!#define ASSERT_ALWAYS(x)
!#endif

use gp_m_assert

#endif
