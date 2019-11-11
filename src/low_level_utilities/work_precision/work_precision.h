#ifndef _work_precision_
#define _work_precision_

#define int4 integer(kind=4)
#define int8 integer(kind=8)

#define real4 real(kind=4)
#define real8 real(kind=8)


#define comp4 complex(kind=4)
#define comp8 complex(kind=8)

use m_work_precision
#define NULL_INT32 (/0.0_i32)
#define NULL_REAL64 (/0.0_:wf64/)


#endif
