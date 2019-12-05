program test_euler_maclaurin
#include "work_precision.h"
use utilities_module, only: muarray_eulermaclaurin, &
              muarray_euler_maclaurin_choice
implicit none

int4 :: mu_num=30
real8 :: mubound = 20.0
real8, dimension(:), pointer :: mus, muweight

allocate(mus(mu_num),muweight(mu_num))

!call muarray_eulermaclaurin(mubound,mu_num,mus,muweight)

call muarray_euler_maclaurin_choice(mubound,mu_num,mus,muweight,1)
print*, "mus", mus
print*, "muweight", muweight

end program test_euler_maclaurin
