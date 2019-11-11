program test_interpo_gradient
#include "work_precision.h"
  use spline_module, only: &
       compute_spl2d_double_per_weight 

  real8,dimension(:,:),pointer :: field_init, field_weight
  int4 :: nodes=4
  int4 :: i,j

  allocate(field_init(4,4))
  allocate(field_weight(4,4))

  do i=1,nodes
     do j=1, nodes
        field_init=1.0_f64
     end do
  end do

  call compute_spl2d_double_per_weight(field_init,field_weight,nodes,nodes)
                                   
  print*, field_weight(:,:) 
end program test_interpo_gradient

