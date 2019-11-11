program mesh_test
#include "work_precision.h"
  use cartesian_mesh, only: cartesian_mesh_1d, &
                            f_new_cartesian_mesh_1d
!use cartesian_mesh, only : test
implicit none

int4 :: nodes = 4
real8 :: eta_min=0.0, eta_max = 4.0
type(cartesian_mesh_1d), pointer :: d_mesh

d_mesh=> f_new_cartesian_mesh_1d(nodes,eta_min,eta_max,"natural")
!  allocate(d_mesh)
!  d_mesh%num_cells=num_cells
  print*, d_mesh%nodes

end program
