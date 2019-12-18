module piclayout
#include "work_precision.h"

  use utilities_module, only: gp_error
  use cartesian_mesh,   only: cartesian_mesh_1d
!  use m_mpilayout,       only: t_layout_2d
  implicit none


  public ::   parameters_set_2d, &
              parameters_array_2d, &
              ful2drank_node, &
              ful2d_node,  &
              ful2dsend_node, & 
              gy2d_node, &
              gy2dsend_node, &
              gy2drank_node,  & 
              gy2dmu_node,   &
              root_precompute_data,&
              pic_field_2d_base, &
              pic_field_2d, &
              initialize_rootdata_structure  

  private

  type :: parameters_set_2d
     type(cartesian_mesh_1d), pointer :: m_x1
     type(cartesian_mesh_1d), pointer :: m_x2
     int4  :: numproc(2)
     int4  :: num_time 
     real8 :: dtgy
     real8 :: dtful
     real8,pointer :: gboxmin(:,:),gboxmax(:,:)  ! global boundary of x1 and x2 for each process
     real8 :: sigma    ! the variance of the gaussian distribution
     character(99) :: geometry, boundary
     int8  :: numequ  ! number of particles for equilibrium sampling
     int8  :: numper       ! numbre of particles for the perturbation
     int8  :: numcircle    ! number of particles on larmor circle
                           ! homogeneously distributed on circle
     
     real8 :: mu
     int4  :: mu_num
     real8 :: mumin,mumax
     int4 :: mu_scheme
     int4 :: mu_tail    !! The number of particles located at the last mu number 
     int4 :: mulast    !! The number of the last mu number which is accepted.
     real8 :: tempt

     int4  :: row   !!! left,right,up,down, the number of lines for the communication for interpolation
     real8 :: gxmin(2),gxmax(2) 
     int4  :: N_points 
     int4  :: cell_per_unit(2)
     int4  :: iter_number
     int4  :: gyroorder   ! gyroorder=1, the fisrt order; gyroorder=2, the 2nd order 

!     real8, dimension(:),pointer :: mu_nodes,mu_weights
!     real8, dimension(:,:),pointer :: temp_i,temp_e
  end type parameters_set_2d

  type parameters_array_2d
     real8, dimension(:),pointer :: mu_nodes,mu_weights
     real8, dimension(:,:),pointer :: temp_i,temp_e
     int4, dimension(:), pointer :: munum_partition
  end type parameters_array_2d
!!                    north
!!              nn*************ee
!!              nn*************ee
!!              **$$$$$$$$$$$$$**
!!              **$           $**
!!        west  **$   box     $**  east      !! WE(6,2)
!!              **$           $**            !! NS(2,13)
!!              **$           $**
!!              **$$$$$$$$$$$$$**
!!              ww*************ss
!!              ww*************ss
!!                    south


  type :: pic_field_2d
    real8, dimension(:,:),    pointer :: Bf01,Bf02,Bf03 !! the magnetic field on the mesh
    real8, dimension(:,:),    pointer :: Bf01wg,Bf02wg,Bf03wg  !! the first index is for the three spatial dimensions
    real8, dimension(:,:),    pointer :: ep  !!  ep is the potential for the
                                             !!  full-particle simulation 
    real8, dimension(:,:),    pointer :: ep_weight
    real8, dimension(:,:),    pointer :: denf  !!! denf: the full density
    real8, dimension(:,:),    pointer :: denfeq 
    real8, dimension(:,:,:),  pointer :: deng,deng_weight  !!! deng : the gyro density
    real8, dimension(:,:,:),  pointer :: dengeq,dengeq_weight   !!! The first index is for the magnetic moment
    real8, dimension(:,:),    pointer :: dengtot,dengeqtot       !!! The total density

!   real8, dimension(:,:),    pointer :: dens_weight  
    real8, dimension(:,:),    pointer :: gep,gep_weight   !! gep is the ful potential for 
                                                          !! the gyrokinetic simulation 
    real8, dimension(:,:,:),    pointer :: epgyro
    real8, dimension(:,:,:),    pointer :: epgy_weight    
    real8, dimension(:,:),    pointer :: epgysq_weight
  end type  pic_field_2d
  
  type,extends(pic_field_2d) :: pic_field_2d_base
    real8, dimension(:,:),    pointer :: bf03_W,bf03_E, bf03_N,bf03_S,bf03_SW,bf03_SE,bf03_NE,bf03_NW
    real8, dimension(:,:),    pointer :: bf03wg_W,bf03wg_E, bf03wg_N,bf03wg_S,bf03wg_SW,bf03wg_SE,bf03wg_NE,bf03wg_NW    

    real8, dimension(:,:),    pointer :: ep_W,ep_E,ep_N,ep_S,ep_SW,ep_SE,ep_NE,ep_NW
    real8, dimension(:,:),    pointer :: epwg_W,epwg_E,epwg_N,epwg_S,epwg_SW,epwg_SE,epwg_NE,epwg_NW   

    real8, dimension(:,:),    pointer :: gep_W,gep_E,gep_N,gep_S,gep_SW,gep_SE,gep_NE,gep_NW
    real8, dimension(:,:),    pointer :: gepwg_W,gepwg_E,gepwg_N,gepwg_S,gepwg_SW,gepwg_SE,gepwg_NE,gepwg_NW

    real8, dimension(:,:),    pointer :: denf_W,denf_E,denf_N, denf_S, denf_SW, denf_SE, &
                                         denf_NE,denf_NW
    real8, dimension(:,:),    pointer :: denfeq_W,denfeq_E,denfeq_N, denfeq_S, denfeq_SW, denfeq_SE, &
                                         denfeq_NE,denfeq_NW

    real8, dimension(:,:,:),  pointer :: dengeq_W,dengeq_E,dengeq_N, dengeq_S, dengeq_SW, dengeq_SE, &
                                         dengeq_NE,dengeq_NW   
    real8, dimension(:,:,:),  pointer :: deng_W,deng_E,deng_N, deng_S, deng_SW, deng_SE, &
                                         deng_NE,deng_NW
    real8, dimension(:,:,:),  pointer :: dengeqwg_W,dengeqwg_E,dengeqwg_N, dengeqwg_S, dengeqwg_SW, dengeqwg_SE, &
                                         dengeqwg_NE,dengeqwg_NW
    real8, dimension(:,:,:),  pointer :: dengwg_W,dengwg_E,dengwg_N,dengwg_S,dengwg_SW, &
                                         dengwg_SE,dengwg_NE,dengwg_NW    
    real8, dimension(:,:),    pointer :: dengtot_W,dengtot_E,dengtot_N,dengtot_S,dengtot_SW,dengtot_SE, &
                                         dengtot_NE,dengtot_NW
    real8, dimension(:,:),    pointer :: dengeqtot_W,dengeqtot_E,dengeqtot_N,dengeqtot_S,dengeqtot_SW, &
                                         dengeqtot_SE,dengeqtot_NE,dengeqtot_NW
       
    real8, dimension(:,:,:),  pointer :: epgy_W,epgy_E,epgy_N,epgy_S,epgy_SW,epgy_SE,epgy_NE,epgy_NW
    real8, dimension(:,:,:),  pointer :: epgywg_W,epgywg_E,epgywg_N,epgywg_S,epgywg_SW,epgywg_SE, &
                                         epgywg_NE,epgywg_NW   
    real8, dimension(:,:),    pointer :: epgysqwg_W,epgysqwg_E,epgysqwg_N,epgysqwg_S,epgysqwg_SW,epgysqwg_SE,  &
                                         epgysqwg_NE,epgysqwg_NW

  end type pic_field_2d_base

  type :: root_precompute_data
    real8, dimension(:,:),    pointer :: ACONTRI,ASPL
    real8, dimension(:),      pointer :: field
    real8, dimension(:,:),    pointer :: prematrix
  end type

  type ful2d_node   !  no parallel velocity
     real8        :: coords(1:4)  ! two spatial coord and two velocity coord
     type(ful2d_node), pointer :: next=>null()
  end type ful2d_node

  type ful2drank_node    ! used to store particles running out of the current box
     real8 :: coords(1:4)
     int4  :: prank
     type(ful2drank_node), pointer :: next=>null()
  end type ful2drank_node

  type gy2d_node    ! no parallel velocity
     real8        :: coords(1:3)
     type(gy2d_node), pointer :: next=>null()
  end type gy2d_node

  type gy2drank_node
     real8 :: coords(1:3)
     int4  :: prank
     type(gy2drank_node), pointer :: next=>null()
  end type gy2drank_node
  
  type pic_sendrecv
     real8, dimension(:),   pointer :: sbuf    ! the sended particle information
     int4, dimension(:),    pointer :: scounts   ! for particle
     int4, dimension(:),    pointer :: sdispls   ! the displace for each process
     real8, dimension(:),   pointer :: rbuf    ! for particle
     int4, dimension(:),    pointer :: rcounts   ! for particle
     int4, dimension(:),    pointer :: rdispls   ! for particle
  end type pic_sendrecv


!  type pic_para_2d_base
!     type(t_layout_2d),       pointer  :: layout2d
!     type(pic_field_2d_base), pointer  :: field2d
!     type(parameters_set_2d),  pointer  :: para2d     
!  end type pic_para_2d_base

  type ful2dsend_node
     type(ful2drank_node),pointer :: ptr
  end type
  
!  type,extends(pic_para_2d_base) :: pic_para_total2d_base
!     type(ful2d_node),        pointer  :: ful2d_head ! store all the particle on one box
!     type(ful2dsend_node), dimension(:),allocatable :: ful2dsend_head  ! store the migration partilces
!     type(gy2d_node),         pointer  :: gy2d_head
!     type(gy2dsend_node), dimension(:), pointer  :: gy2dsend_head   
!   end type pic_para_total2d_base

  type gy2dsend_node
     type(gy2drank_node),pointer :: ptr
  end type

  type gy2dmu_node
     type(gy2d_node), pointer :: ptr
  end type  
  
!  type,extends(pic_para_2d_base) :: pic_para_gy2d_base  
!     type(gy2d_node),         pointer  :: gy2d_head
!     type(gy2dsend_node), dimension(:), pointer  :: gy2dsend_head  
!!     type(gy2d_node),         pointer  :: gy2d_send     
!  end type pic_para_gy2d_base


contains

!  function initialize_pic_para_ful2d_base() result(pic2d)
!    type(pic_para_total2d_base), pointer :: pic2d
!    int4 :: ierr
!    if(.not.associated(pic2d)) then
!       allocate(pic2d, stat=ierr)
!       call gp_error(ierr,"pic_2d")
!    end if
!
!!!$    allocate(pic_ful2d%ful2d_head,  stat=ierr)
!!!$    allocate(pic_ful2d%ful2dr_head, stat=ierr)
!!!$    allocate(pic_ful2d%layout2d,    stat=ierr)
!!!$    allocate(pic_ful2d%field2d,     stat=ierr)
!!!$    allocate(pic_ful2d%para2d,      stat=ierr)
!
!  end function initialize_pic_para_ful2d_base

!  function initialize_pic_para_gy2d_base() result(pic_gy2d)
!    type(pic_para_gy2d_base), pointer :: pic_gy2d
!    int4 :: ierr
!    if(.not.associated(pic_gy2d)) then
!       allocate(pic_gy2d, stat=ierr)
!       call gp_error(ierr,"pic_2d")
!    end if

!!$    allocate(pic_gy2d%gy2d_head,  stat=ierr)
!!$    allocate(pic_gy2d%gy2dr_head, stat=ierr)
!!$    allocate(pic_gy2d%layout2d,    stat=ierr)
!!$    allocate(pic_gy2d%field2d,     stat=ierr)
!!$    allocate(pic_gy2d%para2d,      stat=ierr)

!  end function initialize_pic_para_gy2d_base

   function initialize_rootdata_structure(dimsize) result(rootdata)
      type(root_precompute_data), pointer :: rootdata
      int4, intent(in) :: dimsize
      int4 :: ierr

      allocate(rootdata)
      allocate(rootdata%ACONTRI(dimsize,dimsize), stat=ierr)
      allocate(rootdata%ASPL(dimsize,dimsize), stat=ierr)
      allocate(rootdata%prematrix(dimsize,dimsize),stat=ierr)
      allocate(rootdata%field(dimsize), stat=ierr)

   end function


end module piclayout
