module paradata_layout
#include "work_precision.h"
use cartesian_mesh, only: init_para_cartesian_mesh_1d
use utilities_module, only: muarray_euler_maclaurin_choice
use paradata_type, only: pic_para_2d_base, &
                         pic_para_total2d_base
use piclayout, only: pic_field_2d_base, parameters_array_2d
use paradata_utilities, only: startind_of_process                         
implicit none

public::  initialize_pic_para_2d_base, &
          initialize_pic_para_total2d_base, &
          allocate_memory_to_field_2d_ful, &
          allocate_memory_to_field_2d_gy,  &
          allocate_memory_to_magfield_2d,  &
          allocate_parameters_array_2d,  &
          initialize_parameters_array_2d, &
          computing_mu_number
contains

  function initialize_pic_para_2d_base(size)   result(pic2d)
    type(pic_para_2d_base), pointer :: pic2d
    int4,intent(in) ::size
    allocate(pic2d)
    allocate(pic2d%layout2d)
    allocate(pic2d%layout2d%boxes(0:size-1))
    allocate(pic2d%layout2d%collective)
    allocate(pic2d%field2d)
    allocate(pic2d%para2d)
    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))
  end function


  function initialize_pic_para_total2d_base(size)  result(pic2d)
    class(pic_para_total2d_base), pointer :: pic2d
    int4, intent(in) :: size
    int4 :: i
    allocate(pic2d)
    allocate(pic2d%layout2d)
    allocate(pic2d%layout2d%boxes(0:size-1))
    allocate(pic2d%layout2d%collective)
    allocate(pic2d%field2d)

    allocate(pic2d%ful2d_head)
!    allocate(pic2d%ful2dsend_head(0:size-1))
    allocate(pic2d%gy2d_head)
!    allocate(pic2d%gy2dsend_head(0:size-1))

    allocate(pic2d%para2d)
    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))

!    do i=1,size
!      allocate(pic2d%ful2dsend_head(i-1)%ptr)
!      allocate(pic2d%gy2dsend_head(i-1)%ptr)
!    end do

  end function initialize_pic_para_total2d_base

  function allocate_parameters_array_2d(mu_num,global_sz) result(pamearray)
     class(parameters_array_2d), pointer :: pamearray
     int4, intent(in) :: mu_num, global_sz(2)

     allocate(pamearray)
     allocate(pamearray%temp_i(global_sz(1),global_sz(2)))
     allocate(pamearray%temp_e(global_sz(1),global_sz(2)))
     allocate(pamearray%mu_nodes(mu_num))
     allocate(pamearray%mu_weights(mu_num))
     allocate(pamearray%munum_partition(mu_num))

  end function allocate_parameters_array_2d


  subroutine allocate_memory_to_field_2d_ful(field2d,num1,num2,row)
    type(pic_field_2d_base), pointer,intent(inout) :: field2d
    int4,intent(in) :: num1,num2,row
    int4 :: ierr

   if(.not.associated(field2d)) then
      allocate(field2d,stat=ierr)
   end if

   allocate(field2d%ep(num1,num2),field2d%ep_w(num1,row), &
  field2d%ep_e(num1,row),field2d%ep_n(row,num2),field2d%ep_s(row,num2), &
  field2d%ep_sw(row,row),field2d%ep_se(row,row),field2d%ep_nw(row,row),field2d%ep_ne(row,row),stat=ierr)

  field2d%ep=0.0 
  field2d%ep_w=0.0
  field2d%ep_e=0.0
  field2d%ep_n=0.0
  field2d%ep_s=0.0
  field2d%ep_sw=0.0
  field2d%ep_se=0.0
  field2d%ep_nw=0.0
  field2d%ep_ne=0.0

   allocate(field2d%ep_weight(num1,num2),field2d%epwg_w(num1,row),&
  field2d%epwg_e(num1,row),field2d%epwg_n(row,num2),field2d%epwg_s(row,num2),&
  field2d%epwg_sw(row,row),field2d%epwg_se(row,row),field2d%epwg_nw(row,row),&
  field2d%epwg_ne(row,row), stat=ierr)

  field2d%ep_weight=0.0
  field2d%epwg_w=0.0 
  field2d%epwg_e=0.0
  field2d%epwg_n=0.0
  field2d%epwg_s=0.0
  field2d%epwg_sw=0.0
  field2d%epwg_se=0.0
  field2d%epwg_nw=0.0
  field2d%epwg_ne=0.0  


     allocate(field2d%denf(num1,num2),field2d%denf_w(num1,row),&
  field2d%denf_e(num1,row),field2d%denf_n(row,num2),field2d%denf_s(row,num2), &
  field2d%denf_sw(row,row),field2d%denf_se(row,row),field2d%denf_nw(row,row), &
  field2d%denf_ne(row,row), stat=ierr)

  field2d%denf=0.0
  field2d%denf_w=0.0
  field2d%denf_e=0.0
  field2d%denf_n=0.0
  field2d%denf_s=0.0
  field2d%denf_sw=0.0
  field2d%denf_se=0.0
  field2d%denf_nw=0.0
  field2d%denf_ne=0.0
 

     allocate(field2d%denfeq(num1,num2),field2d%denfeq_w(num1,row), &
  field2d%denfeq_e(num1,row),field2d%denfeq_n(row,num2),field2d%denfeq_s(row,num2), &
  field2d%denfeq_sw(row,row),field2d%denfeq_se(row,row),field2d%denfeq_nw(row,row), &
  field2d%denfeq_ne(row,row), stat=ierr)

  field2d%denfeq=0.0
  field2d%denfeq_w=0.0
  field2d%denfeq_e=0.0
  field2d%denfeq_n=0.0
  field2d%denfeq_s=0.0
  field2d%denfeq_sw=0.0
  field2d%denfeq_se=0.0
  field2d%denfeq_nw=0.0
  field2d%denfeq_ne=0.0

  end subroutine allocate_memory_to_field_2d_ful


  subroutine allocate_memory_to_field_2d_gy(field2d,num1,num2,row,mu_num)
    type(pic_field_2d_base), pointer,intent(inout) :: field2d
    int4,intent(in) :: num1,num2,row,mu_num
    int4 :: ierr

   if(.not.associated(field2d)) then
      allocate(field2d,stat=ierr)
   end if


   allocate(field2d%gep(num1,num2),field2d%gep_w(num1,row), &
  field2d%gep_e(num1,row),field2d%gep_n(row,num2),field2d%gep_s(row,num2), &
  field2d%gep_sw(row,row),field2d%gep_se(row,row),field2d%gep_nw(row,row), &
  field2d%gep_ne(row,row),stat=ierr)

  field2d%gep=0.0
  field2d%gep_w=0.0
  field2d%gep_e=0.0
  field2d%gep_n=0.0
  field2d%gep_s=0.0
  field2d%gep_sw=0.0
  field2d%gep_se=0.0
  field2d%gep_nw=0.0
  field2d%gep_ne=0.0   


   allocate(field2d%gep_weight(num1,num2),field2d%gepwg_w(num1,row), &
  field2d%gepwg_e(num1,row),field2d%gepwg_n(row,num2),field2d%gepwg_s(row,num2),&
  field2d%gepwg_sw(row,row),field2d%gepwg_se(row,row),field2d%gepwg_nw(row,row),&
  field2d%gepwg_ne(row,row), stat=ierr)

  field2d%gep_weight=0.0
  field2d%gepwg_w=0.0
  field2d%gepwg_e=0.0
  field2d%gepwg_n=0.0
  field2d%gepwg_s=0.0
  field2d%gepwg_sw=0.0
  field2d%gepwg_se=0.0  
  field2d%gepwg_nw=0.0
  field2d%gepwg_ne=0.0
 

   allocate(field2d%epgyro(mu_num,num1,num2),field2d%epgy_w(mu_num,num1,row), &
  field2d%epgy_e(mu_num,num1,row),field2d%epgy_n(mu_num,row,num2),field2d%epgy_s(mu_num,row,num2), &
  field2d%epgy_sw(mu_num,row,row),field2d%epgy_se(mu_num,row,row),field2d%epgy_nw(mu_num,row,row), &
  field2d%epgy_ne(mu_num,row,row),stat=ierr)

  field2d%epgyro=0.0
  field2d%epgy_w=0.0
  field2d%epgy_e=0.0
  field2d%epgy_n=0.0
  field2d%epgy_s=0.0
  field2d%epgy_sw=0.0
  field2d%epgy_se=0.0
  field2d%epgy_nw=0.0
  field2d%epgy_ne=0.0

   allocate(field2d%epgy_weight(mu_num,num1,num2),field2d%epgywg_w(mu_num,num1,row), &
  field2d%epgywg_e(mu_num,num1,row),field2d%epgywg_n(mu_num,row,num2),field2d%epgywg_s(mu_num,row,num2), &
  field2d%epgywg_sw(mu_num,row,row),field2d%epgywg_se(mu_num,row,row),field2d%epgywg_nw(mu_num,row,row), &
  field2d%epgywg_ne(mu_num,row,row),stat=ierr)

  field2d%epgy_weight=0.0
  field2d%epgywg_w=0.0
  field2d%epgywg_e=0.0
  field2d%epgywg_n=0.0
  field2d%epgywg_s=0.0
  field2d%epgywg_sw=0.0
  field2d%epgywg_se=0.0
  field2d%epgywg_nw=0.0
  field2d%epgywg_ne=0.0


     allocate(field2d%deng(mu_num,num1,num2),field2d%deng_w(mu_num,num1,row),&
  field2d%deng_e(mu_num,num1,row),field2d%deng_n(mu_num,row,num2),field2d%deng_s(mu_num,row,num2), &
  field2d%deng_sw(mu_num,row,row),field2d%deng_se(mu_num,row,row),field2d%deng_nw(mu_num,row,row), &
  field2d%deng_ne(mu_num,row,row), stat=ierr)

  field2d%deng=0.0
  field2d%deng_w=0.0
  field2d%deng_e=0.0
  field2d%deng_n=0.0
  field2d%deng_s=0.0
  field2d%deng_sw=0.0
  field2d%deng_se=0.0
  field2d%deng_nw=0.0
  field2d%deng_ne=0.0

     allocate(field2d%dengeq(mu_num,num1,num2),field2d%dengeq_w(mu_num,num1,row), &
  field2d%dengeq_e(mu_num,num1,row),field2d%dengeq_n(mu_num,row,num2),field2d%dengeq_s(mu_num,row,num2), &
  field2d%dengeq_sw(mu_num,row,row),field2d%dengeq_se(mu_num,row,row),field2d%dengeq_nw(mu_num,row,row), &
  field2d%dengeq_ne(mu_num,row,row), stat=ierr)

     field2d%dengeq=0.0
     field2d%dengeq_w=0.0
     field2d%dengeq_e=0.0
     field2d%dengeq_n=0.0
  field2d%dengeq_s=0.0
  field2d%dengeq_sw=0.0
  field2d%dengeq_se=0.0
  field2d%dengeq_nw=0.0 
  field2d%dengeq_ne=0.0

     allocate(field2d%dengeq_weight(mu_num,num1,num2),field2d%dengeqwg_w(mu_num,num1,row), &
  field2d%dengeqwg_e(mu_num,num1,row),field2d%dengeqwg_n(mu_num,row,num2),field2d%dengeqwg_s(mu_num,row,num2), &
  field2d%dengeqwg_sw(mu_num,row,row),field2d%dengeqwg_se(mu_num,row,row),field2d%dengeqwg_nw(mu_num,row,row), &
  field2d%dengeqwg_ne(mu_num,row,row), stat=ierr)
     
  field2d%dengeq_weight=0.0
  field2d%dengeqwg_w=0.0
  field2d%dengeqwg_e=0.0
  field2d%dengeqwg_n=0.0
  field2d%dengeqwg_s=0.0
  field2d%dengeqwg_sw=0.0
  field2d%dengeqwg_se=0.0
  field2d%dengeqwg_nw=0.0
  field2d%dengeqwg_ne=0.0 

     allocate(field2d%deng_weight(mu_num,num1,num2),field2d%dengwg_w(mu_num,num1,row), &
  field2d%dengwg_e(mu_num,num1,row),field2d%dengwg_n(mu_num,row,num2),field2d%dengwg_s(mu_num,row,num2), &
  field2d%dengwg_sw(mu_num,row,row),field2d%dengwg_se(mu_num,row,row),field2d%dengwg_nw(mu_num,row,row), &
  field2d%dengwg_ne(mu_num,row,row), stat=ierr)

  field2d%deng_weight=0.0
  field2d%dengwg_w=0.0
  field2d%dengwg_e=0.0
  field2d%dengwg_n=0.0
  field2d%dengwg_s=0.0
  field2d%dengwg_sw=0.0
  field2d%dengwg_se=0.0
  field2d%dengwg_nw=0.0
  field2d%dengwg_ne=0.0

     allocate(field2d%dengtot(num1,num2),field2d%dengtot_w(num1,row), &
  field2d%dengtot_e(num1,row),field2d%dengtot_n(row,num2),field2d%dengtot_s(row,num2), &
  field2d%dengtot_sw(row,row),field2d%dengtot_se(row,row),field2d%dengtot_nw(row,row), &
  field2d%dengtot_ne(row,row), stat=ierr)

  field2d%dengtot=0.0
  field2d%dengtot_w=0.0
  field2d%dengtot_e=0.0
  field2d%dengtot_n=0.0
  field2d%dengtot_s=0.0
  field2d%dengtot_sw=0.0
  field2d%dengtot_se=0.0
  field2d%dengtot_nw=0.0
  field2d%dengtot_ne=0.0

     allocate(field2d%dengeqtot(num1,num2), stat=ierr)
  
  field2d%dengeqtot=0.0


     allocate(field2d%epgysq_weight(num1,num2),field2d%epgysqwg_w(num1,row), &
  field2d%epgysqwg_e(num1,row),field2d%epgysqwg_n(row,num2),field2d%epgysqwg_s(row,num2), &
  field2d%epgysqwg_sw(row,row),field2d%epgysqwg_se(row,row),field2d%epgysqwg_nw(row,row), &
  field2d%epgysqwg_ne(row,row), stat=ierr)

  field2d%epgysq_weight=0.0
  field2d%epgysqwg_w=0.0
  field2d%epgysqwg_e=0.0
  field2d%epgysqwg_n=0.0
  field2d%epgysqwg_s=0.0
  field2d%epgysqwg_sw=0.0
  field2d%epgysqwg_se=0.0
  field2d%epgysqwg_nw=0.0
  field2d%epgysqwg_ne=0.0

  end subroutine allocate_memory_to_field_2d_gy

  subroutine allocate_memory_to_magfield_2d(field2d,num1,num2,row)
    type(pic_field_2d_base), pointer,intent(inout) :: field2d
    int4,intent(in) :: num1,num2,row
    int4 :: ierr

   if(.not.associated(field2d)) then
      allocate(field2d,stat=ierr)
   end if

   allocate(field2d%Bf03(num1,num2),field2d%Bf03_w(num1,row),&
  field2d%Bf03_e(num1,row),field2d%Bf03_n(row,num2),field2d%Bf03_s(row,num2),&
  field2d%Bf03_sw(row,row),field2d%Bf03_se(row,row),field2d%Bf03_nw(row,row),field2d%Bf03_ne(row,row),stat=ierr)

   allocate(field2d%Bf03wg(num1,num2),field2d%Bf03wg_w(num1,row),&
  field2d%Bf03wg_e(num1,row),field2d%Bf03wg_n(row,num2),field2d%Bf03wg_s(row,num2),&
  field2d%Bf03wg_sw(row,row),field2d%Bf03wg_se(row,row),field2d%Bf03wg_nw(row,row),field2d%Bf03wg_ne(row,row),stat=ierr)

  end subroutine

   subroutine initialize_parameters_2d(pic2d,pamearray)
      class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
      class(parameters_array_2d), pointer, intent(inout) :: pamearray
      int4 :: i,j,startind(2),size1,rank,k,size
      real8 :: dmu,delta(2)
      

      size=pic2d%layout2d%collective%size
       rank=pic2d%layout2d%collective%rank

       do i=1,2
         delta(i)=(pic2d%para2d%gxmax(i)-pic2d%para2d%gxmin(i))/  &
                  real(pic2d%para2d%numproc(i)*pic2d%para2d%cell_per_unit(i),8)
       end do

!       allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))

       do i=0,pic2d%para2d%numproc(2)-1
         do j=0,pic2d%para2d%numproc(1)-1
           size1=i*pic2d%para2d%numproc(1)+j
           if(size1==0) then
              pic2d%para2d%gboxmin(0,:)=pic2d%para2d%gxmin(:)
           else
             startind=startind_of_process(size1,pic2d%para2d%numproc,pic2d%layout2d)
             if(rank==0) then
               print*, "size1=",size1,"startind=",startind
             end if
             do k=1,2
               pic2d%para2d%gboxmin(size1,k)=pic2d%para2d%gxmin(k)+(startind(k)-1)*delta(k)
             end do
           end if

           pic2d%para2d%gboxmax(size1,1)=pic2d%para2d%gboxmin(size1,1)+delta(1)*real(pic2d%layout2d%boxes(size1)%i_max &
                         -pic2d%layout2d%boxes(size1)%i_min,8)
           pic2d%para2d%gboxmax(size1,2)=pic2d%para2d%gboxmin(size1,2)+delta(2)*real(pic2d%layout2d%boxes(size1)%j_max &
                        -pic2d%layout2d%boxes(size1)%j_min,8)

         end do
       end do

       pic2d%para2d%m_x1=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%i_max-pic2d%layout2d%boxes(rank)%i_min+1, &
       pic2d%para2d%gboxmin(rank,1),pic2d%para2d%gboxmax(rank,1),delta(1))
       pic2d%para2d%m_x2=>init_para_cartesian_mesh_1d(pic2d%layout2d%boxes(rank)%j_max-pic2d%layout2d%boxes(rank)%j_min+1, &
       pic2d%para2d%gboxmin(rank,2),pic2d%para2d%gboxmax(rank,2),delta(2))

!       dmu=pic2d%para2d%mumax/real(pic2d%para2d%mu_num,8)       
!       do i=1,pic2d%para2d%mu_num
!         pamearray%mu_nodes(i)=real(i-1,8)*dmu
!       end do

!       call muarray_euler_maclaurin_choice(pic2d%para2d%mumax,pic2d%para2d%mu_num, &
!                                           pamearray%mu_nodes,pamearray%mu_weights,pic2d%para2d%mu_scheme)

   end subroutine initialize_parameters_2d


   subroutine initialize_parameters_array_2d(mumax,mu_num,mu_scheme,pamearray,mu_nodes,mu_weights,munum_partition,temp)
     class(parameters_array_2d), pointer :: pamearray
     real8,intent(in) :: mumax,temp
     int4,intent(in)  :: mu_num,mu_scheme
     real8, dimension(:), pointer, intent(in) :: mu_nodes,mu_weights
     int4, dimension(:), pointer,intent(in) :: munum_partition
 
     int4 :: i,j,numdim     
!     allocate(mu_nodes(mu_num+1),mu_weights(mu_num+1))
!     call muarray_euler_maclaurin_choice(mumax,mu_num+1,  &
!                                         mu_nodes,mu_weights,mu_scheme)
     numdim=size(pamearray%temp_i,1)

     pamearray%mu_nodes(1:mu_num)=mu_nodes(1:mu_num)
     pamearray%mu_weights(1:mu_num)=mu_weights(1:mu_num)
     pamearray%munum_partition(1:mu_num)=munum_partition(1:mu_num)
 !    deallocate(mu_nodes,mu_weights) 

     do i=1,numdim
       do j=1, numdim
          pamearray%temp_i(i,j)=temp
          pamearray%temp_e(i,j)=temp
       enddo
     enddo
   end subroutine  initialize_parameters_array_2d

   subroutine computing_mu_number(mu_nods,mu_weigs,munum_part,mu_num,pic2d)
     class(pic_para_total2d_base), pointer :: pic2d
     real8, dimension(:), pointer,intent(inout) :: mu_nods,mu_weigs
     int4, dimension(:), pointer, intent(inout) :: munum_part
     int4, intent(inout) :: mu_num
     int4 :: mutest
     real8, dimension(:), pointer :: mu_nodes,mu_weights   
     int4, dimension(:), pointer :: munum_partition 
     real8 :: integ
     int4 :: i    
 
     mutest=pic2d%para2d%mulast
     allocate(mu_nodes(mutest),mu_weights(mutest),munum_partition(mutest))
     munum_partition=0
     do while(munum_partition(pic2d%para2d%mulast).le.pic2d%para2d%mu_tail)
       deallocate(mu_nodes,mu_weights,munum_partition)
       mutest=mutest+5
  !     print*, mutest
       if(mutest.ge.100) then
         print*, "#error: the numer of munodes exceeds the given upbound."
         stop
       end if
       allocate(mu_nodes(mutest),mu_weights(mutest),munum_partition(mutest))
       call muarray_euler_maclaurin_choice(pic2d%para2d%mumax,mutest,  &
                                         mu_nodes,mu_weights,pic2d%para2d%mu_scheme)
       integ=0.0
       do i=1,mutest
         integ=integ+exp(-mu_nodes(i)/pic2d%para2d%tempt)*mu_weights(i)
       end do

       do i=1,mutest
          munum_partition(i)=NINT(real(pic2d%para2d%numparticles,8)*exp(-mu_nodes(i)/ &
             pic2d%para2d%tempt)*mu_weights(i)/integ)
       end do
     end do

     mu_num=pic2d%para2d%mulast
     do while(mu_num.le.mutest)
       if(munum_partition(mu_num+1).ge.1000) then
         mu_num=mu_num+1
       else
         exit
       end if
     end do
     
     mu_num=mu_num-1   !!! Because the first elements of mu_ndoes and mu_weights equal zero.
     mu_nods(1:mu_num)=mu_nodes(2:mu_num+1)
     mu_weigs(1:mu_num)=mu_weights(2:mu_num+1)
     munum_part(1:mu_num)=munum_partition(2:mu_num+1)
!     pic2d%para2d%mu_num=mu_num

     deallocate(mu_nodes,mu_weights,munum_partition)
  end subroutine computing_mu_number

end module paradata_layout
