module paradata_layout
#include "work_precision.h"
use cartesian_mesh, only: init_para_cartesian_mesh_1d
use utilities_module, only: muarray_euler_maclaurin_choice
use paradata_type, only: pic_para_2d_base, &
                         pic_para_total2d_base
use piclayout, only: pic_field_2d_base
use paradata_utilities, only: startind_of_process                         
implicit none

public::  initialize_pic_para_2d_base, &
          initialize_pic_para_total2d_base, &
          initialize_pic_para_total2d_base_2nd, &
          allocate_memory_to_field_2d_ful, &
          allocate_memory_to_field_2d_gy,  &
          allocate_memory_to_magfield_2d,  &
          initialize_parameters_2d
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
    allocate(pic2d%ful2dsend_head(0:size-1))
    allocate(pic2d%gy2d_head)
    allocate(pic2d%gy2dsend_head(0:size-1))

    allocate(pic2d%para2d)
    allocate(pic2d%para2d%gboxmin(0:size-1,2),pic2d%para2d%gboxmax(0:size-1,2))

    do i=1,size
      allocate(pic2d%ful2dsend_head(i-1)%ptr)
      allocate(pic2d%gy2dsend_head(i-1)%ptr)
    end do

  end function initialize_pic_para_total2d_base

  subroutine initialize_pic_para_total2d_base_2nd(pic2d)
     class(pic_para_total2d_base), pointer, intent(inout) :: pic2d

     allocate(pic2d%para2d%temp_i(pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2))
     allocate(pic2d%para2d%temp_e(pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2))
     allocate(pic2d%para2d%mu_nodes(pic2d%para2d%mu_num))
     allocate(pic2d%para2d%mu_weights(pic2d%para2d%mu_num))

  end subroutine initialize_pic_para_total2d_base_2nd


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

   allocate(field2d%ep_weight(num1,num2),field2d%epwg_w(num1,row),&
  field2d%epwg_e(num1,row),field2d%epwg_n(row,num2),field2d%epwg_s(row,num2),&
  field2d%epwg_sw(row,row),field2d%epwg_se(row,row),field2d%epwg_nw(row,row),&
  field2d%epwg_ne(row,row), stat=ierr)

     allocate(field2d%denf(num1,num2),field2d%denf_w(num1,row),&
  field2d%denf_e(num1,row),field2d%denf_n(row,num2),field2d%denf_s(row,num2), &
  field2d%denf_sw(row,row),field2d%denf_se(row,row),field2d%denf_nw(row,row), &
  field2d%denf_ne(row,row), stat=ierr)

     allocate(field2d%denfeq(num1,num2),field2d%denfeq_w(num1,row), &
  field2d%denfeq_e(num1,row),field2d%denfeq_n(row,num2),field2d%denfeq_s(row,num2), &
  field2d%denfeq_sw(row,row),field2d%denfeq_se(row,row),field2d%denfeq_nw(row,row), &
  field2d%denfeq_ne(row,row), stat=ierr)


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
  field2d%ep_ne(row,row),stat=ierr)

   allocate(field2d%gep_weight(num1,num2),field2d%gepwg_w(num1,row), &
  field2d%gepwg_e(num1,row),field2d%gepwg_n(row,num2),field2d%gepwg_s(row,num2),&
  field2d%gepwg_sw(row,row),field2d%gepwg_se(row,row),field2d%gepwg_nw(row,row),&
  field2d%gepwg_ne(row,row), stat=ierr)

   allocate(field2d%epgyro(mu_num,num1,num2),field2d%epgy_w(mu_num,num1,row), &
  field2d%epgy_e(mu_num,num1,row),field2d%epgy_n(mu_num,row,num2),field2d%epgy_s(mu_num,row,num2), &
  field2d%epgy_sw(mu_num,row,row),field2d%epgy_se(mu_num,row,row),field2d%epgy_nw(mu_num,row,row), &
  field2d%epgy_ne(mu_num,row,row),stat=ierr)

   allocate(field2d%epgy_weight(mu_num,num1,num2),field2d%epgywg_w(mu_num,num1,row), &
  field2d%epgywg_e(mu_num,num1,row),field2d%epgywg_n(mu_num,row,num2),field2d%epgywg_s(mu_num,row,num2), &
  field2d%epgywg_sw(mu_num,row,row),field2d%epgywg_se(mu_num,row,row),field2d%epgywg_nw(mu_num,row,row), &
  field2d%epgywg_ne(mu_num,row,row),stat=ierr)

     allocate(field2d%deng(mu_num,num1,num2),field2d%deng_w(mu_num,num1,row),&
  field2d%deng_e(mu_num,num1,row),field2d%deng_n(mu_num,row,num2),field2d%deng_s(mu_num,row,num2), &
  field2d%deng_sw(mu_num,row,row),field2d%deng_se(mu_num,row,row),field2d%deng_nw(mu_num,row,row), &
  field2d%deng_ne(mu_num,row,row), stat=ierr)

     allocate(field2d%dengeq(mu_num,num1,num2),field2d%dengeq_w(mu_num,num1,row), &
  field2d%dengeq_e(mu_num,num1,row),field2d%dengeq_n(mu_num,row,num2),field2d%dengeq_s(mu_num,row,num2), &
  field2d%dengeq_sw(mu_num,row,row),field2d%dengeq_se(mu_num,row,row),field2d%dengeq_nw(mu_num,row,row), &
  field2d%dengeq_ne(mu_num,row,row), stat=ierr)

     allocate(field2d%dengeq_weight(mu_num,num1,num2),field2d%dengeqwg_w(mu_num,num1,row), &
  field2d%dengeqwg_e(mu_num,num1,row),field2d%dengeqwg_n(mu_num,row,num2),field2d%dengeqwg_s(mu_num,row,num2), &
  field2d%dengeqwg_sw(mu_num,row,row),field2d%dengeqwg_se(mu_num,row,row),field2d%dengeqwg_nw(mu_num,row,row), &
  field2d%dengeqwg_ne(mu_num,row,row), stat=ierr)

     allocate(field2d%deng_weight(mu_num,num1,num2),field2d%dengwg_w(mu_num,num1,row), &
  field2d%dengwg_e(mu_num,num1,row),field2d%dengwg_n(mu_num,row,num2),field2d%dengwg_s(mu_num,row,num2), &
  field2d%dengwg_sw(mu_num,row,row),field2d%dengwg_se(mu_num,row,row),field2d%dengwg_nw(mu_num,row,row), &
  field2d%dengwg_ne(mu_num,row,row), stat=ierr)



     allocate(field2d%dengtot(num1,num2), stat=ierr)
     allocate(field2d%dengeqtot(num1,num2), stat=ierr)



     allocate(field2d%epgyro(mu_num,num1,num2),field2d%epgy_w(mu_num,num1,row), &
  field2d%epgy_e(mu_num,num1,row),field2d%epgy_n(mu_num,row,num2),field2d%epgy_s(mu_num,row,num2), &
  field2d%epgy_sw(mu_num,row,row),field2d%epgy_se(mu_num,row,row),field2d%epgy_nw(mu_num,row,row), &
  field2d%epgy_ne(mu_num,row,row), stat=ierr)

     allocate(field2d%epgy_weight(mu_num,num1,num2),field2d%epgywg_w(mu_num,num1,row), &
  field2d%epgywg_e(mu_num,num1,row),field2d%epgywg_n(mu_num,row,num2),field2d%epgywg_s(mu_num,row,num2), &
  field2d%epgywg_sw(mu_num,row,row),field2d%epgywg_se(mu_num,row,row),field2d%epgywg_nw(mu_num,row,row), &
  field2d%epgywg_ne(mu_num,row,row), stat=ierr)

     allocate(field2d%epgysq_weight(num1,num2),field2d%epgysqwg_w(num1,row), &
  field2d%epgysqwg_e(num1,row),field2d%epgysqwg_n(row,num2),field2d%epgysqwg_s(row,num2), &
  field2d%epgysqwg_sw(row,row),field2d%epgysqwg_se(row,row),field2d%epgysqwg_nw(row,row), &
  field2d%epgysqwg_ne(row,row), stat=ierr)



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

   subroutine initialize_parameters_2d(pic2d)
      class(pic_para_total2d_base), pointer, intent(inout) :: pic2d
      int4 :: i,j,startind(2),size1,rank,k
      real8 :: dmu,delta(2)

       rank=pic2d%layout2d%collective%rank

       do i=1,2
         delta(i)=1._f64/real(pic2d%para2d%cell_per_unit(i),8)
       end do

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

       call muarray_euler_maclaurin_choice(pic2d%para2d%mumax,pic2d%para2d%mu_num, &
                                           pic2d%para2d%mu_nodes,pic2d%para2d%mu_weights,pic2d%para2d%mu_scheme)


   end subroutine initialize_parameters_2d


end module paradata_layout
