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
          allocate_memory_to_field_2d, &
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
    allocate(pic2d%para2d%temp_i(pic2d%layout2d%global_sz1))
    allocate(pic2d%para2d%temp_e(pic2d%layout2d%global_sz2))
    allocate(pic2d%para2d%mu_nodes(pic2d%para2d%mu_num))
    allocate(pic2d%para2d%mu_weights(pic2d%para2d%mu_num))

    do i=1,size
      allocate(pic2d%ful2dsend_head(i-1)%ptr)
      allocate(pic2d%gy2dsend_head(i-1)%ptr)
    end do

  end function


  subroutine allocate_memory_to_field_2d(field2d,num1,num2,row)
    type(pic_field_2d_base), pointer,intent(inout) :: field2d
    int4,intent(in) :: num1,num2,row
    int4 :: ierr

   if(.not.associated(field2d)) then
      allocate(field2d,stat=ierr)
   end if

   allocate(field2d%ep(num1,num2),field2d%ep_w(num1,row),&
  field2d%ep_e(num1,row),field2d%ep_n(row,num2),field2d%ep_s(row,num2),&
  field2d%ep_sw(row,row),field2d%ep_se(row,row),field2d%ep_nw(row,row),field2d%ep_ne(row,row),stat=ierr)


     allocate(field2d%ep_weight(num1,num2),field2d%epwg_w(num1,row),&
  field2d%epwg_e(num1,row),field2d%epwg_n(row,num2),field2d%epwg_s(row,num2),&
  field2d%epwg_sw(row,row),field2d%epwg_se(row,row),field2d%epwg_nw(row,row),&
  field2d%epwg_ne(row,row), stat=ierr)

   allocate(field2d%Bf03(num1,num2),field2d%Bf03_w(num1,row),&
  field2d%Bf03_e(num1,row),field2d%Bf03_n(row,num2),field2d%Bf03_s(row,num2),&
  field2d%Bf03_sw(row,row),field2d%Bf03_se(row,row),field2d%Bf03_nw(row,row),field2d%Bf03_ne(row,row),stat=ierr)

   allocate(field2d%Bf03wg(num1,num2),field2d%Bf03wg_w(num1,row),&
  field2d%Bf03wg_e(num1,row),field2d%Bf03wg_n(row,num2),field2d%Bf03wg_s(row,num2),&
  field2d%Bf03wg_sw(row,row),field2d%Bf03wg_se(row,row),field2d%Bf03wg_nw(row,row),field2d%Bf03wg_ne(row,row),stat=ierr)

     allocate(field2d%den(num1,num2),field2d%den_w(num1,row),&
  field2d%den_e(num1,row),field2d%den_n(row,num2),field2d%den_s(row,num2), &
  field2d%den_sw(row,row),field2d%den_se(row,row),field2d%den_nw(row,row), &
  field2d%den_ne(row,row), stat=ierr)

     allocate(field2d%denequ(num1,num2),field2d%denequ_w(num1,row), &
  field2d%denequ_e(num1,row),field2d%denequ_n(row,num2),field2d%denequ_s(row,num2), &
  field2d%denequ_sw(row,row),field2d%denequ_se(row,row),field2d%denequ_nw(row,row), &
  field2d%denequ_ne(row,row), stat=ierr)

     allocate(field2d%epgyro(num1,num2),field2d%epgy_w(num1,row), &
  field2d%epgy_e(num1,row),field2d%epgy_n(row,num2),field2d%epgy_s(row,num2), &
  field2d%epgy_sw(row,row),field2d%epgy_se(row,row),field2d%epgy_nw(row,row), &
  field2d%epgy_ne(row,row), stat=ierr)

     allocate(field2d%epgy_weight(num1,num2),field2d%epgywg_w(num1,row), &
  field2d%epgywg_e(num1,row),field2d%epgywg_n(row,num2),field2d%epgywg_s(row,num2), &
  field2d%epgywg_sw(row,row),field2d%epgywg_se(row,row),field2d%epgywg_nw(row,row), &
  field2d%epgywg_ne(row,row), stat=ierr)

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
