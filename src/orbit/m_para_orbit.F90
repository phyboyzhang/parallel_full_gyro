module m_para_orbit
#include "work_precision.h"
  use m_para_spline, only: para_spl2d_firstorder_derivatve_point_per_per, &
                            para_compute_spl2d_field_point_per_per
                          
  
  use piclayout, only: ful2d_node,gy2d_node, gy2dmu_node, &
                       ful2dsend_node,gy2dsend_node
  use paradata_type, only: pic_para_total2d_base
  use orbit_data_base, only: pointin_node,pointin_gy_node,fuloutnode,rk4ful2dnode, &
                             gyoutnode,rk4gy2dnode 
  use paradata_utilities, only: prepare_sendrecv_for_mpialltoallv, &
                                prepare_sendrecv_for_mpialltoallv_simple, &
                                compute_process_of_point_per_per, &
                                prep_for_mpi_alltoallv_with_zeroinput, &
                                coordinates_pointoutbound_per_per  
  use m_picutilities, only: mpi2d_alltoallv_send_particle_2d, &
                            sort_particles_among_ranks, &
                            sort_particles_among_ranks_gyallmu, &
                            mpi2d_alltoallv_send_particle_2d_gy 
 implicit none
 include "mpif.h"

  public :: borissolve, &
            fulrk4solve, &
            gyrork4solve, &
            gyrork4solveallmu, &
            boris_single, &
            para_obtain_interpolation_elefield_per_per_ful, &
            sortposition_by_process_ful2d, &
            compute_f_of_points_out_orbit_ful2d, &
            compute_f_of_points_in_orbit_ful2d, &
            fulrkfunc_f_per_per, &
            fulrk4solve_and_sort, &
            borissolve_and_sort,  &
            gyrork4solveallmu_and_sort            

contains


subroutine boris_single(x,v,dtful,magf,elef)
    real8, dimension(:), intent(inout) :: x,v
    real8, intent(in) :: dtful
    real8, dimension(3), intent(in) :: magf, elef
    real8 :: pomat(3,3),nomat(3,3),xco(3),vel(3),vec(3,1)
    real8 :: x1(2)   
    int4 :: N=3, NRHS=1, LDA=3, LDB=3,INFO
    int4 :: IPIV(3)
    int4 i,j

   pomat(1,1)=0.0_f64
   pomat(1,2)=-magf(3)*dtful/2.0_f64
   pomat(1,3)=magf(2)*dtful/2.0_f64
   pomat(2,1)=magf(3)*dtful/2.0_f64
   pomat(2,2)=0.0_f64
   pomat(2,3)=-magf(1)*dtful/2.0_f64
   pomat(3,1)=-magf(2)*dtful/2.0_f64
   pomat(3,2)=magf(1)*dtful/2.0_f64
   pomat(3,3)=0.0_f64
   
   nomat(:,:)=-pomat(:,:)
   nomat(1,1)=nomat(1,1)+1.0_f64
   nomat(2,2)=nomat(2,2)+1.0_f64
   nomat(3,3)=nomat(3,3)+1.0_f64

   pomat(1,1)=pomat(1,1)+1.0_f64
   pomat(2,2)=pomat(2,2)+1.0_f64
   pomat(3,3)=pomat(3,3)+1.0_f64

   do i=1,3,1
    xco(i)=x(i)
    vel(i)=v(i)
   end do

   do i=1,3,1
     vec(i,1)=0.0_f64
     do j=1,3,1
     vec(i,1)=vec(i,1)+nomat(i,j)*vel(j)
     end do
     vec(i,1)=vec(i,1)+dtful*elef(i)
   end do

   call dgesv(N,NRHS,pomat,LDA,IPIV,vec,LDB,INFO)

   do i=1,3,1
     xco(i)=xco(i)+vec(i,1)*dtful
   end do

   do i=1,3,1
    x(i)=xco(i)
    v(i)=vec(i,1)
   end do
   return
 end subroutine boris_single

 subroutine borissolve(ful2d_head,pic2d,numleft)
   class(pic_para_total2d_base), pointer :: pic2d
   class(ful2d_node), pointer, intent(inout) :: ful2d_head
   int4,optional, intent(in) :: numleft
   class(ful2d_node), pointer :: tmp
   character(25) :: geometry,boundary
   real8 :: x(3), v(3),elef(3),magf(3)  
   real8 :: dtful
   int4  :: row

   row=pic2d%para2d%row
   dtful=pic2d%para2d%dtful  
   geometry=pic2d%para2d%geometry
   boundary=pic2d%para2d%boundary
 
  if(present(numleft)) then
    if(numleft==0) then
       goto 20
    else 
       goto 30
    end if
  else

30   tmp=>ful2d_head
   select case (geometry)
   case("cartesian")
     select case (boundary)
       case("double_per")   
       do while(associated(tmp))
          if(.not.associated(tmp%next)) then
             exit
          else
             elef(3)=0.0_f64
             magf(1:2)=0.0_F64
             x(1:2)=tmp%coords(1:2)
             x(3)=0.0
             v(1:2)=tmp%coords(3:4)
             v(3)=0.0
             call para_obtain_interpolation_elefield_per_per_ful(tmp%coords(1:2), elef, pic2d)
             magf(3)=para_compute_spl2d_field_point_per_per(tmp%coords(1:2),pic2d%para2d%m_x1, &
                     pic2d%para2d%m_x2,row, pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,& 
                     pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                     pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                     pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw)   
             call boris_single(x,v,dtful,magf,elef)
             tmp%coords(1:2)=x(1:2)
             tmp%coords(3:4)=v(1:2)
             tmp=>tmp%next
         end if
      end do

      case("nat_per")
        print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
        stop
      case default
        stop
      end select

   case("polar")
     print*, "#ERROR: The current version doesn't include the polar geometry condition."
     stop

   case default
     stop
   end select

   nullify(tmp)
 
  end if

20  end subroutine borissolve


 subroutine borissolve_and_sort(ful2d_head,pic2d)
   class(pic_para_total2d_base), pointer :: pic2d
   class(ful2d_node), pointer, intent(inout) :: ful2d_head
   class(ful2dsend_node), dimension(:), pointer :: ful2dsend_head
   int4, dimension(:), pointer :: num
   int4 :: csize,i

   csize=pic2d%layout2d%collective%size
   allocate(ful2dsend_head(0:csize-1), num(0:csize-1))
   do i=0,csize-1
     allocate(ful2dsend_head(i)%ptr)
   enddo
   num=0
   call borissolve(ful2d_head,pic2d)
   call sort_particles_among_ranks(ful2d_head,ful2dsend_head, pic2d, num)

   call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d)

   deallocate(ful2dsend_head,num)
 end subroutine borissolve_and_sort

 subroutine fulrkfunc_f_per_per(f,vec,elef,magf)
   ! It includes the magnetic field perturbation and the electric field perturbation
   real8, dimension(:), intent(in) :: vec,elef,magf
   real8, dimension(:), intent(inout) :: f
   real8 :: x(3),v(3)

   x(1)=vec(1)
   x(2)=vec(2)
   x(3)=vec(3)
   v(1)=vec(4)
   v(2)=vec(5)
   v(3)=vec(6)

   f(1)=vec(4)
   f(2)=vec(5)
   f(3)=vec(6)
   f(4)=elef(1)+v(2)*magf(3)-v(3)*magf(2)
   f(5)=elef(2)+v(3)*magf(1)-v(1)*magf(3)
   f(6)=elef(3)+v(1)*magf(2)-v(2)*magf(1)

!  print*, "elef=",elef(1:2), "magf=",magf(1:2)
!   print*, "f=",f(1:2),f(4:5)
  return
 end subroutine fulrkfunc_f_per_per

 subroutine fulrk4solve(ful2d_head,pic2d,rk4order,iter_num)
  ! It includes the magnetic field perturbation and the electric field perturbation
   class(pic_para_total2d_base), pointer :: pic2d
   class(ful2d_node), pointer, intent(inout) :: ful2d_head
   int4,optional,intent(in) :: iter_num
   class(ful2d_node), pointer  :: tmp
   class(rk4ful2dnode), pointer :: partlist,partmp
   class(pointin_node), pointer :: inlist, intmp  
   int4, intent(in) :: rk4order        
   int4 :: size,numgr,rank
   int4, dimension(:), pointer :: num,recnum
   real8 :: dt,vec0(6)
   int4 :: i,j,order,h,comm
  
    size=pic2d%layout2d%collective%size
    comm=pic2d%layout2d%collective%comm
    dt=pic2d%para2d%dtful
    allocate(num(0:size-1),recnum(0:size-1))
    allocate(partlist,inlist)
    num=0
    recnum=0 
    rank=pic2d%layout2d%collective%rank
    tmp=>ful2d_head
    partmp=>partlist

    do while(associated(tmp)) 
      if(.not.associated(tmp%next)) then
        exit
      else
         partmp%coords(1:4)=tmp%coords(1:4)
         partmp%vec(1:2)=partmp%coords(1:2)
         partmp%vec(4:5)=partmp%coords(3:4)
         partmp%vec(3)=0.0
         partmp%vec(6)=0.0
 !        partmp%vec0(1:6)=partmp%vec(1:6)
         partmp%at=0 
         tmp=>tmp%next 
         allocate(partmp%next)
         partmp=>partmp%next    
      end if
    end do

    numgr=4

    if(rk4order==2) then
       order=1
       num=0
       recnum=0
       call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)
!     call simple_f(partlist,order)

!if(iter_num==1.or.iter_num==7) then
!partmp=>partlist
!do while(associated(partmp))
!   if(.not.associated(partmp%next)) then
!      exit
!   else
!      print*, "rank1=",rank,"iter_num",iter_num,"coords1=",partmp%coords(1:2)
!      partmp=>partmp%next
!   end if
!end do
!
!end if

!if(iter_num==7) then
!print*, "rank1=",rank
!end if
!   if(iter_num==2) then
!      print*, 1
!   end if

       deallocate(inlist)
       nullify(inlist)
       allocate(inlist)
       recnum=0
       num=0
       order=2
       call sortposition_by_process_ful2d(num,recnum,partlist,inlist,pic2d,dt,order,rk4order,iter_num)
       call compute_f_of_points_out_orbit_ful2d(num,numgr,recnum,inlist,partlist,pic2d,order,iter_num)
       call compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)

       tmp=>ful2d_head
       partmp=>partlist
       do while(associated(partmp))
         if(.not.associated(partmp%next)) then
            exit
         else
           if(.not.associated(tmp)) then
             exit
           else
             vec0(1:2)=partmp%coords(1:2)
             vec0(4:5)=partmp%coords(3:4)
             partmp%vec(:)=vec0(:)+(dt/2.0_f64)*(partmp%f1(:)+partmp%f2(:))
             tmp%coords(1:2)=partmp%vec(1:2)
             tmp%coords(3:4)=partmp%vec(4:5)
!         tmp%coords(1:2)=tmp%coords(1:2)+0.5
!         tmp%coords(3:4)=tmp%coords(3:4)+0.1 
 
             tmp=>tmp%next
           end if
           partmp=>partmp%next
         end if
       end do



     else if(rk4order==4) then

        order=1
        num=0
        recnum=0
        call  compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)

        deallocate(inlist)
        nullify(inlist)
        allocate(inlist)
        recnum=0
        num=0
        order=2
        call sortposition_by_process_ful2d(num,recnum,partlist,inlist,pic2d,dt,order,rk4order,iter_num)
        call compute_f_of_points_out_orbit_ful2d(num,numgr,recnum,inlist,partlist,pic2d,order,iter_num)
        call compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)

        deallocate(inlist)
        nullify(inlist)
        allocate(inlist)
        recnum=0
        num=0
        order=3
        call sortposition_by_process_ful2d(num,recnum,partlist,inlist,pic2d,dt,order,rk4order,iter_num)
        call compute_f_of_points_out_orbit_ful2d(num,numgr,recnum,inlist,partlist,pic2d,order,iter_num)
        call compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)

        deallocate(inlist)
        nullify(inlist)
        allocate(inlist)
        recnum=0
        num=0
        order=4
        call sortposition_by_process_ful2d(num,recnum,partlist,inlist,pic2d,dt,order,rk4order,iter_num)
        call compute_f_of_points_out_orbit_ful2d(num,numgr,recnum,inlist,partlist,pic2d,order,iter_num)
        call compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)    

        tmp=>ful2d_head
        partmp=>partlist
        do while(associated(tmp).and.associated(partmp))
          if(.not.associated(partmp%next).or..not.associated(tmp%next)) then
            exit
          else
            vec0(1:2)=partmp%coords(1:2)
            vec0(4:5)=partmp%coords(3:4)
            partmp%vec(:)=vec0(:)+(dt/6.0_f64)*(partmp%f1(:)+2.0_f64*partmp%f2(:)+2.0_f64*partmp%f3(:) & 
                       +partmp%f4(:))
!            print*,"dvec=", (dt/6.0_f64)*(partmp%f1(:)+2.0_f64*partmp%f2(:)+2.0_f64*partmp%f3(:) &
!                       +partmp%f4(:)) 
!            print*, tmp%coords(1:2)
!            print*, "f1=",partmp%f1(1:2)
!            print*, "f2=",partmp%f2(1:2)
!            print*, "f3=",partmp%f3(1:2)
!            print*, "f4=",partmp%f4(1:2)
            tmp%coords(1:2)=partmp%vec(1:2)
            tmp%coords(3:4)=partmp%vec(4:5)
            tmp=>tmp%next
            partmp=>partmp%next 
          end if
        end do

     end if !!! end the rk4order choice

    deallocate(recnum,num,inlist,partlist)
    nullify(partlist)
    nullify(inlist)
    nullify(partmp)
    nullify(tmp)
  
  end subroutine fulrk4solve

  subroutine fulrk4solve_and_sort(ful2d_head,pic2d,rk4order,iter_num)
    class(pic_para_total2d_base), pointer :: pic2d
    class(ful2d_node), pointer, intent(inout) :: ful2d_head
    int4, intent(in) :: rk4order, iter_num
    class(ful2dsend_node), dimension(:), pointer :: ful2dsend_head
    int4, dimension(:), pointer :: num
    int4 :: csize, i

    csize=pic2d%layout2d%collective%size
    allocate(ful2dsend_head(0:csize-1), num(0:csize-1))
    do i=1,csize
      allocate(ful2dsend_head(i))
    enddo
    num=0  
    call fulrk4solve(ful2d_head,pic2d,rk4order,iter_num)

    call sort_particles_among_ranks(ful2d_head,ful2dsend_head, pic2d, num)
  
    call mpi2d_alltoallv_send_particle_2d(ful2d_head,ful2dsend_head,num,pic2d) 
   
    deallocate(ful2dsend_head,num)
  end subroutine fulrk4solve_and_sort

!  subroutine simple_f_for_test(partlist,order)
!     class(rk4ful2dnode), pointer,intent(inout) :: partlist
!     class(rk4ful2dnode), pointer :: partmp
!     int4,intent(in) :: order     
!     partmp=>partlist
!     do while(associated(partmp))
!       if(.not.associated(partmp%next)) then
!          exit
!       else
!          if(order==1) then
!            partmp%f1(1:2)= 0.1   ! partmp%vec(1:2)
!            partmp%f1(3:6)=0.0 
!          else if(order==2) then
!            partmp%f2(1:2)= 0.1    ! partmp%vec(1:2)
!            partmp%f2(3:6)=0.0
!          end if
!          partmp=>partmp%next
!       end if
!     end do
!     nullify(partmp)
!  end subroutine

  subroutine sortposition_by_process_ful2d(num,recnum,partlist,inlist,pic2d,dt,order,rk4order,iter_num)
     class(pic_para_total2d_base),intent(in), pointer :: pic2d
     class(rk4ful2dnode), pointer, intent(inout) :: partlist
     class(pointin_node), pointer,intent(inout) :: inlist
     int4, dimension(:), pointer, intent(inout) :: recnum, num
     int4, intent(in) :: order,rk4order,iter_num
     real8, intent(in) :: dt 
     class(fuloutnode), dimension(:),pointer :: outhead,outtmp
     class(rk4ful2dnode), pointer :: partmp    
     real8 :: coords(2),vec0(6)
     character(25) :: geometry,boundary
     int4  :: prank,rank,size,row,numgr,comm
     int4 :: i,h

     geometry=pic2d%para2d%geometry
     boundary=pic2d%para2d%boundary 
     rank=pic2d%layout2d%collective%rank
     comm=pic2d%layout2d%collective%comm
     size=pic2d%layout2d%collective%size
     allocate(outhead(0:size-1),outtmp(0:size-1))

     do i=0, size-1
       allocate(outhead(i)%ptr)
       outtmp(i)%ptr=>outhead(i)%ptr
     end do

     numgr=3
     row=pic2d%para2d%row

     select case (geometry)
       case("cartesian")
       select case (boundary)
         case("double_per")    
            partmp=>partlist
             num=0
             h=0
             do while(associated(partmp))
               if(.not.associated(partmp%next)) then
                  exit
               else 
                 h=h+1
                   vec0(1:2)=partmp%coords(1:2)
                   vec0(4:5)=partmp%coords(3:4)                                       
                  if(order==2) then
                     if(rk4order==2) then
                        partmp%vec(:)=vec0(:)+dt*partmp%f1(:)
                     else
                        partmp%vec(:)=vec0(:)+0.5_F64*dt*partmp%f1(:)
                     end if
                  else if(order==3) then
                     partmp%vec(:)=vec0(:)+0.5_F64*dt*partmp%f2(:)
                  else if(order==4) then
                     partmp%vec(:)=vec0(:)+dt*partmp%f3(:)
                  end if
!                  call coordinates_pointoutbound_per_per(partmp%vec(1:2),pic2d)
                  prank=compute_process_of_point_per_per(partmp%vec(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin,&
                        pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
                  if(prank.ne.rank) then
                     num(prank)=num(prank)+1
                     outtmp(prank)%ptr%coords=partmp%vec(1:2)
                     outtmp(prank)%ptr%numpoint=h
                     allocate(outtmp(prank)%ptr%next)
                     outtmp(prank)%ptr=>outtmp(prank)%ptr%next

                     partmp%at= 1
                  else
                     partmp%at= 0
                  end if
                  partmp=>partmp%next      
               end if 
            end do
            call mpi2d_alltoallv_send_points_orbit(num,numgr,recnum,outhead,inlist,pic2d,order)  
  
       case("nat_per")
          print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
          stop

       case default
          print*, "#ERROR: Input the correct boundary condition."
          stop
       end select

     case("polar")
        print*, "#ERROR: The current version doesn't include the polar geometry condition."
        stop

     case default
        print*, "#ERROR: Input the correct geometry condition."
        stop
     end select

     deallocate(outhead,outtmp)
     nullify(partmp)
    end subroutine sortposition_by_process_ful2d


    subroutine compute_f_of_points_in_orbit_ful2d(partlist,pic2d,order)
      class(pic_para_total2d_base), pointer,intent(in) :: pic2d
      class(rk4ful2dnode), pointer, intent(inout) :: partlist
      int4,intent(in) :: order
      class(rk4ful2dnode), pointer :: partmp
      character(25) :: geometry,boundary
      real8 :: magf(3),elef(3)
      int4 :: row, rank

      row=pic2d%para2d%row
      rank=pic2d%layout2d%collective%rank
      geometry=pic2d%para2d%geometry
      boundary=pic2d%para2d%boundary  
      partmp=>partlist

!      magf(1:2)=0._F64
      magf=0.0
      select case (geometry)
        case("cartesian")
        select case (boundary)
          case("double_per")   
          do while(associated(partmp)) 
            if(.not.associated(partmp%next)) then
               exit
            else
               if(partmp%at==0) then
                  call para_obtain_interpolation_elefield_per_per_ful(partmp%vec(1:2), elef(1:3), pic2d)
                  magf(3)=para_compute_spl2d_field_point_per_per(partmp%vec(1:2), &
                          pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, pic2d%field2d%bf03wg, &
                          pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                          pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                          pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw)           
                  if(order==1) then
                     call fulrkfunc_f_per_per(partmp%f1(1:6),partmp%vec(1:6),elef(1:3),magf(1:3))
                  else if(order==2) then
                    call fulrkfunc_f_per_per(partmp%f2(1:6),partmp%vec(1:6),elef(1:3),magf(1:3))
                  else if(order==3) then
                     call fulrkfunc_f_per_per(partmp%f3(1:6),partmp%vec(1:6),elef,magf)
                  else if(order==4) then
                     call fulrkfunc_f_per_per(partmp%f4(1:6),partmp%vec(1:6),elef,magf)
                  end if         
               end if
               partmp=>partmp%next   
            end if
         end do

          case("nat_per")
              print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
              stop
           case default
              stop
         end select

         case("polar")
            print*, "#ERROR: The current version doesn't include the polar geometry condition."
            stop

         case default
            stop
         end select
     
         nullify(partmp)
    end subroutine compute_f_of_points_in_orbit_ful2d

    subroutine compute_f_of_points_out_orbit_ful2d(num,numnode,recnum,inlist,partlist,pic2d,order,iter_num)
      int4, dimension(:), pointer,intent(in) :: num,recnum
      int4, intent(in) :: numnode
      int4, optional,intent(in) :: iter_num
      class(pic_para_total2d_base), pointer,intent(in) :: pic2d
      class(rk4ful2dnode), pointer, intent(inout) :: partlist
      class(pointin_node), pointer, intent(in) :: inlist
      int4, intent(in) :: order
      class(rk4ful2dnode), pointer :: partmp
      class(pointin_node), pointer :: intmp
      real8 :: magf,elef(3),mag(3)
      int4, dimension(:), pointer :: scounts,rcounts,sdispls,rdispls
      int4 :: numsend,numout
      real8, dimension(:), pointer :: sbuf2nd,rbuf2nd
      character(25) :: geometry,boundary
      int4  :: i,size,comm,ierr,h,row,rank

      geometry=pic2d%para2d%geometry
      boundary=pic2d%para2d%boundary  
      row=pic2d%para2d%row
      rank=pic2d%layout2d%collective%rank
      size=pic2d%layout2d%collective%size
      comm=pic2d%layout2d%collective%comm    
      allocate(scounts(0:size-1),rcounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))

      call prep_for_mpi_alltoallv_with_zeroinput(rank,size,numnode,recnum,num,scounts, &
           rcounts,sdispls,rdispls,numsend,numout)      
       
      if(numsend==0) then
         allocate(sbuf2nd(0:0))
       else
         allocate(sbuf2nd(0:numnode*numsend-1))
       end if
       if(numout==0) then
         allocate(rbuf2nd(0:0))
       else
         allocate(rbuf2nd(0:numnode*numout-1))
       end if

      intmp=>inlist
       select case (geometry)
          case("cartesian")
          select case (boundary)
            case("double_per")   
              if(numsend==0) then
              call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)        
              else 
              h=0 
              do while(associated(intmp)) 
                if(.not.associated(intmp%next)) then
                  exit
                else 
                  sbuf2nd(numnode*h)=intmp%numpoint
                  call para_obtain_interpolation_elefield_per_per_ful(intmp%x(1:2), elef, pic2d)
                  magf=para_compute_spl2d_field_point_per_per(intmp%x(1:2), &
                       pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
                       pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                       pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                       pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw) 
                  sbuf2nd(numnode*h+1:numnode*h+2)=elef(1:2)
                  sbuf2nd(numnode*h+3)=magf
                  h=h+1
                  intmp=>intmp%next
                end if
              enddo
              call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts, &
                rdispls,mpi_double_precision,comm,ierr) 
 
10           end if

            case("nat_per")
              print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
              stop
            case default
              stop
          end select

          case("polar")
            print*, "#ERROR: The current version doesn't include the polar geometry condition."
            stop
          case default
            stop
        end select
      

       mag(1)=0.0_f64
       mag(2)=0.0_f64
       elef(3)=0.0_f64
       do i=0,numout-1
         partmp=>partlist
         h=1 
         do while(h.ne.NINT(rbuf2nd(numnode*i))) 
           partmp=>partmp%next
           h=h+1
           if(.not.associated(partmp)) then
             print*, "#ERROR:The information of the current point is lost."
             stop
           end if
         end do 

         elef(1:2)=rbuf2nd(numnode*i+1:numnode*i+2)
         mag(3)=rbuf2nd(numnode*i+3)
         if(order==2) then
           call fulrkfunc_f_per_per(partmp%f2,partmp%vec,elef,mag) 
         else if(order==3) then
           call fulrkfunc_f_per_per(partmp%f3,partmp%vec,elef,mag)
         else if(order==4) then
           call fulrkfunc_f_per_per(partmp%f4,partmp%vec,elef,mag)
         end if
      end do


     deallocate(sbuf2nd,rbuf2nd)
     deallocate(scounts,rcounts,sdispls,rdispls) 
     nullify(partmp,intmp)
  end subroutine compute_f_of_points_out_orbit_ful2d

  subroutine gyrork_f_per_per(f,mu,elef,deri_Bf,magf)
   ! It includes the magnetic field perturbation and the electric field perturbation
   real8, dimension(:), intent(in) :: elef,deri_Bf,magf
   real8, dimension(:), intent(inout) :: f
   real8, intent(in) :: mu

   f(1)=(elef(2)-mu*deri_Bf(2))/magf(3)
   f(2)=(-elef(1)+mu*deri_Bf(1))/magf(3)     
 
   return
 end subroutine gyrork_f_per_per

 subroutine gyrork_f_per_per_2nd(f,mu,elef,deri_Bf,magf,deri_driftsquare)
   ! It includes the magnetic field perturbation and the electric field perturbation
   real8, dimension(:), intent(in) :: elef,deri_Bf,magf,deri_driftsquare
   real8, dimension(:), intent(inout) :: f
   real8, intent(in) :: mu

   f(1)=(elef(2)-mu*deri_Bf(2))/magf(3)-0.5_f64*deri_driftsquare(2)/magf(3)
   f(2)=(-elef(1)+mu*deri_Bf(1))/magf(3)+0.5_f64*deri_driftsquare(1)/magf(3)
   return
 end subroutine gyrork_f_per_per_2nd

 subroutine sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dt,order,muind,iter_num)
   class(pic_para_total2d_base),intent(in), pointer :: pic2d
   class(rk4gy2dnode), pointer,intent(in) :: partlist
   int4, dimension(:), pointer, intent(inout) :: recnum, num
   int4, intent(in) :: order,iter_num,muind
   real8, intent(in) :: dt 
   class(pointin_gy_node), pointer,intent(inout) :: inlist
   class(gyoutnode), dimension(:),pointer :: outtmp, outhead
   class(rk4gy2dnode), pointer :: partmp    
   class(pointin_node), pointer :: intmp
!   int4, dimension(:), pointer :: num
   real8 :: vec(3),magf(3),elef(3),deri_driftsquare(2),deri_bf(2)
   int4  :: prank,rank,size,row,numgr,comm
   character(25) :: geometry,boundary
   int4 :: i,h

    geometry=pic2d%para2d%geometry
    boundary=pic2d%para2d%boundary
    rank=pic2d%layout2d%collective%rank
    comm=pic2d%layout2d%collective%comm
    size=pic2d%layout2d%collective%size
    allocate(outhead(0:size-1),outtmp(0:size-1))
    do i=0, size-1
       allocate(outhead(i)%ptr)
       outtmp(i)%ptr=>outhead(i)%ptr
    end do

    numgr=4
    row=pic2d%para2d%row
   select case (geometry)
   case("cartesian")
   select case (boundary)
       case("double_per")
         num=0 
         h=0  
         partmp=>partlist  
    
          do while(associated(partmp)) 
            if(.not.associated(partmp%next)) then
              exit
            else
              h=h+1    
              if(order==2) then
               partmp%vec(1:2)=partmp%coords(1:2)+0.5_F64*dt*partmp%f1(:)
              else if(order==3) then
               partmp%vec(1:2)=partmp%coords(1:2)+0.5_F64*dt*partmp%f2(:)

!if(order==3.and.muind==3.and.rank==3) then
!print*, "h=",h,partmp%coords(1:2),partmp%f2,partmp%at
!endif

              else if(order==4) then
               partmp%vec(1:2)=partmp%coords(1:2)+dt*partmp%f3(:)
              end if 
              prank=compute_process_of_point_per_per(partmp%vec(1:2),pic2d%para2d%numproc,pic2d%para2d%gxmin, &
                       pic2d%para2d%gxmax, pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
              if(prank.ne.rank) then
                 num(prank)=num(prank)+1
                 outtmp(prank)%ptr%coords(1:2)=partmp%vec(1:2)
                 outtmp(prank)%ptr%coords(3)=partmp%coords(3)
                 outtmp(prank)%ptr%numpoint=h
                 allocate(outtmp(prank)%ptr%next)
                 outtmp(prank)%ptr=>outtmp(prank)%ptr%next

                 partmp%at=1
              else
                 num(prank)=num(prank)+1
                 partmp%at=0
              end if


             partmp=>partmp%next
           end if   
         end do 
!if(muind==3.and.order==2) then
!  print*, "rank=",rank, "num=",sum(num)
!endif
      call mpi2d_alltoallv_send_points_orbit_gy2d(num,numgr,recnum,outhead,inlist,pic2d) 
 
      case("nat_per")
        print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
        stop
      case default
        stop
      end select

   case("polar")
     print*, "#ERROR: The current version doesn't include the polar geometry condition."
     stop

   case default
     stop
   end select

    deallocate(outhead,outtmp)
    nullify(partmp)

   end subroutine sortposition_by_process_gyro2d

   subroutine compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,muind,order,gyroorder,iter_num)
     class(pic_para_total2d_base), pointer,intent(in) :: pic2d
     class(rk4gy2dnode), pointer, intent(inout) :: partlist
     int4, intent(in) :: muind,order,gyroorder,iter_num
     class(rk4gy2dnode), pointer :: partmp
     real8 :: magf(3),elef(3),deri_bf(2),deri_driftsquare(2)
     character(25) :: geometry,boundary
     int4 :: row,rank

     row=pic2d%para2d%row
     rank=pic2d%layout2d%collective%rank
     geometry=pic2d%para2d%geometry
     boundary=pic2d%para2d%boundary
     partmp=>partlist

     deri_bf=0.0
     select case (geometry)
     case("cartesian")
     select case (boundary)
       case("double_per")
     do while(associated(partmp)) 
       if(.not.associated(partmp%next)) then
         exit
       else
         if(partmp%at==0) then
           call para_obtain_interpolation_elefield_per_per_gyro(partmp%vec(1:2), muind,elef, pic2d)
           call compute_deri_of_magfield_per_per(deri_bf,partmp%vec(1:2),pic2d)         
           magf(3)=para_compute_spl2d_field_point_per_per(partmp%vec(1:2), &
                 pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
                 pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                 pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                 pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw)           
           
            select case (gyroorder)
              case (1)
                if(order==1) then
                 call gyrork_f_per_per(partmp%f1,partmp%coords(3),elef,deri_bf,magf) 
!if(iter_num==56.or.iter_num==57.or.iter_num==55) then
!print*, "vec=",partmp%vec(1:2),"#f1=",partmp%f1,"elef=",elef,"magf=",magf(3)
!endif
                else if(order==2) then 
                 call gyrork_f_per_per(partmp%f2,partmp%coords(3),elef,deri_bf,magf)
                else if(order==3) then
                 call gyrork_f_per_per(partmp%f3,partmp%coords(3),elef,deri_bf,magf)
                else if(order==4) then
                 call gyrork_f_per_per(partmp%f4,partmp%coords(3),elef,deri_bf,magf)
                end if         
              case (2)
                call compute_deri_of_sqgyep_per_per(deri_driftsquare,partmp%vec(1:2),pic2d)
                if(order==1) then
                 call gyrork_f_per_per_2nd(partmp%f1,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                else if(order==2) then
                 call gyrork_f_per_per_2nd(partmp%f2,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                else if(order==3) then
                 call gyrork_f_per_per_2nd(partmp%f3,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                else if(order==4) then
                 call gyrork_f_per_per_2nd(partmp%f4,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                end if

               case default
                stop
            end select

          end if
       end if
       partmp=>partmp%next
      end do


      case("nat_per")
        print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
        stop
      case default
        stop
      end select

   case("polar")
     print*, "#ERROR: The current version doesn't include the polar geometry condition."
     stop
   
   case default
     stop
   end select 

   nullify(partmp)

 end subroutine compute_f_of_points_in_orbit_gyro2d

 subroutine compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,muind,order,gyroorder,iter_num)
   int4, dimension(:), pointer, intent(in) :: num,recnum
!   int4, intent(in) :: numnode
   class(pic_para_total2d_base), pointer,intent(in) :: pic2d
   class(rk4gy2dnode), pointer, intent(inout) :: partlist
   class(pointin_gy_node), pointer, intent(in) :: inlist
   int4, intent(in) :: muind,order,gyroorder,iter_num
   class(rk4gy2dnode), pointer :: partmp
   class(pointin_gy_node), pointer :: intmp
   real8 :: mag,elef(3),magf(3),deri_driftsquare(2),deri_bf(2)
   int4, dimension(:), pointer :: scounts,rcounts,sdispls,rdispls
   int4 :: numsend,numout,numnode
   real8, dimension(:), pointer :: sbuf2nd,rbuf2nd
   character(25) :: geometry,boundary
   int4  :: i,dsize,comm,ierr,h,row,rank

   row=pic2d%para2d%row
   geometry=pic2d%para2d%geometry
   boundary=pic2d%para2d%boundary
   dsize=pic2d%layout2d%collective%size
   comm=pic2d%layout2d%collective%comm 
   rank=pic2d%layout2d%collective%rank
   if(gyroorder==1) then
      numnode=6
   else if(gyroorder==2) then
      numnode=8    
   end if
   allocate(scounts(0:dsize-1),rcounts(0:dsize-1),sdispls(0:dsize-1),rdispls(0:dsize-1))
   scounts=0
   rcounts=0
   sdispls=0
   rdispls=0
   call prep_for_mpi_alltoallv_with_zeroinput(rank,dsize,numnode,recnum,num,scounts, &
           rcounts,sdispls,rdispls,numsend,numout)
!if(rank==0) then
!  print*, "num=",num
!  print*, "recnum",recnum
!endif
   if(numsend==0) then
      allocate(sbuf2nd(0:0))
   else
      allocate(sbuf2nd(0:numnode*numsend-1))
   end if
   if(numout==0) then
      allocate(rbuf2nd(0:0))
   else
      allocate(rbuf2nd(0:numnode*numout-1))
   end if
 
   intmp=>inlist
      select case (geometry)
       case("cartesian")
         select case (boundary)
           case("double_per") 

            select case (gyroorder)
            case (1)
             if(numsend==0) then
               call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
             else
               h=0
               do while(associated(intmp))
                 if(.not.associated(intmp%next)) then
                   exit
                 else
                   sbuf2nd(numnode*h)=intmp%numpoint
                   call para_obtain_interpolation_elefield_per_per_gyro(intmp%x(1:2), muind,elef, pic2d)
                   call compute_deri_of_magfield_per_per(deri_bf,intmp%x(1:2),pic2d)
                   mag=para_compute_spl2d_field_point_per_per(intmp%x(1:2),pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
                       pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                       pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                       pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw) 
                   sbuf2nd(numnode*h+1:numnode*h+2)=elef(1:2)
                   sbuf2nd(numnode*h+3:numnode*h+4)=deri_bf(1:2)
                   sbuf2nd(numnode*h+5)=mag
                   h=h+1
                   intmp=>intmp%next
                 end if
               enddo
               call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
             end if

             magf(1)=0.0_f64
             magf(2)=0.0_f64
             elef(3)=0.0_f64
        
             do i=0,numout-1 
               partmp=>partlist
!if(muind==3.and.order==2.and.rank==3) then
!print*, "rank=",rank, "numout=",numout
!print*, NINT(rbuf2nd(numnode*i))
!endif  

   !!! Here, the logic of the  piece of code within the do while loop should be
   !paid attention 
               h=1
               do while(h.ne.NINT(rbuf2nd(numnode*i))) 
                 h=h+1

                 if(associated(partmp).and..not.associated(partmp%next)) then
                     print*, "rank=",rank,"h=",h, "rbuf=",rbuf2nd(numnode*i)
                     print*, "muind=",muind, "order=",order
                     print*, "#ERROR:The information of the current point is lost."
                     stop
                 else
                     partmp=>partmp%next 
                 end if
               end do 

               elef(1:2)=rbuf2nd(numnode*i+1:numnode*i+2)
               deri_bf(1:2)=rbuf2nd(numnode*i+3:numnode*i+4)
               magf(3)=rbuf2nd(numnode*i+5)

               if(order==2) then
                 call gyrork_f_per_per(partmp%f2,partmp%coords(3),elef,deri_bf,magf)

               else if(order==3) then
                 call gyrork_f_per_per(partmp%f3,partmp%coords(3),elef,deri_bf,magf)
               else if(order==4) then
                 call gyrork_f_per_per(partmp%f4,partmp%coords(3),elef,deri_bf,magf)
               end if     
        
!               nullify(partmp)
              end do
! if(muind==3.and.order==2.and.rank==3) then
!partmp=>partlist
!print*, "f2=",partmp%f2,partmp%at
!endif
           
  
            case (2)
                 if(numsend==0) then
                    call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts,rdispls, &
                      mpi_double_precision,comm,ierr)         
                  else
                    h=0
                    do while(associated(intmp)) 
                    if(.not.associated(intmp%next)) then
                      exit
                    else         
                      sbuf2nd(numnode*h)=intmp%numpoint
                      call para_obtain_interpolation_elefield_per_per_gyro(intmp%x(1:2), muind,elef, pic2d)
                      call compute_deri_of_sqgyep_per_per(deri_driftsquare,intmp%x(1:2),pic2d)
                      call compute_deri_of_magfield_per_per(deri_bf,intmp%x(1:2),pic2d)
                      mag=para_compute_spl2d_field_point_per_per(intmp%x(1:2),pic2d%para2d%m_x1,pic2d%para2d%m_x2,row, &
                        pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e,pic2d%field2d%bf03wg_n, &
                        pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se, &
                        pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw) 
                      sbuf2nd(numnode*h+1:numnode*h+2)=elef(1:2)
                      sbuf2nd(numnode*h+3:numnode*h+4)=deri_bf(1:2)
                      sbuf2nd(numnode*h+5:numnode*h+6)=deri_driftsquare(1:2)
                      sbuf2nd(numnode*h+7)=mag
                      h=h+1
                      intmp=>intmp%next
                    end if
                    enddo
                    call mpi_alltoallv(sbuf2nd,scounts,sdispls,mpi_double_precision,rbuf2nd,rcounts,rdispls,&
                      mpi_double_precision,comm,ierr)
                  end if

                  magf(1)=0.0_f64
                  magf(2)=0.0_f64
                  elef(3)=0.0_f64
                  do i=0,numnode-1
                  h=1 
                  partmp=>partlist
                  do while(h.ne.NINT(rbuf2nd(numnode*i))) 
                  partmp=>partmp%next
                  h=h+1
                  if(.not.associated(partmp)) then
                    print*, "#ERROR:The information of the current point is lost."
                    stop
                  end if
                  end do
                  elef(1:2)=rbuf2nd(numnode*i+1:numnode*i+2)
                  deri_bf(1:2)=rbuf2nd(numnode*i+3:numnode*i+4)
                  deri_driftsquare(1:2)=rbuf2nd(numnode*i+5:numnode*i+6)
                  magf(3)=rbuf2nd(numnode*i+7)

                  if(order==2) then
                    call gyrork_f_per_per_2nd(partmp%f2,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                  else if(order==3) then
                    call gyrork_f_per_per_2nd(partmp%f3,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                  else if(order==4) then
                    call gyrork_f_per_per_2nd(partmp%f4,partmp%coords(3),elef,deri_bf,magf,deri_driftsquare)
                  end if

                  end do


             case default
               stop
             end select

        case("nat_per")
           print*, "#ERROR: The current version doesn't include the nat_per boundary condition."
           stop
        
        case default
          stop
        end select

     case("polar")
        print*, "#ERROR: The current version doesn't include the polar geometry condition."
        stop

    case default
       stop
    end select


     deallocate(sbuf2nd,rbuf2nd)
     deallocate(scounts,rcounts,sdispls,rdispls)
     nullify(partmp,intmp)
!call mpi_barrier(comm)

   end subroutine compute_f_of_points_out_orbit_gyro2d

   subroutine gyrork4solve(gy2d_head,pic2d,iter_num)
  ! It includes the magnetic field perturbation and the electric field perturbation
     class(pic_para_total2d_base), pointer :: pic2d
     class(gy2d_node), pointer, intent(inout) :: gy2d_head
     class(gy2d_node), pointer  :: tmp
     class(rk4gy2dnode), pointer :: partlist,partmp
     class(pointin_gy_node), pointer :: inlist, intmp
     int4, intent(in) :: iter_num
     int4 :: prank,size,rank,numgr,gyroorder
     int4, dimension(:), pointer :: num,recnum
     real8 :: dtgy,coords(2)
     int4 :: i,j,h,order,comm

     rank=pic2d%layout2d%collective%rank
     dtgy=pic2d%para2d%dtgy
     comm=pic2d%layout2d%collective%comm
   
    gyroorder=pic2d%para2d%gyroorder
    size=pic2d%layout2d%collective%size
    dtgy=pic2d%para2d%dtgy
    allocate(num(0:size-1),recnum(0:size-1))
    allocate(partlist,inlist) 

    tmp=>gy2d_head
    partmp=>partlist

    do while(associated(tmp)) 
       if(.not.associated(tmp%next)) then
         exit
       else
          partmp%coords(1:3)=tmp%coords(1:3)
          partmp%vec(1:2)=tmp%coords(1:2)
          partmp%at=0
          tmp=>tmp%next
          allocate(partmp%next)
          partmp=>partmp%next     
       end if
    end do

    order=1
    call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,1,order,gyroorder,iter_num)

    num=0
    recnum=0
    order=2
    call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,1,iter_num)
    call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,1,order,gyroorder,iter_num)
    call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,1,order,gyroorder,iter_num)

    deallocate(inlist)
    nullify(inlist)
    allocate(inlist)    
    num=0
    recnum=0
    order=3
    call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,1,iter_num)
    call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,1,order,gyroorder,iter_num)
    call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,1,order,gyroorder,iter_num)



    deallocate(inlist)
    nullify(inlist)
    allocate(inlist)
    num=0
    recnum=0
    order=4
    call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,1,iter_num)
    call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,1,order,gyroorder,iter_num)
    call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,1,order,gyroorder,iter_num)    


    tmp=>gy2d_head
    partmp=>partlist

    do while(associated(tmp).and.associated(partmp))
       if(.not.associated(tmp%next).or..not.associated(partmp%next)) then
         exit
       else
          tmp%coords(1:2)=partmp%coords(1:2)+(dtgy/6.0_f64)*(partmp%f1(1:2)+2.0_f64*partmp%f2(1:2) &
                          +2.0_f64*partmp%f3(1:2)+partmp%f4(1:2))
          tmp=>tmp%next
          partmp=>partmp%next
       end if
    end do

      deallocate(recnum,num,inlist,partlist)
      nullify(partlist)
      nullify(inlist)
      nullify(partmp)
      nullify(tmp)
  
60  end subroutine gyrork4solve


   subroutine gyrork4solveallmu(gy2dmu_head,pic2d,iter_num)
  ! It includes the magnetic field perturbation and the electric field
  ! perturbation
        class(pic_para_total2d_base), pointer :: pic2d
        class(gy2dmu_node),dimension(:), pointer, intent(inout) :: gy2dmu_head
        int4, intent(in) :: iter_num
        class(gy2dmu_node), dimension(:), pointer :: gy2dmutmp
        class(rk4gy2dnode), pointer :: partlist,partmp
        class(pointin_gy_node), pointer :: inlist, intmp
        int4 :: prank,size,rank,numgr,gyroorder
        int4 :: mu_num
        int4, dimension(:), pointer :: num,recnum
        real8 :: dtgy,coords(2)
        int4 :: i,j,h,order,comm,k,parnum=0

        size=pic2d%layout2d%collective%size
        rank=pic2d%layout2d%collective%rank
        dtgy=pic2d%para2d%dtgy
        comm=pic2d%layout2d%collective%comm
        mu_num=pic2d%para2d%mu_num

        gyroorder=pic2d%para2d%gyroorder
        dtgy=pic2d%para2d%dtgy
        allocate(num(0:size-1),recnum(0:size-1))
        allocate(gy2dmutmp(1:mu_num))
        do k=1, mu_num
           if(.not.associated(gy2dmu_head(k)%ptr)) then
             print*, "#error: gy2dmu_head(k)%ptr is not initialized, k=",k
             stop
           else
             gy2dmutmp(k)%ptr=>gy2dmu_head(k)%ptr 
           endif
        enddo 

        do k=1,mu_num
          allocate(partlist,inlist)
 !         parnum=0
          partmp=>partlist
          do while(associated(gy2dmutmp(k)%ptr))
            if(.not.associated(gy2dmutmp(k)%ptr%next)) then
              exit
            else
 !         parnum=parnum+1 
              partmp%coords(1:3)=gy2dmutmp(k)%ptr%coords(1:3)
              partmp%vec(1:2)=gy2dmutmp(k)%ptr%coords(1:2)
              partmp%at=0
              gy2dmutmp(k)%ptr=>gy2dmutmp(k)%ptr%next
              allocate(partmp%next)
              partmp=>partmp%next
            end if
          end do

          order=1
          call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,k,order,gyroorder,iter_num)

          num=0
          recnum=0
          order=2
          call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,k,iter_num)
          call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,k,order,gyroorder,iter_num)
          call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,k,order,gyroorder,iter_num)

          deallocate(inlist)
          nullify(inlist)
          allocate(inlist)
          num=0
          recnum=0
          order=3
 
!          partmp=>partlist
!          do while(associated(partmp))
!            if(.not.associated(partmp%next)) then
!              exit
!            else
!              print*,partmp%f2(1:2)
!              partmp=>partmp%next
!            endif
!          enddo         
          call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,k,iter_num)
          call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,k,order,gyroorder,iter_num)
          call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,k,order,gyroorder,iter_num)

          deallocate(inlist)
          nullify(inlist)
          allocate(inlist)
          num=0
          recnum=0
          order=4
          call sortposition_by_process_gyro2d(num,recnum,partlist,inlist,pic2d,dtgy,order,k,iter_num)

          call compute_f_of_points_out_orbit_gyro2d(num,recnum,inlist,partlist,pic2d,k,order,gyroorder,iter_num)
          call compute_f_of_points_in_orbit_gyro2d(partlist,pic2d,k,order,gyroorder,iter_num)

          gy2dmutmp(k)%ptr=>gy2dmu_head(k)%ptr
          partmp=>partlist

          do while(associated(gy2dmutmp(k)%ptr).and.associated(partmp))
            if(.not.associated(gy2dmutmp(k)%ptr%next).or..not.associated(partmp%next)) then
              exit
            else
              gy2dmutmp(k)%ptr%coords(1:2)=partmp%coords(1:2)+(dtgy/6.0_f64)*(partmp%f1(1:2)+2.0_f64*partmp%f2(1:2) &
                        +2.0_f64*partmp%f3(1:2)+partmp%f4(1:2))

!if(iter_num==56.or.iter_num==55.or.iter_num==54.or.iter_num==57) then
!print*, "iter_num=",iter_num,"coords=",partmp%f1,partmp%f2
!print*, "f1=",partmp%f1,"f2=",partmp%f2,"f3=",partmp%f3, "f4=",partmp%f4
!endif
              gy2dmutmp(k)%ptr=>gy2dmutmp(k)%ptr%next
              partmp=>partmp%next
            end if
          end do

          deallocate(inlist,partlist)
          nullify(partlist)
          nullify(inlist)
        nullify(partmp)
 
!if(k==4) then
!print*,"rank=",rank
!endif

        end do
!print*, rank

        deallocate(recnum,num)
        deallocate(gy2dmutmp)

   end subroutine gyrork4solveallmu

   subroutine gyrork4solveallmu_and_sort(gy2dmu_head,pic2d,iter_num)
      class(pic_para_total2d_base), pointer :: pic2d
      class(gy2dmu_node),dimension(:), pointer, intent(inout) :: gy2dmu_head
      int4, intent(in) :: iter_num
      class(gy2dsend_node), dimension(:),pointer :: gy2dsend_head
      int4 :: csize, mu_num, i,j
      int4, dimension(:), pointer :: num

      csize=pic2d%layout2d%collective%size
      mu_num=pic2d%para2d%mu_num
      allocate(num(0:csize-1))
      
      do i=1,mu_num
        allocate(gy2dsend_head(0:csize-1))
        num=0
        do j=0,csize-1
          allocate(gy2dsend_head(j)%ptr)
        end do
        call sort_particles_among_ranks_gyallmu(gy2dmu_head,gy2dsend_head,i,pic2d, num) 
        call mpi2d_alltoallv_send_particle_2d_gy(gy2dmu_head,gy2dsend_head,num,i,pic2d)
      enddo 

     deallocate(num,gy2dsend_head) 
   end subroutine gyrork4solveallmu_and_sort

   subroutine para_obtain_interpolation_elefield_per_per_ful(x1, elef, pic2d)

    class(pic_para_total2d_base), pointer,intent(in) :: pic2d
    real8, dimension(2),intent(in) :: x1
    real8, dimension(3), intent(inout) :: elef
    int4 :: row,rank
    real8 :: deri_firstorder(2)
    row=pic2d%para2d%row 
    rank=pic2d%layout2d%collective%rank
    call para_spl2d_firstorder_derivatve_point_per_per(deri_firstorder,x1, &
      pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
      row, pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e, &
      pic2d%field2d%epwg_n,pic2d%field2d%epwg_s,pic2d%field2d%epwg_sw,pic2d%field2d%epwg_se,&
      pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw) 

    elef(1)=-deri_firstorder(1)
    elef(2)=-deri_firstorder(2)
    elef(3)= 0.0_f64


   end subroutine para_obtain_interpolation_elefield_per_per_ful

   subroutine para_obtain_interpolation_elefield_per_per_gyro(x1, muind, elef, pic2d)

     class(pic_para_total2d_base), pointer,intent(in) :: pic2d
     real8, dimension(2),intent(in) :: x1
     real8, dimension(3), intent(inout) :: elef
     int4, intent(in) :: muind
     int4 :: row, num1,num2
     real8 :: deri_firstorder(2)
     real8,dimension(:,:),pointer :: buf,box,re,rs,rw,rn,rne,rse,rsw,rnw
     int4 :: ierr

       num1=size(pic2d%field2d%ep,1)
       num2=size(pic2d%field2d%ep,2)
       row=pic2d%para2d%row
       allocate(buf(num1,num2))
       allocate(box(num1,num2),rw(num1,row),re(num1,row),rn(row,num2),rs(row,num2),&
             rsw(row,row),rse(row,row),rnw(row,row),rne(row,row),stat=ierr)

     row=pic2d%para2d%row
     box=pic2d%field2d%epgy_weight(muind,:,:)
     rw=pic2d%field2d%epgywg_w(muind,:,:)
     re=pic2d%field2d%epgywg_e(muind,:,:)
     rn=pic2d%field2d%epgywg_n(muind,:,:)
     rs=pic2d%field2d%epgywg_s(muind,:,:)
     rsw=pic2d%field2d%epgywg_sw(muind,:,:)
     rse=pic2d%field2d%epgywg_se(muind,:,:)
     rne=pic2d%field2d%epgywg_ne(muind,:,:)
     rnw=pic2d%field2d%epgywg_nw(muind,:,:)

     call para_spl2d_firstorder_derivatve_point_per_per(deri_firstorder,x1, &
          pic2d%para2d%m_x1,pic2d%para2d%m_x2, row, &
          box,rw,re,rn,rs,rsw,rse,rne,rnw)


!     call para_spl2d_firstorder_derivatve_point_per_per(deri_firstorder,x1, &
!      pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
!      row, pic2d%field2d%epgy_weight(muind,:,:),pic2d%field2d%epgywg_w(muind,:,:),pic2d%field2d%epgywg_e(muind,:,:), &
!      pic2d%field2d%epgywg_n(muind,:,:),pic2d%field2d%epgywg_s(muind,:,:),pic2d%field2d%epgywg_sw(muind,:,:), &
!      pic2d%field2d%epgywg_se(muind,:,:),pic2d%field2d%epgywg_ne(muind,:,:),pic2d%field2d%epgywg_nw(muind,:,:)) 

     elef(1)=-deri_firstorder(1)
     elef(2)=-deri_firstorder(2)
     elef(3)= 0.0_f64
!print*, elef
  
     deallocate(buf,box,re,rs,rw,rn,rne,rse,rsw,rnw)
   end subroutine para_obtain_interpolation_elefield_per_per_gyro

   subroutine compute_deri_of_magfield_per_per(deri_bf,coords,pic2d)
     class(pic_para_total2d_base), pointer, intent(in) :: pic2d
     real8, intent(in) :: coords(2)
     real8, intent(inout) :: deri_bf(2)
     int4 :: row

     row=pic2d%para2d%row
     call para_spl2d_firstorder_derivatve_point_per_per(deri_bf,coords,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
      row, pic2d%field2d%bf03wg,pic2d%field2d%bf03wg_w,pic2d%field2d%bf03wg_e, &
      pic2d%field2d%bf03wg_n,pic2d%field2d%bf03wg_s,pic2d%field2d%bf03wg_sw,pic2d%field2d%bf03wg_se,&
      pic2d%field2d%bf03wg_ne,pic2d%field2d%bf03wg_nw)
     return 
  end subroutine

  subroutine compute_deri_of_sqgyep_per_per(deri_sqgyep,coords,pic2d)
     class(pic_para_total2d_base), pointer, intent(in) :: pic2d
     real8, intent(in) :: coords(2)
     real8, intent(inout) :: deri_sqgyep(2)
     int4 :: row

     row=pic2d%para2d%row
     call para_spl2d_firstorder_derivatve_point_per_per(deri_sqgyep,coords,pic2d%para2d%m_x1,pic2d%para2d%m_x2, &
      row, pic2d%field2d%epgysq_weight,pic2d%field2d%epgysqwg_w,pic2d%field2d%epgysqwg_e, &
      pic2d%field2d%epgysqwg_n,pic2d%field2d%epgysqwg_s,pic2d%field2d%epgysqwg_sw,pic2d%field2d%epgysqwg_se,&
      pic2d%field2d%epgysqwg_ne,pic2d%field2d%epgysqwg_nw)
     return
  end subroutine

    subroutine mpi2d_alltoallv_send_points_orbit(num,numgr,numrecv,outhead,pointhead,pic2d,order)
       class(pic_para_total2d_base), pointer, intent(in) :: pic2d
       class(fuloutnode),dimension(:), pointer,intent(in) :: outhead
       class(fuloutnode),dimension(:), pointer :: sendtmp
       class(pointin_node), pointer, intent(inout) :: pointhead
       class(pointin_node), pointer :: pointmp
       int4, intent(in) :: numgr,order
       int4, dimension(:), pointer, intent(in) :: num
       int4, dimension(:), pointer, intent(inout) :: numrecv
       int4 :: comm,ierr,numout,numsend,size,rank
       int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls,rbuf0
       real8, dimension(:), pointer :: sbuf,rbuf
       int4 :: i,j,h

       rank=pic2d%layout2d%collective%rank
       size=pic2d%layout2d%collective%size
       comm=pic2d%layout2d%collective%comm
       allocate(sendtmp(0:size-1))
       do i=0,size-1
          sendtmp(i)%ptr=>outhead(i)%ptr
       end do
       allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1),rbuf0(0:size-1))
       rcounts=0
       scounts=0
       sdispls=0
       rdispls=0
       rbuf0=0
       call mpi_alltoall(num,1,mpi_integer,rbuf0,1,mpi_integer,comm,ierr)
       numrecv=rbuf0  !!!????
       numrecv(rank)=0
       call prep_for_mpi_alltoallv_with_zeroinput(rank,size,numgr,num,rbuf0,scounts, &
             rcounts,sdispls,rdispls,numsend,numout)
       if(numsend==0) then
         allocate(sbuf(0:0))
       else
         allocate(sbuf(0:numgr*numsend-1))
       end if
       sbuf=0.0
       if(numout==0) then
         allocate(rbuf(0:0))
       else
         allocate(rbuf(0:numgr*numout-1))
       end if
       rbuf=0.0
       if(numsend==0) then
         goto 80
       else
          h=0
          do i=0, size-1
            if(size==rank) then
              goto 70
            else 
              do while(associated(sendtmp(i)%ptr))
                if(.NOT.associated(sendtmp(i)%ptr%next)) then
                  exit
                else
                  sbuf(numgr*h)=real(sendtmp(i)%ptr%numpoint,8)
                  sbuf(numgr*h+1:numgr*h+numgr-1)=sendtmp(i)%ptr%coords(1:numgr-1)
                  sendtmp(i)%ptr=>sendtmp(i)%ptr%next
                  h=h+1
                end if
              end do
70          end if
          end do
80      end if

       call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
       if(.not.associated(pointhead)) then
         print*, "pointheadtype is not allocated"
         stop
       end if
       pointmp=>pointhead
       do while(associated(pointmp))
          if(associated(pointmp%next)) then
            pointmp=>pointmp%next
          else
            exit
          end if
       end do
       
       if(numout==0) then
         goto 20
       else      
         do i=0,numout-1
           pointmp%numpoint=rbuf(numgr*i)
           pointmp%x(1:numgr-1)=rbuf(numgr*i+1:numgr*i+numgr-1)
           allocate(pointmp%next)
           pointmp=>pointmp%next
         end do
20       end if

   end subroutine mpi2d_alltoallv_send_points_orbit


subroutine mpi2d_alltoallv_send_points_orbit_gy2d(num,numgr,numrecv,outhead,pointhead,pic2d)
  class(pic_para_total2d_base), pointer, intent(in) :: pic2d
  class(gyoutnode),dimension(:), pointer,intent(in) :: outhead
  class(gyoutnode),dimension(:), pointer :: sendtmp
  class(pointin_gy_node), pointer, intent(in) :: pointhead
  class(pointin_gy_node), pointer :: pointmp
  int4, intent(in) :: numgr
  int4, dimension(:), pointer, intent(in) :: num
  int4, dimension(:), pointer, intent(inout) :: numrecv
  int4 :: comm,ierr,numout,numsend,size,rank
  int4, dimension(:),pointer :: rcounts,scounts,sdispls,rdispls,rbuf0
  real8, dimension(:), pointer :: sbuf,rbuf
  int4 :: i,j,h

  rank=pic2d%layout2d%collective%rank
  size=pic2d%layout2d%collective%size
  comm=pic2d%layout2d%collective%comm
  allocate(sendtmp(0:size-1))
  do i=0,size-1
  sendtmp(i)%ptr=>outhead(i)%ptr
  end do
  allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1), &
           rbuf0(0:size-1))
  rcounts=0
  scounts=0
  sdispls=0
  rdispls=0
  rbuf0=0
  call mpi_alltoall(num,1,mpi_integer,rbuf0,1,mpi_integer,comm,ierr)

  numrecv=rbuf0
  numrecv(rank)=0
  call prep_for_mpi_alltoallv_with_zeroinput(rank,size,numgr,num,rbuf0,scounts, &
             rcounts,sdispls,rdispls,numsend,numout)
  if(numsend==0) then
    allocate(sbuf(0:0))
  else
    allocate(sbuf(0:numgr*numsend-1))
  end if
  sbuf=0.0
  if(numout==0) then
    allocate(rbuf(0:0))
  else
    allocate(rbuf(0:numgr*numout-1))
  end if
  rbuf=0.0

  if(numsend==0) then
    goto 110
  else
    h=0
    do i=0, size-1
       if(size==rank) then
         goto 120
       else 
         do while(associated(sendtmp(i)%ptr))
           if(.not.associated(sendtmp(i)%ptr%next)) then
             exit
           else
             sbuf(numgr*h)=real(sendtmp(i)%ptr%numpoint,8)
             sbuf(numgr*h+1:numgr*h+numgr-1)=sendtmp(i)%ptr%coords(1:numgr-1)
             sendtmp(i)%ptr=>sendtmp(i)%ptr%next
             h=h+1
           end if
         end do
120    end if
    end do
110 end if

   call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls, &
                        mpi_double_precision,comm,ierr)
    if(.not.associated(pointhead)) then
       print*, "pointheadtype is not allocated"
       stop
    end if
    pointmp=>pointhead
    do while(associated(pointmp))
      if(associated(pointmp%next)) then
        pointmp=>pointmp%next
      else
        exit
      end if
    end do

    if(numout==0) then
       goto 130
    else    
      do i=0,numout-1
          pointmp%numpoint=rbuf(numgr*i)
          pointmp%x(1:numgr-1)=rbuf(numgr*i+1:numgr*i+numgr-1)
          allocate(pointmp%next)
          pointmp=>pointmp%next
      end do
130    end if

  end subroutine mpi2d_alltoallv_send_points_orbit_gy2d



end module m_para_orbit
