!!!! don't use pic2d as dummy variable
module paradata_utilities
#include "work_precision.h"
use paradata_type,only: pic_para_total2d_base
use m_mpilayout, only: t_layout_2d,& 
                       get_layout_2d_box_index
use piclayout,   only: parameters_set_2d
implicit none      
include "mpif.h"

public:: globalind_from_localind_2d, &
         local_index_grid_2d,  &
 !        coordinates_in_para_mesh_per_per, &
         localind_from_globalind_2d, &
         local_index_process_2d,  &
         compute_rank_from_globalind_2d, &
         compute_process_of_point_per_per, &
         startind_of_process, &
         globnum_for_scatter_from_rootprocess, &
         coordinates_pointoutbound_per_per, &
         coordinate_pointoutbound_per, &
         get_rank_from_processcoords,  &
         get_coords_from_processrank,  &
         element_number_in_rank_per_per, &
         dimsize_of_rank_per_per, &
         point_location_for_cubic_spline_2d, &
         prepare_sendrecv_for_mpialltoallv, &
         prepare_sendrecv_for_mpialltoallv_simple, &
         coords_from_localind, &
         prepare_recv_for_gatherv, &
         prep_for_mpi_alltoallv_with_zeroinput

contains

   function globalind_from_localind_2d(localindex,numproc,rank,layout2d,boundary) result(globind)
      class(t_layout_2d), pointer, intent(in) :: layout2d
      int4, intent(in) :: numproc(2)
      int4, intent(in) :: rank,localindex(2)
      character(len=*), intent(in) :: boundary
      int4 ::  globind(2)   !globalind_from_localind_2d(2)
      int4 :: h,i,nx1,nx2,startind(2)


     nx1=modulo(rank,numproc(1))
     nx2=(rank-nx1)/numproc(1)
     startind=1
     i=0  
     do while(i.lt.nx1)
        startind(1)=startind(1)+(layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min)
        i=i+1
     end do
     i=0
     do while(i.lt.nx2)
        startind(2)=startind(2)+(layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min)
        i=i+1
     end do

     globind(1)=startind(1)+localindex(1)-1
     globind(2)=startind(2)+localindex(2)-1
     select case(boundary)
     case ("double_per")
     if(globind(1)==layout2d%global_sz1+1) then
       globind(1)=1
     end if
     if(globind(2)==layout2d%global_sz2+1) then
       globind(2)=1
     end if

     case ("nat_per")
       if(globind(2)==layout2d%global_sz2+1) then
         globind(2)=1
       end if

     case default
       stop

     end select

   end function globalind_from_localind_2d

   function localind_from_globalind_2d(globind,rank,numproc,layout2d) result(localind)
      int4, intent(in) :: globind(2),rank,numproc(2)
  !    real8,dimension(:,:), pointer, intent(in) :: gboxmin,gboxmax
      class(t_layout_2d),pointer, intent(in) :: layout2d
      int4 :: localind(2)
      int4 :: nx1,nx2, i,j,startind(2)
      
  !    nx1=modulo(rank,numproc(1))
  !    nx2=(rank-nx1)/numproc(1)
      startind=startind_of_process(rank,numproc,layout2d)
      localind=globind-startind+1

   end function localind_from_globalind_2d

   function local_index_process_2d(num,length)
      int4, intent(in) :: length(2),num
      int4 :: local_index_process_2d(2)
      local_index_process_2d(1)=modulo(num,length(1))
      local_index_process_2d(2)=(num-local_index_process_2d(1))/length(1)

   end function local_index_process_2d

   !!!!! compute the localindex within the local process
   function local_index_grid_2d(num,length) !!! num begins from 1 and ends with length(1)*length(2)
      int4, intent(in) :: length(2),num
      int4 :: local_index_grid_2d(2)
      local_index_grid_2d(1)=modulo(num,length(1))
      if(local_index_grid_2d(1)==0) then
        local_index_grid_2d=length(1)
      end if
      local_index_grid_2d(2)=(num-local_index_grid_2d(1))/length(1)+1
    
   end function local_index_grid_2d

   function compute_rank_from_globalind_2d(nx1,nx2,numproc,layout2d) result(rank)
      int4, intent(in) :: nx1,nx2,numproc(2)
    !  real8,dimension(:,:), pointer, intent(in) :: gboxmin,gboxmax
      class(t_layout_2d), pointer, intent(in) :: layout2d
      int4 :: rank
      int4 :: mx1(2),mx2(2),rankone,h
     
      if(nx1==1) then
         if(nx2==1) then
           rank=0
         else
           mx2(1)=1
           mx2(2)=1
           do h=0,numproc(2)-1
              mx2(1)=mx2(2)
              mx2(2)=mx2(2)+(layout2d%boxes(h*numproc(1))%j_max-layout2d%boxes(h*numproc(1))%j_min)   
              if(nx2.gt.mx2(1).and.nx2.le.mx2(2)) then
                 rank=h*numproc(1)
                 goto 150
              end if           
           end do
150      end if
     else
         mx1(1)=1
         mx1(2)=1
         do h=0, numproc(1)-1
            mx1(1)=mx1(2)
            mx1(2)=mx1(1)+(layout2d%boxes(h)%i_max-layout2d%boxes(h)%i_min)
            if(nx1.gt.mx1(1).and.nx1.le.mx1(2)) then
               rankone=h
               goto 250
            end if
         end do
250      if(nx2==1) then
             rank=rankone
         else
            mx2(1)=1
            mx2(2)=1
            do h=0, numproc(2)-1
              mx2(1)=mx2(2)
              mx2(2)=mx2(2)+(layout2d%boxes(h*numproc(1))%j_max-layout2d%boxes(h*numproc(1))%j_min)
              if(nx2.gt.mx2(1).and. nx2.le.mx2(2)) then
                rank=h*numproc(1)+rankone
                goto 350
              end if
            end do
350       end if
        end if

  end function compute_rank_from_globalind_2d


  !!! compute the rank of the location x
  function compute_process_of_point_per_per(x,numproc,gxmin,gxmax,gboxmin,gboxmax) result(rank)
      real8, intent(inout) :: x(2)
      int4, intent(in) :: numproc(2)
      real8,dimension(:,:), pointer, intent(in) :: gboxmin,gboxmax
      real8,intent(in) ::  gxmin(2),gxmax(2)
      int4 :: h,rankone  
      int4 :: rank

     call coordinates_pointoutbound_per_per(x,gxmin,gxmax) 

            do h=0, numproc(1)-1
              if(gboxmin(h,1).le.x(1).and.gboxmax(h,1).gt.x(1)) then
                 rankone=h
                 exit
              end if
            end do

           do h=0, numproc(2)-1
                 if(gboxmin(h*numproc(1),2).le.x(2).and.gboxmax(h*numproc(2),2).gt.x(2)) then
                     rank=h*numproc(1)+rankone
                     exit
                 end if
           end do

   end function compute_process_of_point_per_per

   !!!The start global index of the process
   function startind_of_process(rank,numproc,layout2d) result(startind)
     int4, intent(in) :: rank, numproc(2)
     class(t_layout_2d), pointer, intent(in) :: layout2d
     int4 :: nx1,nx2,startind(2),i
     nx1=modulo(rank,numproc(1))
     nx2=(rank-nx1)/numproc(1)

     startind=(/1,1/)
     i=0
     do while(i.lt.nx1)
        startind(1)=startind(1)+(layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min)
        i=i+1
     end do
     i=0
     do while(i.lt.nx2)
        startind(2)=startind(2)+(layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min)
        i=i+1
     end do
    

   end function startind_of_process

!!!! Find the 
   function globnum_for_scatter_from_rootprocess(num,global_sz,gboxmin,gboxmax,numproc,layout2d) result(globnum)
     class(t_layout_2d), pointer, intent(in) :: layout2d
     int4, intent(in) :: num,global_sz(2), numproc(2)
     real8, dimension(:,:), pointer, intent(in) :: gboxmin,gboxmax
     int4 :: globnum
     int4 :: globind(2),rank,numcol,localind(2),startind(2)
     int4 :: i,h,nx1,nx2

       globind(1)=modulo(num,global_sz(1))
       globind(2)=(num-globind(1))/global_sz(1)       
       if(globind(1)==0) then
         globind(1)=global_sz(1)
       end if
       rank=compute_rank_from_globalind_2d(globind(1),globind(2),numproc,layout2d)       
       localind= localind_from_globalind_2d(globind,rank,numproc,layout2d)

       nx1=modulo(rank,numproc(1))
       nx2=(rank-nx1)/numproc(1)
         numcol=1
         h=0
          do while(h.lt.nx2)
               numcol=numcol+layout2d%boxes(h*numproc(1))%j_max-layout2d%boxes(h*numproc(1))%j_min
               h=h+1
          end do
       globnum=numcol*global_sz(1)
       h=0
       do while(h.lt.nx1)
           globnum=globnum+(layout2d%boxes(rank-h-1)%i_max-layout2d%boxes(rank-h-1)%i_min) &
                   *(layout2d%boxes(rank-h-1)%j_max-layout2d%boxes(rank-h-1)%j_min)
           h=h+1
       end do
       globnum=globnum+(localind(2)-1)*(layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min)+localind(1)
       
   end function globnum_for_scatter_from_rootprocess

   subroutine coordinates_pointoutbound_per_per(x,gxmin,gxmax)
      real8,intent(inout) :: x(2)
      real8, intent(in) :: gxmin(2),gxmax(2)
      real8 :: gboxmin(2),gboxmax(2),xbar

!      gboxmin=pic2d%para2d%gxmin
!      gboxmax=pic2d%para2d%gxmax
      if(x(1).lt.gxmin(1)) then
         xbar=(gxmin(1)-x(1))-real(floor((gxmin(1)-x(1))/(gxmax(1)-gxmin(1))),8)*(gxmax(1)-gxmin(1))
         x(1)=gxmax(1)-xbar  !!!! Here, we use all the mesh in positive number
      else if(x(1).ge.gxmax(1)) then
         xbar=x(1)-gxmax(1)-real(floor((x(1)-gxmax(1))/(gxmax(1)-gxmin(1))),8)*(gxmax(1)-gxmin(1))
         x(1)=gxmin(1)+xbar
      end if
      if(x(2).lt.gxmin(2)) then
         xbar=(gxmin(2)-x(2))-real(floor((gxmin(2)-x(2))/(gxmax(2)-gxmin(2))),8)*(gxmax(2)-gxmin(2))
         x(2)=gxmax(2)-xbar  !!!! Here, we use all the mesh in positive number
      else if(x(2).ge.gxmax(2)) then
         xbar=x(2)-gxmax(2)-real(floor((x(2)-gxmax(2))/(gxmax(2)-gxmin(2))),8)*(gxmax(2)-gxmin(2))
         x(2)=gxmin(2)+xbar
      end if

      if(x(1)==gxmax(1)) then
         x(1)=gxmin(1)
      end if

      if(x(2)==gxmax(2)) then
         x(2)=gxmin(2)
      end if

   end subroutine coordinates_pointoutbound_per_per

 function coordinate_pointoutbound_per(x,gboxmin,gboxmax)
     real8, intent(in) :: x, gboxmin,gboxmax
     real8 :: coordinate_pointoutbound_per,xbar
      if(x.lt.gboxmin) then
         xbar=(gboxmin-x)-real(floor((gboxmin-x)/(gboxmax-gboxmin)),8)*(gboxmax-gboxmin)
        coordinate_pointoutbound_per =gboxmax-xbar  !!!! Here, we use all the mesh in positive number
      else if(x.ge.gboxmax) then
         xbar=x-gboxmax-real(floor((x-gboxmax)/(gboxmax-gboxmin)),8)*(gboxmax-gboxmin)
         coordinate_pointoutbound_per=gboxmin+xbar
      end if
   end function coordinate_pointoutbound_per

  function get_rank_from_processcoords(coords,numproc)
     int4 :: coords(2),numproc(2)
     int4 :: get_rank_from_processcoords,coordsone(2)

     coordsone(1)=modulo(coords(1),numproc(1))
     coordsone(2)=modulo(coords(2),numproc(2))
     get_rank_from_processcoords=coordsone(2)*numproc(1)+coordsone(1)
  end function

  function get_coords_from_processrank(rank,numproc) result(coords)
     int4 :: rank,numproc(2)
     int4 :: coords(2)

     coords(1)=modulo(rank,numproc(1))
     coords(2)=(rank-coords(1))/numproc(1)
  end function


  function element_number_in_rank_per_per(rank,numproc,layout2d) result(num)
     int4,intent(in) :: rank, numproc(2)
     class(t_layout_2d), pointer, intent(in) :: layout2d
     int4 :: coords(2),num

        coords=get_coords_from_processrank(rank,numproc)
        if(coords(1)==numproc(1)-1.and.coords(2).ne.numproc(2)-1) then
           num= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min)*  &
                (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)
        else if(coords(1).ne.numproc(1)-1.and.coords(2)==numproc(2)-1) then
           num= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1)*  &
                (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min)
        else if(coords(1)==numproc(1)-1.and.coords(2)==numproc(2)-1)   then
           num= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min)*  &
                (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min)
        else
           num= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1)*  &
                (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)
        end if
     
  end function

  function dimsize_of_rank_per_per(rank,numproc,layout2d) result(dimsize)
     int4,intent(in) :: rank, numproc(2)
     class(t_layout_2d), pointer, intent(in) :: layout2d
     int4 :: coords(2),dimsize(2)
      
        coords=get_coords_from_processrank(rank,numproc)
        if(coords(1)==numproc(1)-1.and.coords(2).ne.numproc(2)-1) then
           dimsize(1)= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min)
           dimsize(2)= (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)
        else if(coords(1).ne.numproc(1)-1.and.coords(2)==numproc(2)-1) then
           dimsize(1) = (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1)
           dimsize(2) = (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min)
        else if(coords(1)==numproc(1)-1.and.coords(2)==numproc(2)-1)   then
           dimsize(1) = (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min)
           dimsize(2) = (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min)
        else
           dimsize(1) = (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1)
           dimsize(2) = (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)
        end if     

  end function


  subroutine copy_boundary_value_per_per(box,myrank,numproc,layout2d)
    class(t_layout_2d), pointer, intent(in) :: layout2d
    real8, dimension(:,:), pointer, intent(in) :: box
    int4, intent(in) :: myrank,numproc(2)
    real8, dimension(:), pointer :: sbuf,rbuf
    int4 :: comm, coords(2),num,ierr,dest,source,coords1(2),status
    int4 :: i,j

    comm=layout2d%collective%comm
    coords=get_coords_from_processrank(myrank,numproc)
    if(coords(2)==0) then
      num=layout2d%boxes(myrank)%i_max-layout2d%boxes(myrank)%i_min+1
      allocate(sbuf(0:num-1))
      do j=1,num
         sbuf(j-1)=box(j,1)
      end do
      coords1(1)=coords(1)
      coords1(2)=numproc(2)-1
      dest=get_rank_from_processcoords(coords1,numproc)
  !  print*,"myrank=",myrank, 0
      call mpi_send(sbuf,num,mpi_double_precision,dest,coords(1),comm,ierr)
      deallocate(sbuf)

    else if(coords(2)==numproc(2)-1) then
      num=layout2d%boxes(myrank)%i_max-layout2d%boxes(myrank)%i_min+1 
      allocate(rbuf(0:num-1)) 
      coords1(1)=coords(1)
      coords1(2)=0
      source=get_rank_from_processcoords(coords1,numproc)      
      call mpi_recv(rbuf,num,mpi_double_precision,source,coords(1),comm,status,ierr)
      do i=1,layout2d%boxes(myrank)%i_max-layout2d%boxes(myrank)%i_min+1
         box(i,num)=rbuf(i-1)
      end do
      deallocate(rbuf)
    end if 
!print*, 3
    if(coords(1)==0) then
      num=layout2d%boxes(myrank)%j_max-layout2d%boxes(myrank)%j_min+1
      allocate(sbuf(0:num-1))
      do j=1,num
         sbuf(j-1)=box(1,j)
      end do
      coords1(1)=numproc(1)-1
      coords1(2)=coords(2)
      dest=get_rank_from_processcoords(coords1,numproc)
      call mpi_send(sbuf,num,mpi_double_precision,dest,coords(2),comm,ierr)
      deallocate(sbuf)

    else if(coords(1)==numproc(1)-1) then
      num=layout2d%boxes(myrank)%j_max-layout2d%boxes(myrank)%j_min+1
      allocate(rbuf(0:num-1))
      coords1(1)=0
      coords1(2)=coords(2)
      source=get_rank_from_processcoords(coords1,numproc)
      call mpi_recv(rbuf,num,mpi_double_precision,source,coords(2),comm,status,ierr)
      do i=1,layout2d%boxes(myrank)%j_max-layout2d%boxes(myrank)%j_min+1
         box(num,i)=rbuf(i-1)
      end do
      deallocate(rbuf)
    end if    
!print*, 4
  end subroutine copy_boundary_value_per_per

!!!! get the rank and local index of a interpolation point, which is used in
!gyroaverage 
  subroutine point_location_for_cubic_spline_2d(rankpoint,coords,myrank,numproc,ind,para2d,layout2d)
    class(t_layout_2d),pointer, intent(in) :: layout2d
    class(parameters_set_2d), pointer, intent(in) :: para2d
    int4, intent(in) :: myrank,ind(2),numproc(2)
    int4, intent(out) :: rankpoint,coords(2)
    int4 :: index(4), numone(2),numtwo(2),nd(2),rankone,ndone

    nd(1)=para2d%m_x1%nodes
    nd(2)=para2d%m_x2%nodes
    numone(1)=modulo(myrank,numproc(1))
    numone(2)=(myrank-numone(1))/numproc(2)

    if(ind(1).gt.nd(1)) then
      if(nd(1).le.2) then
         print*, "the nodes of this box is not large than 2."
         stop
      end if
      numtwo(1)=modulo(numone(1)+1,numproc(1))
      coords(1)=ind(1)-nd(1)
    else if(ind(1).lt.1) then
      numtwo(1)=modulo(numone(1)-1,numproc(1))
      rankone=get_rank_from_processcoords((/numtwo(1),numone(2)/),numproc)
      call get_layout_2d_box_index(layout2d,rankone,index)
      ndone=index(2)-index(1)+1
      coords(1)= ndone-(1-ind(1))
    else
      numtwo(1)=numone(1)
      coords(1)=ind(1) 
    end if 
      
    if(ind(2).gt.nd(2)) then
      if(nd(2).le.2) then
         print*, "the nodes of this box is not large than 2."
         stop
      end if
      numtwo(2)=modulo(numone(2)+1,numproc(2))
      coords(2)=ind(2)-nd(2)
    else if(ind(2).lt.2) then
      numtwo(2)=modulo(numone(2)-1,numproc(2))
      rankone=get_rank_from_processcoords((/numone(1),numtwo(2)/),numproc)
      call get_layout_2d_box_index(layout2d,rankone,index)
      ndone=index(4)-index(3)+1
      coords(2)= ndone-(1-ind(2))
    else 
      numtwo(2)=numone(2)
      coords(2)=ind(2)
    end if    

    rankpoint=get_rank_from_processcoords(numtwo,numproc)    

  end subroutine   point_location_for_cubic_spline_2d


  subroutine  prepare_sendrecv_for_mpialltoallv(num,numgr,size,comm,scounts,rcounts,sdispls,rdispls, &
              numsend,numrecv)
     int4, dimension(:), pointer,intent(in) :: num
     int4, intent(in) :: numgr,size,comm
     int4, dimension(:), pointer,intent(inout) :: scounts,rcounts,sdispls,rdispls
     int4, intent(inout) :: numsend,numrecv

     int4, dimension(:),pointer :: rcountmp
     int4 :: i,numout,ierr

     allocate(rcountmp(0:size-1))
    rcountmp=0
    call mpi_alltoall(num,1,mpi_integer,rcountmp,1,mpi_integer,comm,ierr)
    numsend=0
    do i=0,size-1
       numsend=numsend+num(i)
    end do
    numsend=numgr*numsend
    numrecv=0
    do i=0, size-1
      numrecv=numrecv+rcountmp(i)
    end do
!    allocate(rbuf(0:numout*numgr-1))
    numrecv=numgr*numrecv
   do i=0,size-1
      scounts(i)=num(i)*numgr
      rcounts(i)=rcountmp(i)*numgr
       if(i==0) then
          sdispls(i)=0
          rdispls(i)=0
       else
          sdispls(i)=sdispls(i-1)+scounts(i-1)
          rdispls(i)=rdispls(i-1)+rcounts(i-1)
       end if
    end do      

    deallocate(rcountmp)
  end subroutine prepare_sendrecv_for_mpialltoallv

  subroutine prepare_sendrecv_for_mpialltoallv_simple(size,numgr,snum,rnum,scounts,rcounts,sdispls,rdispls, &
             numsend,numrecv)
    int4, dimension(:),pointer, intent(in) :: snum,rnum
    int4, dimension(:),pointer, intent(inout) :: scounts,rcounts,sdispls,rdispls
    int4, intent(in) :: size,numgr
    int4, intent(inout) :: numsend,numrecv
    int4 :: i
    numsend=0
    do i=0,size-1
       numsend=numsend+snum(i)
    end do
    numsend=numgr*numsend
    numrecv=0
    do i=0, size-1
      numrecv=numrecv+rnum(i)
    end do
   numrecv=numgr*numrecv
   do i=0,size-1
      scounts(i)=snum(i)*numgr
      rcounts(i)=rnum(i)*numgr
       if(i==0) then
          sdispls(i)=0
          rdispls(i)=0
       else
          sdispls(i)=sdispls(i-1)+scounts(i-1)
          rdispls(i)=rdispls(i-1)+rcounts(i-1)
       end if
    end do      

  end subroutine prepare_sendrecv_for_mpialltoallv_simple

  subroutine prep_for_mpi_alltoallv_with_zeroinput(rank,size,numgr,snum,rnum,scounts, &
             rcounts,sdispls,rdispls,numsend,numout)
    int4, dimension(:),pointer, intent(in) :: snum,rnum
    int4, dimension(:),pointer, intent(inout) :: scounts,rcounts,sdispls,rdispls
    int4, intent(in) :: size,numgr,rank
!    real8, dimension(:), pointer, intent(out) :: sbuf,rbuf
    int4,intent(inout) :: numsend,numout
    int4 :: i

   numsend=0
   do i=0,size-1
      numsend=numsend+snum(i)
   end do
   numsend=numsend-snum(rank)
   numout=0
   do i=0, size-1
      numout=numout+rnum(i)
   end do
   numout=numout-rnum(rank)

   if(numsend==0) then
!     allocate(sbuf(0:0))
     do i=0,size-1
        scounts(i)=0
        sdispls(i)=0
     end do
   else
!     allocate(sbuf(0:numsend*numgr-1))
     do i=0,size-1
        if(i==rank) then
          scounts(i)=0
        else
          scounts(i)=snum(i)*numgr
        end if
        if(i==0) then
          sdispls(i)=0
        else
          sdispls(i)=sdispls(i-1)+scounts(i-1)
        end if
     end do
   endif

   if(numout==0) then
!     allocate(rbuf(0:0))
     do i=0,size-1
        rcounts(i)=0
        rdispls(i)=0
     end do
   else
!     allocate(rbuf(0:numout*numgr-1))
     do i=0,size-1
       if(i==rank) then
         rcounts(i)=0
       else
         rcounts(i)=rnum(i)*numgr
       end if

       if(i==0) then
         rdispls(i)=0
       else
         rdispls(i)=rdispls(i-1)+rcounts(i-1)
       end if
     end do
  end if

  end subroutine  prep_for_mpi_alltoallv_with_zeroinput 



  subroutine prepare_recv_for_gatherv(numtot,numrecv,rdispls, &
             numleft,numgr,size,comm,rank)
     int4,intent(in) :: numleft,numgr,size,comm,rank
     int4, dimension(:),pointer,intent(out) :: numrecv
     int4, dimension(:),pointer, intent(inout) :: rdispls
     int4, intent(inout) :: numtot
     int4, dimension(:), pointer :: rbuf
     int4 :: root
     int4 :: ierr, i
  
     root=0
     allocate(rbuf(0:size-1))
     rbuf=0
     call mpi_gather(numleft,1,mpi_integer,rbuf,1,mpi_integer,root,comm,ierr)    
     call mpi_Bcast(rbuf,size,mpi_integer,root,comm,ierr)
!print*, "rbuf=",rbuf
     numrecv=numgr*rbuf
     rdispls(0)=0
     do i=1, size-1
        rdispls(i)=rdispls(i-1)+numrecv(i-1)          
     end do
     numtot=0 
     do i=0, size-1
        numtot=numtot+numrecv(i)
     end do

     deallocate(rbuf)
!     end if
  end subroutine prepare_recv_for_gatherv


  function coords_from_localind(localind,rank,gboxmin,delta) result(coords)
    real8 :: coords(2)
    int4,intent(in) :: rank, localind(2)
    real8,dimension(:,:),pointer, intent(in) :: gboxmin
    real8, intent(in) :: delta(2)
    
    coords(1)=gboxmin(rank,1)+real(localind(1)-1,8)*delta(1)
    coords(2)=gboxmin(rank,2)+real(localind(2)-1,8)*delta(2) 

  end function

end module paradata_utilities


