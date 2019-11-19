!!!! don't use pic2d dummy varaible
module m_parautilities
#include "work_precision.h"
#include "gp_assert.h"
  use m_collective
  use m_mpilayout, only: t_layout_2d, &
                         get_layout_2d_box_index
  use utilities_module, only: gp_error
!  use MPI, only: mpi_alltoall, mpi_alltoallv, mpi_integer, &
!                 mpi_double_precision, mpi_cart_coords, &
!                 mpi_cart_rank

  use paradata_utilities, only: local_index_grid_2d, &
                               globalind_from_localind_2d, &
                               globnum_for_scatter_from_rootprocess, &
                               get_coords_from_processrank, &
                               get_rank_from_processcoords, & 
                               element_number_in_rank_per_per, &
                               dimsize_of_rank_per_per, &
                               copy_boundary_value_per_per, &
                               element_number_in_rank_per_per
  implicit none
include  "mpif.h"

  public :: mpi2d_alltoallv_box_per_per, &
            gather_field_to_rootprocess, &
            gather_field_to_rootprocess_per_per, &
            scatter_field_from_rootprocess, &
            scatter_field_from_rootprocess_per_per, &
            combine_boundfield_between_neighbor_alltoall

!!!!the order of the four edge of the 2d box
!!! the arrangement of the multiple lines needed in the communication
!!       3     
!!      *********** vertice
!!      ***********  2 
!!      ***********  n2
!!      *********** 
!!      *********** 
!!           1 n1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! The data is sent in the unti-clockwise direction



contains

!!!! rbuf and recounts are not right.
 subroutine mpi2d_alltoall_box_nat_per_2d(box,row,commold,numproc,boxindex,rank)
   int4, intent(in) :: row  ! the number of the rows at the boundary of box
                            ! needed to be sent
   real8, dimension(:,:),pointer, intent(inout) :: box
   int4,  intent(in) :: commold
   int4,  intent(in) :: numproc(2)
   int4,  intent(in) :: boxindex(4)
   int4,  intent(in) :: rank
   real8, dimension(:), pointer :: sendbuf   
   int4,  dimension(:), pointer :: sendcount  
   real8, dimension(:), pointer :: rbuf
   real8, dimension(:), pointer :: recvcount
   int4,   dimension(:), pointer :: sdispl, rdispl
   int4 :: cart_com, rankone
   int4 :: coords(2),h,k,i,j
   int4 :: rank1,rank2
   int4 :: comm,numproc1,numproc2
   int4 :: ierr, n1, n2,dimsize(2),reorder=1
   logical ::periods(0:1)
   
   periods(0)=.false.
   periods(0)=.true.
   numproc1=numproc(1)
   numproc2=numproc(2)
   allocate(sendcount(numproc1*numproc2))
   allocate(recvcount(numproc1*numproc2))
   allocate(sdispl(numproc1*numproc2))
   allocate(rdispl(numproc1*numproc2))
   
   call MPI_Cart_create(commold,2,dimsize,periods,reorder,comm,ierr)
   call MPI_Cart_coords(comm, rank, 2,coords, ierr)
   h=coords(1)
   k=coords(2)

   n1=boxindex(2)-boxindex(1)+1
   n2=boxindex(4)-boxindex(3)+1


 !!!==================  
   if(h==0.and.k==0) then   !!!! start here fi
      allocate(sendbuf(0:row*n1+row*n2), stat=ierr)  !! here, the extra 1 is from the vertice
      call gp_error(ierr,'sendbuf')
      allocate(rbuf(0:row*n1+row*n2))
      call  gp_error(ierr,'rbuf')
      do i=0,row-1
         do j=0,n2-1
            sendbuf(j+n2*i)=box(n1-i-1,j)
         end do
      end do

      do j=0,row-1
         do i=0,n1-1
            sendbuf(row*n2+j*n1+i)=box(i,n2-j-1)
         end do
      end do
      sendbuf(row*1+row*n2)=box(n1-1,n2-1)
         
      do i=0,numproc1-1
         do j=0, numproc2-1
            coords=(/i,j/)            
            if(i==0.and.j==0) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=0
               recvcount(rankone)=0
            else if(i==0.and.j==1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else if(i==1.and.j==0) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==1.and.j==1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=1
               recvcount(rankone)=1
            else
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if
        end do
     end do

     if(rank<=1) then
        sdispl(rank)=0
        rdispl(rank)=0
     else if(2<=rank.and.rank<=numproc1) then
        sdispl(rank)=n2*row
        rdispl(rank)=n2*row
     else if(rank==numproc1+1) then
        sdispl(rank)=n2*row+n1*row
        rdispl(rank)=n2*row+n1*row
     else if(rank>=numproc1+2) then
        sdispl(rank)=n2*row+n1*row+1
        rdispl(rank)=n2*row+n1*row+1
     end if   
      
      call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_PRECISION,rbuf, &
           recvcount,rdispl,MPI_DOUBLE_PRECISION,commold)

      do i=0,row-1
         do j=0,n2-1
            box(n1-i-1,j)=box(n1-i-1,j)+rbuf(j+n2*i)
         end do
      end do
      do j=0,row-1
         do i=0,n1-1
            box(i,n2-j-1)=box(i,n2-j-1)+rbuf(row*n2+j*n1+i)
         end do
      end do      
      box(n1-1,n2-1)=box(n1-1,n2-1)+rbuf(row*1+row*n2)
      
      deallocate(rbuf)
      deallocate(sendbuf)
      
!!!=======================
   else if(h==numproc1-1.and.k==0) then
      allocate(sendbuf(0:row*n1+row*n2),stat=ierr)  !! here, the extra 1 is from the vertice
      call gp_error(ierr,'sendbuf')
      allocate(rbuf(0:row*n1+row*n2),stat=ierr)
      call gp_error(ierr,'rbuf')

      do j=0,row-1
         do i=0,n1-1
            sendbuf(i+n1*j)=box(i,n2-j-1)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(i,j)
         end do
      end do     
      sendbuf(row*1+row*n2)=box(0,n2-1)
      
      do i=0,numproc1-1
         do j=0, numproc2-1
             coords=(/i,j/)
            if(i==numproc1-2.and.j==0) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=n2*row
               recvcount(rankone)=n2*row
            else if(i==numproc1-1.and.j==1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==numproc1-2.and.j==1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if   
        end do
     end do
     
     if(rank<=numproc1-2) then
        sdispl(rank)=0
        rdispl(rank)=0
     else if(numproc1-1<=rank.and.rank<=2*numproc1-3) then
        sdispl(rank)=n2*row
        rdispl(rank)=n2*row
     else if(rank==2*numproc1-2) then
        sdispl(rank)=n2*row+1
        rdispl(rank)=n2*row+1
     else if(rank>numproc1-2) then
        sdispl(rank)=n2*row+n1*row+1
        rdispl(rank)=n2*row+n1*row+1        
     end if
     
     call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
          recvcount,rdispl,MPI_DOUBLE_precision,commold)

      do j=0,row-1
         do i=0,n1-1
            box(i,n2-j-1)=box(i,n2-j-1)+rbuf(i+n1*j)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            box(i,j)=box(i,j)+sendbuf(row*n1+i*n2+j)
         end do
      end do     
      box(0,n2-1)=box(0,n2-1)+sendbuf(row*1+row*n2)

     deallocate(sendbuf)
     deallocate(rbuf)
     
!!!!!==========================
     else if(h==numproc1-1.and.k==numproc2-1) then
        allocate(sendbuf(0:row*n1+row*n2),stat=ierr)  !! here, the extra 1 is from the vertice
        call gp_error(ierr,'sendbuf')
        allocate(rbuf(0:row*n1+row*n2))
        call gp_error(ierr,'rbuf')
      
      do j=0,row-1
         do i=0,n1-1
            sendbuf(i+n1*j)=box(i,j)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(i,j)
         end do
      end do     
      sendbuf(row*1+row*n2)=box(0,0)

      do i=0,numproc1-1
         do j=0, numproc2-1
            coords=(/i,j/)
            if(i==numproc1-2.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=1
               recvcount(rankone)=1
            else if(i==numproc1-1.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==numproc1-2.and.j==numproc2-1) then
               call MPI_Cart_rank(cart_com, coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if   
        end do
     end do      
         
     if(rank<=numproc1*(numproc2-1)-2) then
        sdispl(rank)=0
        rdispl(rank)=0
     else if(rank==(numproc1*(numproc2-1)-1)) then
        sdispl(rank)=1
        rdispl(rank)=1
     else if(numproc1*(numproc2-1)<=rank.and.rank  &
             <=numproc1*numproc2-2)  then
        sdispl(rank)=n1*row+1
        rdispl(rank)=n1*row+1
     else if(rank==(numproc1*numproc2-1)) then
        sdispl(rank)=n2*row+n1*row+1
        rdispl(rank)=n2*row+n1*row+1
     end if
     
     call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
           recvcount,rdispl,MPI_DOUBLE_precision,commold)      

      do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(i+n1*j)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            box(i,j)=box(i,j)+rbuf(row*n1+i*n2+j)
         end do
      end do     
      box(0,0)=box(0,0)+rbuf(row*1+row*n2)
      
     deallocate(sendbuf)
     deallocate(rbuf)

!!!!!==========================
    else if(h==0.and.k==numproc2-1) then
      allocate(sendbuf(0:row*n1+row*n2),stat=ierr)  !! here, the extra 1 is from the vertice
      call gp_error(ierr,'sendbuf')
      allocate(rbuf(0:row*n1+row*n2),stat=ierr)
       call gp_error(ierr,'rbuf')     
      do j=0,row-1
         do i=0,n1-1
            sendbuf(i+n1*j)=box(i,j)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(n1-i-1,j)
         end do
      end do     
      sendbuf(row*1+row*n2)=box(n1-1,0)

      do i=0,numproc1-1
         do j=0, numproc2-1
            coords=(/i,j/)
            if(i==0.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=n1*row
               recvcount(rankone)=n1*row
            else if(i==1.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==1.and.j==numproc2-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if   
        end do
     end do      
         
     if(rank<=numproc1*(numproc2-2)) then
        sdispl(rank)=0
        rdispl(rank)=0
     else if(rank==(numproc1*(numproc2-2)+1)) then
        sdispl(rank)=n1*row
        rdispl(rank)=n1*row
     else if((numproc1*(numproc2-2)+2)<=rank.and. &
             rank<=(numproc1*(numproc2-1)+1)) then
        sdispl(rank)=n1*row+1
        rdispl(rank)=n1*row+1
     else if(rank>=numproc1*(numproc2-1)+2) then
        sdispl(rank)=n2*row+n1*row+1
        rdispl(rank)=n2*row+n1*row+1
     end if
     
     call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
          recvcount,rdispl,MPI_DOUBLE_precision,commold)

      do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(i+n1*j)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            box(n1-i-1,j)=box(n1-i-1,j)+rbuf(row*n1+i*n2+j)
         end do
      end do     
      box(n1-1,0)=box(n1-1,0)+rbuf(row*1+row*n2)

     deallocate(sendbuf)
     deallocate(rbuf)     

!!!===============================
  else if(k==0.and.(h.ne.0).and.(h.ne.numproc1-1)) then
     allocate(sendbuf(0:2*n2*row+n1*row+1),stat=ierr)
      call gp_error(ierr,'sendbuf')     
     allocate(rbuf(0:2*n2*row+n1*row+1),stat=ierr)
      call gp_error(ierr,'rbuf')
      do i=0,row-1
         do j=0,n2-1
            sendbuf(n2*i+j)=box(numproc1-i-1,j)
         end do
      end do
      do j=0,row-1
         do i=0,n1-1
            sendbuf(row*n2+j*n1+i)=box(i,numproc2-j-1)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n2+row*n1+i*n2+j)=box(i,j)
         end do
      end do     
      sendbuf(row*1+2*row*n2)=box(n1-1,n2-1)
      sendbuf(row*1+2*row*n2+1)=box(0,n2-1)

      do i=0,numproc1-1
         do j=0, numproc2-1
            coords=(/i,j/)      
            if(i==h+1.and.j==k) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=n2*row
               recvcount(rankone)=n2*row
            else if(i==h+1.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==h.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row  
            else if(i==h-1.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==h-1.and.j==k) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row  
            else
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if   
        end do
     end do

     if(rank<=h) then
        sdispl(rank)=0
        rdispl(rank)=0
     else if(rank==h.or.rank==h+1) then
        sdispl(rank)=n2*row
        rdispl(rank)=n2*row
     else if(rank.ge.h+2.and.rank.le.numproc1+h-1)  then
        sdispl(rank)=n1*row+n1*row
        rdispl(rank)=n1*row+n1*row
     else if(rank==numproc1+h) then
        sdispl(rank)=2*n2*row+1
        rdispl(rank)=2*n2*row+1
     else if(rank==numproc1+h+1) then
        sdispl(rank)=2*n2*row+n1*row+1
        rdispl(rank)=2*n2*row+n1*row+1
     else if(rank.ge.numproc1+h+2) then
        sdispl(rank) = 2*n2*row+n1*row+2
        rdispl(rank) = 2*n2*row+n1*row+2
     end if
     
     call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
           recvcount,rdispl,MPI_DOUBLE_precision,commold)      

      do i=0,row-1
         do j=0,n2-1
            box(numproc1-i-1,j)=box(numproc1-i-1,j)+rbuf(n2*i+j)
         end do
      end do
      do j=0,row-1
         do i=0,n1-1
            box(i,numproc2-j-1)=box(i,numproc2-j-1)+rbuf(row*n2+j*n1+i)
         end do
      end do
      do i=0,row-1
         do j=0,n2-1
            box(i,j)=box(i,j)+rbuf(row*n2+row*n1+i*n2+j)
         end do
      end do     
      box(n1-1,n2-1)=box(n1-1,n2-1)+rbuf(row*1+2*row*n2)
      box(0,n2-1)=box(0,n2-1)+rbuf(row*1+2*row*n2+1)
     
     deallocate(sendbuf)
     deallocate(rbuf)     
     

!!!===============================
   else if(h==numproc1-1.and.(k.ne.0).and.(k.ne.numproc2-1)) then
      allocate(sendbuf(0:2*n1*row+n2*row+1),stat=ierr)
      call gp_error(ierr,'sendbuf')
      allocate(rbuf(0:2*n1*row+n2*row+1),stat=ierr)
      call gp_error(ierr,'rbuf')

       do j=0,row-1
         do i=0,n1-1
            sendbuf(n1*j+i)=box(i,j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(numproc1-i-1,j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            sendbuf(row*n2+row*n1+j*n1+i)=box(i,numproc2-j-1)
         end do
       end do     
       sendbuf(row*1+2*row*n2)=box(0,n2-1)
       sendbuf(row*1+2*row*n2+1)=box(0,0)

       do i=0,numproc1-1
          do j=0, numproc2-1
            coords=(/i,j/)
            if(i==numproc1-2.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=1
               recvcount(rankone)=1
            else if(i==numproc1-1.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==numproc1-2.and.j==k) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else if(i==numproc1-2.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==numproc1-1.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=n1*row
               recvcount(rankone)=n1*row
            else 
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0 
            end if   
         end do
       end do

       if(rank<=numproc1*k-2) then
          sdispl(rank)=0
          rdispl(rank)=0
       else if(rank==numproc1*k-1) then
          sdispl(rank)=1
          rdispl(rank)=1
       else if(rank.gt.numproc1*k-1.and.rank.le.numproc1*(k+1)-2)  then
          sdispl(rank)=n1*row+1
          rdispl(rank)=n1*row+1
       else if(rank.ge.numproc1*(k+1)-1.and.rank.le.numproc1*(k+2)-2) then
          sdispl(rank)=n1*row+n2*row+1
          rdispl(rank)=n1*row+n2*row+1
       else if(rank==numproc1*(k+2)-1) then
          sdispl(rank)=n2*row+n1*row+2
          rdispl(rank)=n2*row+n1*row+2
       else if(rank.gt.numproc1*(k+2)-1) then
          sdispl(rank) = n2*row+2*n1*row+2
          rdispl(rank) = n2*row+2*n1*row+2
       end if
       
       call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
            recvcount,rdispl,MPI_DOUBLE_precision,commold)
 
       do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(n1*j+i)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(numproc1-i-1,j)=box(numproc1-i-1,j)+rbuf(row*n1+i*n2+j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            box(i,numproc2-j-1)=box(i,numproc2-j-1)+rbuf(row*n2+row*n1+j*n1+i)
         end do
       end do         
       box(0,n2-1)=box(0,n2-1)+rbuf(row*1+2*row*n2)
       box(0,0)=box(0,0)+rbuf(row*1+2*row*n2+1)
       
       deallocate(sendbuf)
       deallocate(rbuf)


!!!===============================
   else if(k==numproc2-1.and.(h.ne.0).and.(h.ne.numproc1-1)) then
      allocate(sendbuf(0:2*n2*row+n1*row+1),stat=ierr)
      call gp_error(ierr,'rbuf')      
      allocate(rbuf(0:2*n2*row+n1*row+1),stat=ierr)
      call gp_error(ierr,'rbuf')      

       do j=0,row-1
         do i=0,n1-1
            sendbuf(n1*j+i)=box(i,j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(numproc1-i-1,j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n2+row*n1+i*n2+j)=box(i,j)
         end do
       end do     
       sendbuf(row*1+2*row*n2)=box(n1-1,0)
       sendbuf(row*1+2*row*n2+1)=box(0,0)

       do i=0,numproc1-1 
          do j=0, numproc2-1
            coords=(/i,j/)
            if(i==h-1.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=1
               recvcount(rankone)=1
            else if(i==h.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==h+1.and.j==numproc2-2) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==h+1.and.j==numproc2-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else if(i==h-1.and.j==numproc2-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=n2*row
               recvcount(rankone)=n2*row
            else 
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0 
            end if   
         end do
       end do

       if(rank<=numproc1*(numproc2-2)+h-1) then
          sdispl(rank)=0
          rdispl(rank)=0
       else if(rank==numproc1*(numproc2-2)+h) then
          sdispl(rank)=1
          rdispl(rank)=1
       else if(rank==numproc1*(numproc2-2)+h+1)  then
          sdispl(rank)=n1*row+1
          rdispl(rank)=n1*row+1
       else if(rank.ge.numproc1*(numproc2-2)+h+2.and.rank.le.numproc1*(numproc2-1)+h-1) then
          sdispl(rank)=n1*row+2
          rdispl(rank)=n1*row+2
       else if(rank==numproc1*(numproc2-1)+h.or.rank==numproc1*(numproc2-1)+h+1) then
          sdispl(rank)=n2*row+n1*row+2
          rdispl(rank)=n2*row+n1*row+2
       else if(rank.ge.numproc1*(numproc2-1)+h+2) then
          sdispl(rank) =2*n2*row+n1*row+2
          rdispl(rank) =2*n2*row+n1*row+2
       end if
       
       call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
           recvcount,rdispl,MPI_DOUBLE_precision,commold)      

       do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(n1*j+i)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(numproc1-i-1,j)=box(numproc1-i-1,j)+rbuf(row*n1+i*n2+j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(i,j)=box(i,j)+rbuf(row*n2+row*n1+i*n2+j)
         end do
       end do     
       box(n1-1,0)=box(n1-1,0)+rbuf(row*1+2*row*n2)
       box(0,0)=box(0,0)+rbuf(row*1+2*row*n2+1)      
       
       deallocate(sendbuf)     
       deallocate(rbuf)       



!!!===============================
    else if(h==0.and.(k.ne.0).and.(k.ne.numproc2-1)) then
       allocate(sendbuf(0:2*n1*row+n2*row+1),stat=ierr)
       call gp_error(ierr,'sendbuf')       
       allocate(rbuf(0:2*n1*row+n2*row+1),stat=ierr)
       call gp_error(ierr,'rbuf')       

       do j=0,row-1
         do i=0,n1-1
            sendbuf(n1*j+i)=box(i,j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(numproc1-i-1,j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            sendbuf(row*n2+row*n1+j*n1+i)=box(i,numproc2-j-1)
         end do
       end do     
       sendbuf(row*1+2*row*n2)=box(n1-1,0)
       sendbuf(row*1+2*row*n2+1)=box(n1-1,n2-1)

       do i=0,numproc1-1
          do j=0, numproc2-1
            coords=(/i,j/)           
            if(i==h.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=n1*row
               recvcount(rankone)=n1*row
            else if(i==h+1.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==h+1.and.j==k) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else if(i==h.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==h+1.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=1
               recvcount(rankone)=1
            else 
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=0
               recvcount(rankone)=0
            end if   
         end do
       end do

       if(rank<=numproc1*(k-1)) then
          sdispl(rank)=0
          rdispl(rank)=0
       else if(rank==numproc1*(k-1)+1) then
          sdispl(rank)=n1*row
          rdispl(rank)=n1*row
       else if(rank.ge.numproc1*(k-1)+2.and.rank.le.numproc1*k+1)  then
          sdispl(rank)=n1*row+1
          rdispl(rank)=n1*row+1
       else if(rank.ge.numproc1*k+2.and.rank.le.numproc1*(k+1)) then
          sdispl(rank)=n1*row+n2*row+1
          rdispl(rank)=n1*row+n2*row+1
       else if(rank==numproc1*(k+1)+1) then
          sdispl(rank)=n2*row+2*n1*row+1
          rdispl(rank)=n2*row+2*n1*row+1
       else if(rank.ge.numproc1*(k+1)+2) then
          sdispl(rank) =n2*row+2*n1*row+2
          rdispl(rank) =n2*row+2*n1*row+2
       end if
       
       call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
            recvcount,rdispl,MPI_DOUBLE_precision,commold)

       do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(n1*j+i)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(numproc1-i-1,j)=box(numproc1-i-1,j)+rbuf(row*n1+i*n2+j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            box(i,numproc2-j-1)=box(i,numproc2-j-1)+rbuf(row*n2+row*n1+j*n1+i)
         end do
       end do     
       box(n1-1,0)=box(n1-1,0)+rbuf(row*1+2*row*n2)     
       box(n1-1,n2-1)=box(n1-1,n2-1)+rbuf(row*1+2*row*n2+1)
       
       deallocate(sendbuf)     
       deallocate(rbuf)


     !!!===============================
    else 
       allocate(sendbuf(0:2*n1*row+2*n2*row+3),stat=ierr)
       call gp_error(ierr,'sendbuf')       
       allocate(rbuf(0:2*n1*row+2*n2*row+3),stat=ierr)
       call gp_error(ierr,'rbuf') 
       do j=0,row-1
         do i=0,n1-1
            sendbuf(n1*j+i)=box(i,j)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(row*n1+i*n2+j)=box(numproc1-i-1,j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            sendbuf(row*n2+row*n1+j*n1+i)=box(i,numproc2-j-1)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            sendbuf(2*n1*row+n2*row+i*n2+j)=box(i,j)
         end do
       end do     
       sendbuf(2*row*1+2*row*n2)=box(0,0)
       sendbuf(2*row*1+2*row*n2+1)=box(n1-1,0)
       sendbuf(2*row*1+2*row*n2+2)=box(n1-1,n2-1)
       sendbuf(2*row*1+2*row*n2+3)=box(0,n2-1)      

       do i=0,numproc1-1
          do j=0, numproc2-1
            coords=(/i,j/)           
            if(i==h-1.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)=0
               recvcount(rankone)=0
            else if(i==h.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n1*row
               recvcount(rankone)= n1*row
            else if(i==h+1.and.j==k-1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= 1
               recvcount(rankone)= 1
            else if(i==h+1.and.j==k) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)
               sendcount(rankone)= n2*row
               recvcount(rankone)= n2*row
            else if(i==h+1.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=1
               recvcount(rankone)=1
            else if(i==h.and.j==k+1) then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=n1*row
               recvcount(rankone)=n1*row
            else if(i==h-1.and.j==k+1)  then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=1
               recvcount(rankone)=1
            else if(i==h-1.and.j==k)  then
               call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
               sendcount(rankone)=n2*row
               recvcount(rankone)=n2*row
            else
              call MPI_Cart_rank(cart_com,coords,rankone,ierr)      
              sendcount(rankone)= 0
              recvcount(rankone)= 0
            end if   
         end do
      end do

      if(rank<=numproc1*(k-1)+h-1) then
         sdispl(rank)=0
         rdispl(rank)=0
      else if(rank==numproc1*(k-1)+h) then
         sdispl(rank)=1
         rdispl(rank)=1
      else if(rank==numproc1*(k-1)+h+1) then
         sdispl(rank)=n1*row+1
         rdispl(rank)=n1*row+1 
      else if(rank.ge.numproc1*(k-1)+h+2.and.rank.le.numproc1*k+h-1)  then
         sdispl(rank)=n1*row+2
         rdispl(rank)=n1*row+2
      else if(rank==numproc1*k+h.or.rank==numproc1*k+h+1) then
         sdispl(rank)=n1*row+n2*row+2
         rdispl(rank)=n1*row+n2*row+2
      else if(rank.ge.numproc1*k+h+2.and.rank.le.numproc1*(k+1)+h-1)  then
         sdispl(rank)=n1*row+2*n2*row+2
         rdispl(rank)=n1*row+2*n2*row+2
      else if(rank==numproc1*(k+1)+h) then
         sdispl(rank) =n1*row+2*n2*row+3
         rdispl(rank) =n1*row+2*n2*row+3
      else if(rank==numproc1*(k+1)+h+1) then
         sdispl(rank) =2*n1*row+2*n2*row+3
         rdispl(rank) =2*n1*row+2*n2*row+3
      else
         sdispl(rank) = 2*n1*row+2*n2*row+4
         rdispl(rank) = 2*n1*row+2*n2*row+4
      end if
      
      call MPI_Alltoallv(sendbuf,sendcount,sdispl,MPI_DOUBLE_precision,rbuf, &
           recvcount,rdispl,MPI_DOUBLE_precision,commold)

       do j=0,row-1
         do i=0,n1-1
            box(i,j)=box(i,j)+rbuf(n1*j+i)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(numproc1-i-1,j)=box(numproc1-i-1,j)+rbuf(row*n1+i*n2+j)
         end do
       end do
       do j=0,row-1
         do i=0,n1-1
            box(i,numproc2-j-1)=box(i,numproc2-j-1)+rbuf(row*n2+row*n1+j*n1+i)
         end do
       end do
       do i=0,row-1
         do j=0,n2-1
            box(i,j)=box(i,j)+rbuf(2*n1*row+n2*row+i*n2+j)
         end do
       end do     
       box(0,0)=box(0,0)+rbuf(2*row*1+2*row*n2)
       box(n1-1,0)=box(n1-1,0)+rbuf(2*row*1+2*row*n2+1)
       box(n1-1,n2-1)=box(n1-1,n1-1)+rbuf(2*row*1+2*row*n2+2)
       sendbuf(2*row*1+2*row*n2+3)=box(0,n2-1)
       box(0,n2-1)=box(0,n2-1)+rbuf(2*row*1+2*row*n2+3)

      deallocate(sendbuf)     
      deallocate(rbuf)    

    end if   ! end here fi


    deallocate(sendcount)
    deallocate(recvcount)
    deallocate(sdispl)
    deallocate(rdispl)

  end subroutine  mpi2d_alltoall_box_nat_per_2d 

subroutine mpi2d_alltoallv_box_per_per(row,commold,rank,numproc, &
                                        box,rw,re,rn,rs,rsw,rse,rnw,rne,boxindex)
 !  int4 :: row  ! the number of the rows at the boundary of box
    ! needed to be sent
   int4,intent(in) :: commold,rank,row
   int4,dimension(:),intent(in) :: boxindex
   int4,dimension(2), intent(in) ::numproc
   real8, dimension(:,:),pointer, intent(inout) :: box, rw,re,rn,rs,rsw,rse,rne,rnw
   int4, dimension(:), pointer :: scounts,rcounts
   real8, dimension(:), pointer :: rbuf,sbuf
   int4,   dimension(:), pointer :: sdispls, rdispls
   int4 :: comm=1, rankone,size
   int4 :: coords(2),h,l,k,i,j,h1,l1,i1,j1
   int4 :: rank1,rank2,numproc1,numproc2
   int4 :: n1,n2,ierr,num1,num2

   if(numproc(1).lt.3.or.numproc(2).lt.3) then
      print*, "The number of the processes for each dimension shoud be equal to or great than 3 "
      stop
   end if
   numproc1=numproc(1)
   numproc2=numproc(2)
   size=numproc(1)*numproc(2)
   
   allocate(scounts(0:size-1))
   allocate(rcounts(0:size-1))
   allocate(sdispls(0:size-1))
   allocate(rdispls(0:size-1))
   scounts=0
   rcounts=0
   sdispls=0
   rdispls=0

  coords=get_coords_from_processrank(rank,numproc)

   h=coords(1)
   k=coords(2)
   n1=boxindex(2)-boxindex(1)+1   ! the nodes number of x1 dimension of box
   n2=boxindex(4)-boxindex(3)+1   ! the nodes number of x2 dimension of box
   
   scounts(rank)=0
   coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/) 
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=row**2

   coords=(/h,modulo(k-1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=row*n1

   coords=(/modulo(h+1,numproc1),modulo(k-1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc) 
   scounts(rank1)=row**2

   coords=(/modulo(h-1,numproc1),k/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=n2*row

   coords=(/modulo(h+1,numproc1),k/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=n2*row

   coords=(/modulo(h-1,numproc1),modulo(k+1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=row**2

   coords=(/h,modulo(k+1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=n1*row

   coords=(/modulo(h+1,numproc1),modulo(k+1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc)
   scounts(rank1)=row**2

   call  mpi_alltoall(scounts,1,mpi_integer, &
                      rcounts,1,mpi_integer, &
                      commold,ierr)


!   if(rank==0) then
!      print*, "scounts=", scounts
!      print*, "rcounts=", rcounts
!   end if
                   
   do i=0, size-1
      if(i==0) then
         sdispls(0)=0
         rdispls(0)=0
      else
         sdispls(i)=sdispls(i-1)+scounts(i-1)
         rdispls(i)=rdispls(i-1)+rcounts(i-1)
      end if
   end do

!if(rank==0) then
!   print*, "sdispls=",sdispls
!   print*, "rdispls=", rdispls
!end if 
!!!!!!!!!!!==================
   
       allocate(sbuf(0:2*n1*row+2*n2*row+4*row**2-1),stat=ierr)
       call gp_error(ierr,'sendbuf')       
       allocate(rbuf(0:2*n1*row+2*n2*row+4*row**2-1),stat=ierr)
       call gp_error(ierr,'rbuf')
       sbuf=0.0
       rbuf=0.0

!!!!!! Here, we need to introduce the cases of "numproc2=2 or !numproc1=2" in the future.
       coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/)
       rank1=get_rank_from_processcoords(coords,numproc)
       do j=0,row-1
         do i=0,row-1
            sbuf(sdispls(rank1)+row*j+i)=box(i+1,j+1)
         end do
      end do
     coords=(/h,modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)    
      do j=0,row-1
         do i=0,n1-1
          sbuf(sdispls(rank1)+j*n1+i)=box(i+1,j+1)
         end do
      end do

    coords=(/modulo(h+1,numproc1),modulo(k-1,numproc2)/)
    rank1=get_rank_from_processcoords(coords,numproc)   
    do j=0,row-1
         do i=0,row-1
           sbuf(sdispls(rank1)+j*row+i)=box(n1-row+1+i,j+1)
         end do
      end do
     coords=(/modulo(h-1,numproc1),k/)
     rank1=get_rank_from_processcoords(coords,numproc)  
     do j=0,n2-1
         do i=0,row-1
           sbuf(sdispls(rank1)+j*row+i)=box(i+1,j+1)
        end do
      end do
     coords=(/modulo(h+1,numproc1),k/)
     rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,n2-1
         do i=0,row-1
           sbuf(sdispls(rank1)+j*row+i)=box(n1-row+1+i,j+1)
        end do
      end do
     coords=(/modulo(h-1,numproc1),modulo(k+1,numproc2)/)
    rank1=get_rank_from_processcoords(coords,numproc)      
     do j=0,row-1
         do i=0,row-1
           sbuf(sdispls(rank1)+j*row+i)=box(i+1,n2-row+1+j)
        end do
      end do
     coords=(/modulo(h,numproc1),modulo(k+1,numproc2)/)
    rank1=get_rank_from_processcoords(coords,numproc)   
   do j=0,row-1
         do i=0,n1-1
          sbuf(sdispls(rank1)+j*n1+i)=box(i+1,n2-row+1+j)
        end do
      end do      
     coords=(/modulo(h+1,numproc1),modulo(k+1,numproc2)/)
   rank1=get_rank_from_processcoords(coords,numproc)   
   do j=0,row-1
         do i=0,row-1
           sbuf(sdispls(rank1)+j*row+i)=box(n1-row+1+i,n2-row+1+j)
        end do
      end do    

     call MPI_Alltoallv(sbuf,scounts,sdispls,MPI_DOUBLE_precision, &
                         rbuf,rcounts,rdispls,MPI_DOUBLE_precision, &
                         commold,ierr)
    if(numproc(1).lt.2) then
         print*, "The number of processor in the first dimension is ", numproc(1)
         stop
     else if(numproc(2).lt.2) then
         print*, "The number of processor in the second dimension is ",numproc(2)
         stop
     end if        

!!********============= 
     if(numproc(1)==2) then
    !!!!!!!!=======
        if(numproc(2)==2) then
        coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/)
      rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,row-1
         do i=0,row-1
            rne(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do
     num1=rdispls(rank1)+j*row+i+1
     do j1=0,row-1
         do i1=0,row-1
            rse(i+1,j+1)=rbuf(num1+j1*row+i1)
         end do
      end do
      num1=num1+j1*row+i1+1
      do h=0,row-1   
         do l=0,row-1
            rsw(l+1,h+1)=rbuf(num1+h*row+l)
        end do
      end do
      num1=num1+h*row+l+1
      do h1=0,row-1
         do l1=0,row-1
            rnw(l+1,h+1)=rbuf(num1+h1*row+l1)
         end do
      end do

     coords=(/h,modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)    
     do j=0,row-1
         do i=0,n1-1
            re(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
        end do
      end do
     num1=rdispls(rank1)+j*n1+i+1
     do h=0,row-1
        do l=0,n1-1
            rw(l+1,h+1)=rbuf(num1+h*n1+l)           
        end do
     end do 

!     coords=(/modulo(h+1,numproc1),modulo(k-1,numproc2)/)

     coords=(/modulo(h-1,numproc1),k/)
    rank1=get_rank_from_processcoords(coords,numproc) 
      do j=0,n2-1
         do i=0,row-1
           rn(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do
      do h=0,n2-1
         do l=0,row-1
           rs(l+1,h+1)=rbuf(rdispls(rank1)+j*row+i+h*row+l)
         end do
      end do
  !!!!!!==================
     else if(numproc(2).gt.2) then
      coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/)
      rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,row-1   
         do i=0,row-1
            rnw(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do
      num1=rdispls(rank1)+j*row+i+1
      do j=0,row-1
         do i=0,row-1
            rsw(i+1,j+1)=rbuf(num1+j*row+i)
         end do
      end do
     coords=(/modulo(h-1,numproc(1)),modulo(k+1,numproc2)/)
      do j=0,row-1
         do i=0,row-1
            rne(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do
      num1= rdispls(rank1)+j*row+i+1
     do j=0,row-1
         do i=0,row-1
            rse(i+1,j+1)=rbuf(num1+j*row+i)
         end do
      end do

    coords=(/modulo(h-1,numproc1),k/)
    rank1=get_rank_from_processcoords(coords,numproc) 
    do j=0,n2-1
         do i=0,row-1
           rn(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do
    num1=rdispls(rank1)+j*row+i+1
    do j=0,n2-1
         do i=0,row-1
           rs(i+1,j+1)=rbuf(num1+j*row+i)
         end do
      end do

     coords=(/h,modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,n1-1
            rw(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
         end do
      end do

     coords=(/h,modulo(k+1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,n1-1
            re(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
         end do
      end do

     end if


!!!!!********=========
  else if(numproc(1).gt.2) then
    !!!!!!=========
    if(numproc(2)==2) then
      coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/)
      rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,row-1   
         do i=0,row-1
            rse(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do              
      num1=rdispls(rank1)+j*row+i+1
      do j=0,row-1
         do i=0,row-1
            rsw(i+1,j+1)=rbuf(num1+j*row+i)
        end do
      end do     
         
      coords=(/modulo(h+1,numproc1),modulo(k-1,numproc2)/)
      rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,row-1
         do i=0,row-1
            rne(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do
      num1=rdispls(rank1)+j*row+i+1
      do j=0,row-1
         do i=0,row-1
            rnw(i+1,j+1)=rbuf(num1+j*row+i)
        end do
      end do

     coords=(/h,modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)     
     do j=0,row-1
         do i=0,n1-1
            re(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
        end do
     end do
     num1=rdispls(rank1)+j*n1+i+1
     do j=0,row-1
         do i=0,n1-1
            rw(i+1,j+1)=rbuf(num1+j*n1+i)
        end do
     end do
 
    coords=(/modulo(h-1,numproc1),k/)
    rank1=get_rank_from_processcoords(coords,numproc) 
    do j=0,n2-1
         do i=0,row-1
           rs(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do
    
     coords=(/modulo(h+1,numproc1),k/)
     rank1=get_rank_from_processcoords(coords,numproc)   
     do j=0,n2-1
         do i=0,row-1
          rn(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do

 !!!!!!!!==============  
     else if(numproc(2).gt.2) then
      coords=(/modulo(h-1,numproc1),modulo(k-1,numproc2)/)
      rank1=get_rank_from_processcoords(coords,numproc)
      do j=0,row-1   !rcounts(rank1)-2
         do i=0,row-1
            rsw(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do

     coords=(/h,modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)     
     do j=0,row-1
         do i=0,n1-1
            rw(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
        end do
      end do

     coords=(/modulo(h+1,numproc1),modulo(k-1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,row-1
            rnw(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do

     coords=(/modulo(h-1,numproc1),k/)
    rank1=get_rank_from_processcoords(coords,numproc) 
    do j=0,n2-1
         do i=0,row-1
           rs(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do

     coords=(/modulo(h+1,numproc1),k/)
     rank1=get_rank_from_processcoords(coords,numproc)   
     do j=0,n2-1
         do i=0,row-1
        rn(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
        end do
      end do

     coords=(/modulo(h-1,numproc1),modulo(k+1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,row-1
            rse(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
      end do

     coords=(/h,modulo(k+1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,n1-1
            re(i+1,j+1)=rbuf(rdispls(rank1)+j*n1+i)
         end do
      end do

     coords=(/modulo(h+1,numproc1),modulo(k+1,numproc2)/)
     rank1=get_rank_from_processcoords(coords,numproc)
     do j=0,row-1
         do i=0,row-1
            rne(i+1,j+1)=rbuf(rdispls(rank1)+j*row+i)
         end do
     end do
  
   end if

end if

    deallocate(sbuf,rbuf,scounts,rcounts,sdispls,rdispls)       
  end subroutine mpi2d_alltoallv_box_per_per    

!  function get_rank_from_processcoords(coords,numproc)
!     int4 :: coords(2),numproc(2)
!     int4 :: get_rank_from_processcoords
!      
!     get_rank_from_processcoords=coords(2)*numproc(1)+coords(1)
!  end function
!
!  function get_coords_from_processrank(rank,numproc) result(coords)
!     int4 :: rank,numproc(2)
!     int4 :: coords(2)
! 
!     coords(1)=modulo(rank,numproc(1))
!     coords(2)=(rank-coords(1))/numproc(1)
!  end function   


  subroutine gather_field_to_rootprocess(rootfield,boxfield,rank,size,boxindex,numproc,layout2d,boundary)
 !    class(pic_para_total2d_base), pointer, intent(in) :: pic2d
     class(t_layout_2d), pointer,intent(in) :: layout2d
     int4, intent(in) :: size,boxindex(4),numproc(2),rank
     real8, dimension(:), pointer, intent(inout)     :: rootfield
     real8, dimension(:,:), pointer,intent(in) :: boxfield
     character(len=*),intent(in) :: boundary
     int4,  dimension(:), pointer :: rbuf, rdispls,rcounts
     real8, dimension(:), pointer :: sbuf2,rbuf2
     int4 :: i,j,h,ierr,numout,numout1
     int4 :: nx1,nx2,numproc1,numproc2,startgridindex(2),globalind(2),localind(2)
     int4 :: length(2), globnum, num,root, comm
     
     comm=layout2d%collective%comm
     numout=0
     root=0

        numout= (layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1)*  &
                (layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)
     allocate(rbuf(0:size-1), stat=ierr)
      rbuf=0
     call mpi_gather(3*numout,1,mpi_integer,rbuf,1,mpi_integer,root,comm,ierr)

     allocate(sbuf2(0:3*numout-1))
     sbuf2=0._F64
     allocate(rcounts(0:size-1),rdispls(0:size-1))
     rcounts=0
     h=0
     do i=1, boxindex(4)-boxindex(3)+1
          do j=1, boxindex(2)-boxindex(1)+1
             globalind=globalind_from_localind_2d((/j,i/),numproc,rank,layout2d,boundary)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h-1)=boxfield(i,j)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h)=real(globalind(1),8)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h+1)=real(globalind(2),8)
             h=h+2
          end do
     end do

        numout1=0
     do i=0,size-1
        numout1=numout1+rbuf(i)
        rcounts(i)=rbuf(i)
        if(i==0) then
           rdispls(i)=0
        else
           rdispls(i)=rdispls(i-1)+rcounts(i-1)
        end if
    end do
     allocate(rbuf2(0:numout1-1),stat=ierr)
     call mpi_gatherv(sbuf2,3*numout,mpi_double_precision,rbuf2,rcounts,rdispls,mpi_double_precision, &
                     root,comm,ierr)

     do i=0,size-1
        numout=numout+(layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min+1)*  &
                (layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min+1)
     end do
  
     if(rank==0) then
          h=0
          do i=0,numout-1
             globalind(1)=NINT(rbuf2(i+h+1))
             globalind(2)=NINT(rbuf2(i+h+2))
             globnum=(globalind(2)-1)*layout2d%global_sz1+globalind(1)
             rootfield(globnum)=rbuf2(i+h)
             h=h+2
          end do
      end if
   
      deallocate(rbuf,rbuf2,rcounts,sbuf2)
   end subroutine gather_field_to_rootprocess


subroutine gather_field_to_rootprocess_per_per(rootfield,boxfield,rank,size,boxindex,numproc,layout2d)
 !    class(pic_para_total2d_base), pointer, intent(in) :: pic2d
     class(t_layout_2d), pointer,intent(in) :: layout2d
     int4, intent(in) :: size,boxindex(4),numproc(2),rank
     real8, dimension(:), pointer, intent(inout)     :: rootfield
     real8, dimension(:,:), pointer,intent(in) :: boxfield
     int4,  dimension(:), pointer :: rbuf, rdispls,rcounts
     real8, dimension(:), pointer :: sbuf2,rbuf2
     character(len=25) :: boundary="double_per"
     int4 :: i,j,h,ierr,numout,numout1
     int4 :: nx1,nx2,numproc1,numproc2,startgridindex(2),globalind(2),localind(2)
     int4 :: length(2), globnum, num,root, comm,coords(2)
     
     comm=layout2d%collective%comm
     numout=0
     root=0
     coords=get_coords_from_processrank(rank,numproc)
     numout=element_number_in_rank_per_per(rank,numproc,layout2d)
     allocate(rbuf(0:size-1), stat=ierr)
      rbuf=0
     call mpi_gather(3*numout,1,mpi_integer,rbuf,1,mpi_integer,root,comm,ierr)

     allocate(sbuf2(0:3*numout-1))
     sbuf2=0._F64
     allocate(rcounts(0:size-1),rdispls(0:size-1))
     rcounts=0
     if(coords(1)==numproc(1)-1.and.coords(2).ne.numproc(2)-1) then
     h=0
     do i=1, boxindex(4)-boxindex(3)+1
          do j=1, boxindex(2)-boxindex(1)
             globalind=globalind_from_localind_2d((/j,i/),numproc,rank,layout2d,boundary)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h-1)=boxfield(j,i)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h)=real(globalind(1),8)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h+1)=real(globalind(2),8)
             h=h+2
          end do
     end do
     else if(coords(1).ne.numproc(1)-1.and.coords(2)==numproc(2)-1) then
     h=0
     do i=1, boxindex(4)-boxindex(3)
          do j=1, boxindex(2)-boxindex(1)+1
             globalind=globalind_from_localind_2d((/j,i/),numproc,rank,layout2d,boundary)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h-1)=boxfield(j,i)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h)=real(globalind(1),8)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h+1)=real(globalind(2),8)
             h=h+2
          end do
     end do  
     else if(coords(1)==numproc(1)-1.and.coords(2)==numproc(2)-1)   then
     h=0
     do i=1, boxindex(4)-boxindex(3)
          do j=1, boxindex(2)-boxindex(1)
             globalind=globalind_from_localind_2d((/j,i/),numproc,rank,layout2d,boundary)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h-1)=boxfield(j,i)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h)=real(globalind(1),8)
             sbuf2((i-1)*(boxindex(2)-boxindex(1))+j+h+1)=real(globalind(2),8)
             h=h+2
          end do
     end do
     else
     h=0
     do i=1, boxindex(4)-boxindex(3)+1
          do j=1, boxindex(2)-boxindex(1)+1
             globalind=globalind_from_localind_2d((/j,i/),numproc,rank,layout2d,boundary)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h-1)=boxfield(j,i)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h)=real(globalind(1),8)
             sbuf2((i-1)*(boxindex(2)-boxindex(1)+1)+j+h+1)=real(globalind(2),8)
             h=h+2
          end do
     end do
     end if

  !    print*, "rank=",rank,"sbuf2=",sbuf2
     
        numout1=0
     do i=0,size-1
        numout1=numout1+rbuf(i)
        rcounts(i)=rbuf(i)
        if(i==0) then
           rdispls(i)=0
        else
           rdispls(i)=rdispls(i-1)+rcounts(i-1)
        end if
    end do
     allocate(rbuf2(0:numout1-1),stat=ierr)
     call mpi_gatherv(sbuf2,3*numout,mpi_double_precision,rbuf2,rcounts,rdispls,mpi_double_precision, &
                     root,comm,ierr)
  
     if(rank==0) then
          h=0
          do i=0,numout1/3-1
             globalind(1)=NINT(rbuf2(i+h+1))
             globalind(2)=NINT(rbuf2(i+h+2))
             globnum=(globalind(2)-1)*layout2d%global_sz1+globalind(1)
             rootfield(globnum)=rbuf2(i+h)
             h=h+2
          end do
      end if
   
      deallocate(rbuf,rbuf2,rcounts,sbuf2)
   end subroutine gather_field_to_rootprocess_per_per



!   subroutine coordinates_in_para_mesh_per_per(x,pic2d)
!      class(pic_para_total2d_base), intent(in) :: pic2d
!      real8,intent(2) :: x(2)
!      real8 :: gboxmin(2),gboxmax(2),xbar
!
!      gboxmin(:)=pic2d%para2d%gboxmin(0,:)
!      gboxmax(1)=pic2d%para2d%gboxmax(numproc-1,1)
!      gboxmax(2)=pic2d%para2d%gboxmax(size-1,2)
!      if(x(1).lt.gboxmin(1)) then
!         xbar=(gboxmin(1)-x(1))-real(floor((gboxmin(1)-x(1))/(gboxmax(1)-gboxmin(1))),8)*(gboxmax(1)-gboxmin(1))
!         x(1)=gboxmax(1)-xbar  !!!! Here, we use all the mesh in positive number
!      else if(x(1).gt.gboxmax(1)) then
!         xbar=x(1)-gboxmax(1)-real(floor((gboxmin(1)-x(1))/(gboxmax(1)-gboxmin(1))),8)*(gboxmax(1)-gboxmin(1))
!         x(1)=gboxmin(1)+xbar
!      end if
!      if(x(2).lt.gboxmin(2)) then
!         xbar=(gboxmin(2)-x(2))-real(floor((gboxmin(2)-x(2))/(gboxmax(2)-gboxmin(2))),8)*(gboxmax(2)-gboxmin(2))
!         x(2)=gboxmax(2)-xbar  !!!! Here, we use all the mesh in positive number
!      else if(x(2).gt.gboxmax(2)) then
!         xbar=x(2)-gboxmax(2)-real(floor((gboxmin(2)-x(2))/(gboxmax(2)-gboxmin(2))),8)*(gboxmax(2)-gboxmin(2))
!         x(2)=gboxmin(2)+xbar
!      end if
!
!   end subroutine coordinates_in_para_mesh_per_per



   subroutine  scatter_field_from_rootprocess(rootfield,boxfield,size,numproc,global_sz,layout2d,boundary)
     class(t_layout_2d), pointer,intent(in) :: layout2d
     int4, intent(in) :: size,numproc(2),global_sz(2)
     real8, dimension(:), pointer, intent(in)     :: rootfield
     real8, dimension(:,:), pointer,intent(inout) :: boxfield
     character(len=*),intent(in) :: boundary
     real8,  dimension(:), pointer :: sbuf,rbuf
     int4,   dimension(:), pointer :: scounts,sdispls
     int4 :: i,j,h,ierr,num,myrank,rcounts,globnum
     int4 :: nx1,nx2,numproc1,numproc2,startgridindex(2),globalind(2),boxindex(4)
     int4 :: rank,localind(2),rankone,root,comm,boxsize(2), length(2)
     
     num=0
     do i=0,size-1
        num=num+(layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min+1)*  &
                (layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min+1)
     end do
     allocate(sbuf(0:num-1))
     allocate(scounts(0:size-1),sdispls(0:size-1))
     do i=0,size-1
        scounts(i)=(layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min+1)* &
                   (layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min+1)
        if(i==0) then
          sdispls(i)=0
        else 
          sdispls(i)=sdispls(i-1)+scounts(i-1)
        end if
     end do
     myrank=layout2d%collective%rank
     rcounts=(layout2d%boxes(myrank)%i_max-layout2d%boxes(myrank)%i_min+1)* &
                   (layout2d%boxes(myrank)%j_max-layout2d%boxes(myrank)%j_min+1)
     allocate(rbuf(0:rcounts-1))  
     
     j=0
     do i=0,size-1
        length(1)=layout2d%boxes(i)%i_max-layout2d%boxes(i)%i_min+1
        length(2)=layout2d%boxes(i)%j_max-layout2d%boxes(i)%j_min+1
        h=1
        do while(h.le.scounts(i))
            localind=local_index_grid_2d(h,length) 
            globalind=globalind_from_localind_2d(localind,numproc,i,layout2d,boundary)
            globnum=global_sz(1)*(globalind(2)-1)+globalind(1)  
            sbuf(j)=rootfield(globnum)
            h=h+1
            j=j+1
        end do
     end do
     root=0
     comm=layout2d%collective%comm
     call mpi_scatterv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts, &
                       mpi_double_precision,root,comm,ierr) 
     boxsize(1)=layout2d%boxes(myrank)%i_max-layout2d%boxes(myrank)%i_min+1
     boxsize(2)=layout2d%boxes(myrank)%j_max-layout2d%boxes(myrank)%j_min+1
     do j=1,  boxsize(2)  
        do i=1, boxsize(1)
          boxfield(i,j)=rbuf((j-1)*boxsize(1)+i-1)
        end do
     end do
   
   end subroutine scatter_field_from_rootprocess 
  

 !!!! For per_per boundary condition, the global_sz here equals the origional global_sz minus (/1,1/) 
   subroutine  scatter_field_from_rootprocess_per_per(rootfield,boxfield,size,numproc,global_sz,layout2d)
     class(t_layout_2d), pointer,intent(in) :: layout2d
     int4, intent(in) :: size,numproc(2),global_sz(2)
     real8, dimension(:), pointer, intent(in)     :: rootfield
     real8, dimension(:,:), pointer,intent(inout) :: boxfield
     real8,  dimension(:), pointer :: sbuf,rbuf,sbuf1,rbuf1,sbuf2,rbuf2
     int4,   dimension(:), pointer :: scounts,sdispls
     character(len=25) :: boundary="double_per"
     int4 :: i,j,h,ierr,num,numrank,myrank,rcounts,globnum
     int4 :: nx1,nx2,numproc1,numproc2,startgridindex(2),globalind(2),boxindex(4)
     int4 :: rank,localind(2),rankone,root,comm,boxsize(2), length(2),coords(2),votex,coords1(2)
     int4 :: dest,status
     real8 :: vortex

     root=0
     comm=layout2d%collective%comm    
     num=0
     numrank=0
     do i=0,size-1
        numrank=element_number_in_rank_per_per(i,numproc,layout2d)
        num=num+numrank
     end do

     allocate(sbuf(0:num-1))
     allocate(scounts(0:size-1),sdispls(0:size-1))
     scounts=0
     sdispls=0
     sbuf=0.0
     do i=0,size-1 
        scounts(i)=element_number_in_rank_per_per(i,numproc,layout2d)
        if(i==0) then
          sdispls(i)=0
        else 
          sdispls(i)=sdispls(i-1)+scounts(i-1)
        end if
     end do
     myrank=layout2d%collective%rank
     rcounts=scounts(myrank)
     allocate(rbuf(0:rcounts-1))  
     rbuf=0.0 
     j=0
     do i=0,size-1
       length=dimsize_of_rank_per_per(i,numproc,layout2d)     
        h=1
        do while(h.le.scounts(i))
            localind=local_index_grid_2d(h,length) 
            globalind=globalind_from_localind_2d(localind,numproc,i,layout2d,boundary)
            globnum=global_sz(1)*(globalind(2)-1)+globalind(1)  
            sbuf(j)=rootfield(globnum)
            h=h+1
            j=j+1
        end do
     end do

     call mpi_scatterv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,mpi_double_precision,root,comm,ierr) 
     boxsize=dimsize_of_rank_per_per(myrank,numproc,layout2d)
     do j=1,  boxsize(2)  
        do i=1, boxsize(1)
          boxfield(i,j)=rbuf((j-1)*boxsize(1)+i-1)
        end do
     end do
     coords=get_coords_from_processrank(rank,numproc)
     if(rank==0) then
        coords1(1)=0
        coords1(2)=numproc(2)-1
        dest=get_rank_from_processcoords(coords1,numproc)
        call mpi_send(boxfield(1,1),1,mpi_double_precision,dest,11,comm,ierr)         
        
        coords1(1)=numproc(1)-1
        coords1(2)=0
        dest=get_rank_from_processcoords(coords1,numproc)
        call mpi_send(boxfield(1,1),1,mpi_double_precision,dest,11,comm,ierr)
     else if(coords(1)==0.and.coords(2)==numproc(2)-1) then
        call mpi_recv(votex,1,mpi_double_precision,0,11,comm,status,ierr)
        boxfield(1,layout2d%boxes(rank)%j_max-layout2d%boxes(rank)%j_min+1)=votex      
     else if(coords(1)==numproc(1)-1.and.coords(2)==0) then
        call mpi_recv(votex,1,mpi_double_precision,0,11,comm,status,ierr)
        boxfield(layout2d%boxes(rank)%i_max-layout2d%boxes(rank)%i_min+1,1)=votex
     end if
!print*,"myrank=",myrank, boxfield
    call copy_boundary_value_per_per(boxfield,myrank,numproc,layout2d)
        
    deallocate(scounts,sdispls,sbuf,rbuf) 
   end subroutine scatter_field_from_rootprocess_per_per 
 
!!!!! This is used to accumulate the value on the interpolation point to the grid point
   subroutine  combine_boundfield_between_neighbor_alltoall(box,numproc,layout2d)
      class(t_layout_2d),pointer,intent(in) :: layout2d
      real8, dimension(:,:), pointer, intent(inout) :: box
      int4, intent(in) :: numproc(2)
      real8,dimension(:,:), pointer :: re,rs,rw,rn,rne,rse,rsw,rnw
      int4 :: index(4),myrank,IERR,comm,length(2)
      int4 :: row=1
      int4 :: i,j

      comm=layout2d%collective%comm
      call MPI_COMM_RANK(comm,myrank,ierr)
      call get_layout_2d_box_index(layout2d,myrank,index)
      length(1)=index(2)-index(1)+1
      length(2)=index(4)-index(3)+1
      allocate(re(length(1),row),rs(row,length(2)),rw(length(1),row), rn(row,length(2)))
      allocate(rne(row,row),rse(row,row),rsw(row,row),rnw(row,row))
      
      call mpi2d_alltoallv_box_per_per(row,comm,myrank,numproc, &
                                        box,rw,re,rn,rs,rsw,rse,rnw,rne,index)    

      box(:,1)=box(:,1)+rw(:,1)
      box(:,length(2))=box(:,length(2))+re(:,1)
      box(1,:)=box(1,:)+rs(1,:)
      box(length(1),:)=box(length(1),:)+rn(1,:)
      box(1,1)=box(1,1)+rsw(1,1)
      box(length(1),1)=box(length(1),1)+rnw(1,1)
      box(length(1),length(2))=box(length(1),length(2))+rne(1,1)
      box(1,length(1))=box(1,length(1))+rse(1,1)


   end subroutine combine_boundfield_between_neighbor_alltoall

end module m_parautilities
