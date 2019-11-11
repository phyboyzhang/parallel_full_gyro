!For periodic boundary condition, the maxiumal line in both dimension records the same values of the first line in both
!dimensions.
module m_para_spline
#include "work_precision.h"

use spline_module
use m_mpilayout,only: t_layout_2d
use paradata_utilities,only:  get_rank_from_processcoords,  &
                              get_coords_from_processrank
use m_parautilities, only: scatter_field_from_rootprocess_per_per
use cartesian_mesh, only: cartesian_mesh_1d

implicit none
include "mpif.h"

public ::  para_compute_spl2d_field_point_nat_per, &
     para_compute_spl2d_field_point_per_per, &
     para_spl2d_firstorder_derivatve_point_per_per, &
     para_compute_spl2d_point_per_per_weight


contains
  
!!              *******************
!!    p         *  4  *  7  *  3  *
!!    e         *     *     *     * 
!!    r         *******************                
!!    i         *  8  *  9  *  6  *
!!    o         *     *     *     *
!!    d         *******************                
!!    i         *     *     *     *
!!    c         *  1  *  5  *  2  *          
!!              *******************        
!!                      natural

!!!!! This function needs modification
                      
 function para_compute_spl2d_field_point_nat_per(x,eta_min,eta_max,row,rank,NC,numproc, &
   box,rw,re,rn,rs,rsw,rse,rne,rnw,commold) result(field_value)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
   real8,intent(in) :: x(2)
   int4, intent(in) :: row,commold, rank, NC(2),numproc(2)  ! comm is for the cartecian communicator
   real8, dimension(:),intent(in) :: eta_min, eta_max
    real8,dimension(:,:), intent(in) :: box, rw,re,rn,rs,rsw,rse,rne,rnw
    real8 :: field_value

    real8 :: eta_star(2) 
    real8 :: val(-1:2,-1:2), weight
    int4 :: ii(2),ind(2),ell_1,ell_2
    int4 :: coords(2), ierr, nproc, numproc1, numproc2
    int4 :: flag, comm=2
    

    numproc1=numproc(1)
    numproc2=numproc(2)
    nproc=numproc1*numproc2
    
    coords=get_coords_from_processrank(rank,numproc) 
    call s_localize_new(x,eta_min,eta_max,ii,eta_star,NC,flag)
    call s_contribution_spl(eta_star, val)

    field_value = 0._f64
    do ell_2=-1,2
       ind(2)=ii(2)+ell_2
       do ell_1=-1,2           
          ind(1)=ii(1)+ell_1

          if(coords(1)==0 ) then
             ! natural boundary for x1,and just on the boundary points of the simulated domain
             if(ind(1)<0) then   
                ind(1)=0
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=row+ind(2)
                   weight=rs(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rn(ind(1),ind(2))
                ! within the box
                else
                   weight=box(ind(1),ind(2))
                end if

             ! out of the eastern box boundary
             else if(ind(1).gt.NC(1)-1) then
                ind(1)=ind(1)-NC(1)
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=ind(2)+row
                   weight=rse(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rne(ind(1), ind(2))
                ! within the box
                else
                   weight=re(ind(1),ind(2))
                end if                
                   
             else
                if(ind(2)==-1) then
                   ind(2)=ind(2)+row
                   weight=rs(ind(1),ind(2))
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rn(ind(1),ind(2))
                else
                   weight=box(ind(1),ind(2))
                end if
  
             end if

       else if(coords(1)==numproc1-1) then 

             if(ind(1)<0) then   
                ind(1)=ind(1)+row
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=row+ind(2)
                   weight=rsw(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rnw(ind(1),ind(2))
                ! within
                else
                   weight=rw(ind(1),ind(2))
                end if

             ! out of the eastern box boundary
             else if(ind(1).gt.NC(1)-1) then
                ind(1)=NC(1)-1
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=row+ind(2)
                   weight=rs(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rn(ind(1), ind(2))
                ! within the box
                else
                   weight=box(ind(1),ind(2))
                end if                
                   
             else
                if(ind(2)==-1) then
                   ind(2)=abs(ind(2))-1
                   weight=rs(ind(1),ind(2))
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rn(ind(1),ind(2))
                else
                   weight=box(ind(1),ind(2))
                end if
  
             end if

             
       else if(coords(1).lt.numproc1-1.and.coords(1).gt.0) then              
             if(ind(1)<0) then   
                ind(1)=ind(1)+row
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=ind(2)+row
                   weight=rsw(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rnw(ind(1),ind(2))
                ! within the box
                else
                   weight=rw(ind(1),ind(2))
                end if

             ! out of the eastern box boundary
             else if(ind(1).gt.NC(1)-1) then
                ind(1)=ind(1)-NC(1)
                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
                if(ind(2)==-1) then
                   ind(2)=ind(2)+row
                   weight=rse(ind(1),ind(2))
                ! out of the northern box boundary
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rne(ind(1), ind(2))
                ! within the box
                else
                   weight=re(ind(1),ind(2))
                end if                
                   
             else
                if(ind(2)==-1) then
                   ind(2)=ind(2)+row
                   weight=rs(ind(1),ind(2))
                else if(ind(2).gt.NC(2)-1) then
                   ind(2)=ind(2)-NC(2)
                   weight=rn(ind(1),ind(2))
                else
                   weight=box(ind(1),ind(2))
                end if
  
             end if 


!!$    else if(coords(0)==0.and. &
!!$          (coords(1).gt.0 .and.coords(1).lt.numproc2-1)  &
!!$           ) then              
!!$             if(ind(1)<0) then   
!!$                ind(1)=0
!!$                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rs(ind(1),ind(2)
!!$                ! out of the northern box boundary
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rn(ind(1),ind(2))
!!$                ! within the box
!!$                else
!!$                   weight=box(ind(1),ind(2))
!!$                end if
!!$
!!$             ! out of the eastern box boundary
!!$             else if(ind(1).gt.NC(1)-1) then
!!$                ind(1)=ind(1)-NC(1)
!!$                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rse(ind(1),ind(2)
!!$                ! out of the northern box boundary
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rne(ind(1), ind(2))
!!$                ! within the box
!!$                else
!!$                   weight=re(ind(1),ind(2))
!!$                end if                
!!$                   
!!$             else
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rs(ind(1),ind(2))
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rn(ind(1),ind(2))
!!$                else
!!$                   weight=box(ind(1),ind(2))
!!$                end if
!!$  
!!$             end if             
          
!!$          else    !!! if(coords(0)=numproc1-1.and. &
!!$                  !!!   (coords(1).gt.0.and.coords(1).lt.pic2dtype%para2d%numproc2-1) &
!!$                  !!!    ) then              
!!$             if(ind(1)<0) then   
!!$                ind(1)=-ind(1)-1
!!$                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rsw(ind(1),ind(2)
!!$                ! out of the northern box boundary
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rnw(ind(1),ind(2))
!!$                ! within the box
!!$                else
!!$                   weight=rw(ind(1),ind(2))
!!$                end if
!!$
!!$             ! out of the eastern box boundary
!!$             else if(ind(1).gt.NC(1)-1) then
!!$                ind(1)=NC(1)-1
!!$                ! periodic boudnary condition for x2; out of the sourthern simulated boundary
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rs(ind(1),ind(2)
!!$                ! out of the northern box boundary
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rn(ind(1), ind(2))
!!$                ! within the box
!!$                else
!!$                   weight=box(ind(1),ind(2))
!!$                end if                
!!$                   
!!$             else
!!$                if(ind(2)==-1) then
!!$                   ind(2)=abs(ind(2))-1
!!$                   weight=rs(ind(1),ind(2))
!!$                else if(ind(2).gt.NC(2)-1) then
!!$                   ind(2)=ind(2)-NC(2)
!!$                   weight=rn(ind(1),ind(2))
!!$                else
!!$                   weight=box(ind(1),ind(2))
!!$                end if
!!$  
!!$             end if

       end if  !!(all 9 kinds of bos are enumerated)
          
       field_value=field_value+val(ell_1,ell_2)*weight      
    end do
  end do
 
  end function para_compute_spl2d_field_point_nat_per  

!! compute the weight of the point  in the parallel environment
!!! Here, NC is the cell of the box
!!! Here, x is located within the  box labeled by the current rank
  subroutine para_compute_spl2d_point_per_per_weight(weight,m_x1,m_x2,x,row, &
    box,rw,re,rn,rs,rsw,rse,rne,rnw) 
    real8,dimension(:,:),pointer, intent(inout) :: weight
    class(cartesian_mesh_1d), pointer, intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    int4, intent(in) :: row  ! comm is for the cartesian communicator
    real8,dimension(:,:), intent(in) :: box, rw,re,rn,rs,rsw,rse,rne,rnw
    real8 :: eta_star(2), eta_min(2),eta_max(2)
    real8 :: val(-1:2,-1:2)
    int4 :: ii(2),ind(2),ell_1,ell_2
    int4 :: coords(2), ierr,NC(2)
    int4 :: flag,i,j   
 
    eta_min=(/m_x1%eta_min,m_x2%eta_min/)
    eta_max=(/m_x1%eta_max,m_x2%eta_max/)
    NC=(/m_x1%nodes,m_x2%nodes/)
 
!    coords=get_coords_from_processrank(rank,numproc)   
    call s_localize_new(x,eta_min,eta_max,ii,eta_star,NC-(/1,1/),flag)
    call s_contribution_spl(eta_star, val)

    do ell_1=-1,2
       ind(1)=ii(1)+ell_1
      !!! IN the sourthern part
       if(ind(1).le.0) then   
         ind(1)=row+ind(1)-1 
           do ell_2=-1,2           
               ind(2)=ii(2)+ell_2
              ! in the southwestern box
                if(ind(2).le.0) then
                   ind(2)=row+ind(2)-1
                   weight(ell_1,ell_2)=rsw(ind(1),ind(2))
                ! in the sourthest box
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight(ell_1,ell_2)=rse(ind(1),ind(2))
                ! in the sourthern box
                else
                  weight(ell_1,ell_2)=rs(ind(1),ind(2))
                end if
            end do
             ! in the northern part
         else if(ind(1).gt.NC(1)) then
            ind(1)=ind(1)-NC(1)+1
            do ell_2=-1,2           
                ind(2)=ii(2)+ell_2
                ! periodic boudnary condition for x2; 
                ! in the northwestern part
                if(ind(2).le.0) then
                   ind(2)=ind(2)+row-1
                   weight(ell_1,ell_2)=rnw(ind(1),ind(2))
                ! in the northwestern part
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight(ell_1,ell_2)=rne(ind(1), ind(2))
                ! in the northern part
                else
                  weight(ell_1,ell_2)=rn(ind(1),ind(2))
                end if                
             end do          
           else 
                ! the first dimension withen in box
                ! in the western part
             do ell_2=-1,2           
                ind(2)=ii(2)+ell_2
         
                if(ind(2).le.0) then
                   ind(2)=ind(2)+row-1
                   weight(ell_1,ell_2)=rw(ind(1),ind(2))
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight(ell_1,ell_2)=re(ind(1),ind(2))
                else
                   weight(ell_1,ell_2)=box(ind(1),ind(2))
               end if 
              end do 
          end if
         
    end do

!    do i=-1,2
!      do j=-1, 2
!         print*, i,j,weight(i,j)
!      end do
!    end do  

!   print*, "box=",box
!   print*, "re=", re
!   print*, "rn=", rn
!   print*, "rw=", rn
!   print*,  "rs=", rs
!   print*, "rsw=", rsw
!   print*, "rse=", rse
!   print*, "rnw=", rnw
!   print*, "rne=", rne
 
 end subroutine para_compute_spl2d_point_per_per_weight

!!! compute the field value at a point in the parallel enviroment
!!! Here, NC is the nodes number of the boxes
 function para_compute_spl2d_field_point_per_per(x,m_x1,m_x2,row, &
   box,rw,re,rn,rs,rsw,rse,rne,rnw) result(field_value)
 !   class(cartesian_mesh_1d), intent(in) :: m_x1,m_x2
    real8,intent(in) :: x(2)
    class(cartesian_mesh_1d), pointer, intent(in) :: m_x1,m_x2
    int4, intent(in) :: row ! comm is for the cartesian communicator
    real8,dimension(:,:), intent(in) :: box, rw,re,rn,rs,rsw,rse,rne,rnw
    real8 :: field_value
    real8 :: eta_star(2),eta_min(2),eta_max(2)
    real8 :: val(-1:2,-1:2), weight
    int4 :: ii(2),ind(2),ell_1,ell_2
    int4 :: ierr,Nc(2)
    int4 :: flag

    eta_min=(/m_x1%eta_min,m_x2%eta_min/)
    eta_max=(/m_x1%eta_max,m_x2%eta_max/)
    NC=(/m_x1%nodes,m_x2%nodes/)
    
    call s_localize_new(x,eta_min,eta_max,ii,eta_star,NC-(/1,1/),flag)
    call s_contribution_spl(eta_star, val)

    field_value = 0._f64
    do ell_1=-1,2           
       ind(1)=ii(1)+ell_1

          if(ind(1).le.0) then   
            ind(1)=row+ind(1)-1
            do ell_2=-1,2
                 ind(2)=ii(2)+ell_2
                 if(ind(2).le.0) then
                   ind(2)=row+ind(2)-1
                   weight=rsw(ind(1),ind(2))
                 else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=rse(ind(1),ind(2))
                 else
                   weight=rs(ind(1),ind(2))
                 end if

              field_value=field_value+val(ell_1,ell_2)*weight
             end do

           else if(ind(1).gt.NC(1)) then
             ind(1)=ind(1)-NC(1)+1
             do ell_2=-1,2
                ind(2)=ii(2)+ell_2
                if(ind(2).le.0) then
                   ind(2)=ind(2)+row-1
                   weight=rnw(ind(1),ind(2))
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=rne(ind(1), ind(2))
                else
                   weight=rn(ind(1),ind(2))
                end if      

              field_value=field_value+val(ell_1,ell_2)*weight             
             end do      

           else
             do ell_2=-1,2
                ind(2)=ii(2)+ell_2
                if(ind(2).le.0) then
                   ind(2)=ind(2)+row-1
                   weight=rw(ind(1),ind(2))
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=re(ind(1),ind(2))
                else
                   weight=box(ind(1),ind(2))
                end if

             field_value=field_value+val(ell_1,ell_2)*weight
             end do 
          end if

    end do
end function para_compute_spl2d_field_point_per_per

!!! compute the derivative at a point in the parallel environment
 subroutine para_spl2d_firstorder_derivatve_point_per_per(deri_firstorder,x,m_x1,m_x2, &
    row,box, rw,re,rn,rs,rsw,rse,rne,rnw)
    real8,intent(in) :: x(2)
    int4, intent(in) :: row
    class(cartesian_mesh_1d), pointer, intent(in) :: m_x1,m_x2 
    real8, dimension(:,:),pointer, intent(in) :: box, rw,re,rn,rs,rsw,rse,rne,rnw
    real8, dimension(:), intent(inout) :: deri_firstorder
    real8 :: field_value
    real8 :: eta_star(2),eta_min(2),eta_max(2),delta_eta(2)
    real8 :: val(-1:2,-1:2), weight
    real8 :: val_1st_deri(-1:2,-1:2), val_2nd_deri(-1:2,-1:2)   
    int4 :: ii(2),ind(2),ell_1,ell_2,NC(2)
    int4 :: flag,ierr

!    int4 :: rank
!
!    call mpi_comm_rank(mpi_comm_world,rank,ierr)

    eta_min=(/m_x1%eta_min,m_x2%eta_min/)
    eta_max=(/m_x1%eta_max,m_x2%eta_max/)
    delta_eta=(/m_x1%delta_eta,m_x2%delta_eta/)
    NC=(/m_x1%nodes,m_x2%nodes/) 

    deri_firstorder=0.0_f64
      
    call s_localize_new(x,eta_min,eta_max,ii,eta_star,NC-(/1,1/),flag)
    call compute_spl_derivative_spline_2d(eta_star,val_1st_deri,val_2nd_deri,delta_eta)     
    call s_contribution_spl(eta_star, val)


    field_value = 0._f64
    do ell_1=-1,2           
       ind(1)=ii(1)+ell_1
       if(ind(1).le.0) then   
           ind(1)=row+ind(1)-1
           do ell_2=-1,2
                ind(2)=ii(2)+ell_2
                if(ind(2).le.0) then
                   ind(2)=row+ind(2)-1
                   weight=rsw(ind(1),ind(2))
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=rse(ind(1),ind(2))
                else
                   weight=rs(ind(1),ind(2))
                end if
                deri_firstorder(1)=deri_firstorder(1)+val_1st_deri(ell_1,ell_2)*weight
                deri_firstorder(2)=deri_firstorder(2)+val_2nd_deri(ell_1,ell_2)*weight
            end do

        else if(ind(1).gt.NC(1)) then
           ind(1)=ind(1)-NC(1)+1
           do ell_2=-1,2
                ind(2)=ii(2)+ell_2
                if(ind(2).le.0) then
                   ind(2)=row+ind(2)-1
                   weight=rnw(ind(1),ind(2))
                else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=rne(ind(1), ind(2))
                else
                   weight=rn(ind(1),ind(2))
                end if                
                deri_firstorder(1)=deri_firstorder(1)+val_1st_deri(ell_1,ell_2)*weight
                deri_firstorder(2)=deri_firstorder(2)+val_2nd_deri(ell_1,ell_2)*weight 
           end do
        
        else
           do ell_2=-1,2
              ind(2)=ii(2)+ell_2
              if(ind(2).le.0) then
                   ind(2)=row+ind(2)-1
                   weight=rw(ind(1),ind(2))
              else if(ind(2).gt.NC(2)) then
                   ind(2)=ind(2)-NC(2)+1
                   weight=re(ind(1),ind(2))
              else
                   weight=box(ind(1),ind(2))
              end if
              deri_firstorder(1)=deri_firstorder(1)+val_1st_deri(ell_1,ell_2)*weight
              deri_firstorder(2)=deri_firstorder(2)+val_2nd_deri(ell_1,ell_2)*weight
           end do  
        end if
 
    end do

 end subroutine para_spl2d_firstorder_derivatve_point_per_per

!!! compute the weight with the known field and ASPL.
  subroutine para_compute_spl2D_weight(ASPL,rootfield,boxweight,numproc,layout2d,boundary)
     class(t_layout_2d),pointer,intent(in) :: layout2d
     real8, dimension(:,:), pointer, intent(in) :: ASPL
     real8, dimension(:),   pointer, intent(in) :: rootfield
     real8, dimension(:,:),   pointer, intent(inout) :: boxweight
     character(len=*), intent(in) :: boundary
     int4, intent(in) :: numproc(2)
     real8, dimension(:,:), pointer :: buf
     int4 :: global_sz(2), myrank
     int4 :: i,j,size,ierr
     global_sz(1)=layout2d%global_sz1
     global_sz(2)=layout2d%global_sz2
     size=layout2d%collective%size
     allocate(buf(global_sz(1)*global_sz(2),1))
 !    allocate(buf1(global_sz(1)*global_sz(2)))
     call MPI_COMM_RANK(layout2d%collective%comm,myrank,ierr) 
     if(myrank==0) then
       buf(:,1)=rootfield(:)
       buf(:,1)=matmul(ASPL,buf(:,1))  
     end if 
     select case(boundary)
        case ("double_per")
           call scatter_field_from_rootprocess_per_per(buf(:,1),boxweight,size,numproc,global_sz,layout2d)
        case ("nat_per")
           stop
        case default
          stop
      end select
 
      deallocate(buf) 
  end subroutine para_compute_spl2D_weight  




end module m_para_spline
