module para_write_file
#include "work_precision.h"
use m_parautilities, only: gather_field_to_rootprocess_per_per
use paradata_utilities, only: prepare_recv_for_gatherv
use paradata_type, only: pic_para_total2d_base
use piclayout, only: root_precompute_data
use m_mpilayout,only: get_layout_2d_box_index
use orbit_data_base, only: tp_ful2d_node,tp_gy2d_node,tp_gy2dallmu_node
implicit none
include "mpif.h"

  public :: open_file, &
            close_file, &
            para_write_field_file_2d, &
            para_write_orbit_file_2d, &
            para_write_orbit_file_2d_gy, &
            para_write_orbit_file_2d_gy_allmu, &
            allmu_para_write_field_file_2d            
contains


  subroutine open_file(fileitem,filepath,rank)
    int4,intent(in) :: fileitem,rank
    character(len=*), intent(in) :: filepath
    logical :: alive

    if(rank==0) then
    inquire(file=trim(filepath),exist=alive) 
    if(.not.alive) then
       open(fileitem,file=trim(filepath),status='new')
    else
       open(fileitem,file=trim(filepath),status='replace')
    end if
    end if

  end subroutine

  subroutine close_file(fileitem,rank)
    int4, intent(in) :: fileitem,rank

    if(rank==0) then
      close(fileitem)
    end if
  end subroutine


  subroutine para_write_field_file_2d(boxfield,fileitem,rootdata,pic2d)
    real8, dimension(:,:),pointer,intent(in) :: boxfield
    class(root_precompute_data), pointer, intent(in) :: rootdata
    int4,intent(in) :: fileitem
    class(pic_para_total2d_base), pointer, intent(in) :: pic2d
    int4 :: rank,size,boxindex(4),global_sz(2),globind
    real8 :: gxmin(2),delta(2),x(2)
    int4 :: i,j

    rank=pic2d%layout2d%collective%rank
    size=pic2d%layout2d%collective%size
    call get_layout_2d_box_index(pic2d%layout2d,rank,boxindex)
    global_sz=(/pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2/)
    gxmin=pic2d%para2d%gxmin
    delta(1)=pic2d%para2d%m_x1%delta_eta
    delta(2)=pic2d%para2d%m_x2%delta_eta

    call gather_field_to_rootprocess_per_per(rootdata%field,boxfield,rank,size,boxindex, &
         pic2d%para2d%numproc,pic2d%layout2d)    

    if(rank==0) then
       do j=1,global_sz(2)
          x(2)=gxmin(2)+real(j-1,8)*delta(2)     
          do i=1,global_sz(1)
             if(j==1.and.i==1) then
                write(fileitem,*)
                write(fileitem,*)
             else if(j.ge.2.and.i==1) then
                write(fileitem, *)
             end if
             x(1)=gxmin(1)+real(i-1,8)*delta(1)
             globind=(j-1)*global_sz(1)+i
             write(fileitem,"(3E19.10)") x(2), rootdata%field(globind),x(1)           
          end do
       end do

    end if
 
    end subroutine para_write_field_file_2d

   
    subroutine allmu_para_write_field_file_2d(boxfield,fileitem,muind,rootdata,pic2d)
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d      
      real8, dimension(:,:,:), pointer, intent(in) :: boxfield
      int4, intent(in) :: fileitem,muind
      class(root_precompute_data), pointer, intent(in) :: rootdata
      real8, dimension(:,:), pointer :: buff
      int4 :: num1,num2

      num1=size(boxfield,2)
      num2=size(boxfield,3)
      allocate(buff(num1,num2))
      
      buff=boxfield(muind,:,:)

      call para_write_field_file_2d(buff,fileitem,rootdata,pic2d) 
   
      deallocate(buff)
    end subroutine

    subroutine para_write_orbit_file_2d(tp_ful2d_head,numleft,numgr,fileitem,pic2d,iter_num)
      class(tp_ful2d_node), pointer,intent(in) :: tp_ful2d_head
      int4, intent(in) :: fileitem,numleft,numgr,iter_num
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      int4, dimension(:),pointer :: rcounts,rdispls 
      real8, dimension(:), pointer :: outarray,outarrayone
      int4 :: numtot,size,rank,comm,i
      real8 :: vperp,rholength,gc(2)      

      comm=pic2d%layout2d%collective%comm
      size=pic2d%layout2d%collective%size
      rank=pic2d%layout2d%collective%rank
      allocate(rcounts(0:size-1),rdispls(0:size-1))
      call prepare_recv_for_gatherv(numtot,rcounts,rdispls,numleft,numgr,size,comm,rank) 
      if(rank==0) then
        allocate(outarray(1:(numtot/numgr)*(numgr-1)))
        allocate(outarrayone(1:(numtot/numgr)*(numgr+1)))
      else
        allocate(outarray(0:0))
      end if

      call tp_gather_coords_to_rootprocess(outarray,tp_ful2d_head,numleft,numtot,numgr,rcounts, &
                                               rdispls,pic2d)
!!!!!! HERE, numger=5
      if(rank==0) then
        do i=0,numtot/numgr-1
          outarrayone((numgr+1)*i+1:(numgr+1)*i+4)=outarray((numgr-1)*i+1:(numgr-1)*i+4)
          vperp=sqrt(outarray((numgr-1)*i+3)**2+outarray((numgr-1)*i+4)**2)
          rholength=vperp
          outarrayone((numgr+1)*i+5)=outarray((numgr-1)*i+1)+rholength*outarray((numgr-1)*i+4)/vperp
          outarrayone((numgr+1)*i+6)=outarray((numgr-1)*i+2)-rholength*outarray((numgr-1)*i+3)/vperp          
        end do 

         write(fileitem, *) iter_num,outarrayone(:)
!         print*, outarrayone(1:5),outarrayone((numtot/numgr)*(numgr+1))
         deallocate(outarrayone)
      end if

      deallocate(rcounts,rdispls)
        deallocate(outarray)
    end subroutine para_write_orbit_file_2d 


    subroutine para_write_orbit_file_2d_gy(tp_gy2dallmu_head,numleft,numgr,muind,fileitem,pic2d,iter_num)
      class(tp_gy2dallmu_node), dimension(:),pointer,intent(in) :: tp_gy2dallmu_head
      int4, intent(in) :: fileitem,numleft,numgr,iter_num
      int4, intent(in) :: muind
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      int4, dimension(:),pointer :: rcounts,rdispls 
      real8, dimension(:), pointer :: outarray
      class(tp_gy2dallmu_node), dimension(:), pointer :: tpgy2dmutmp
      class(tp_gy2d_node), pointer :: tpgy2d_head, tpgy2dtmp
      int4 :: numtot,size,rank,comm,mu_num
      int4 :: mu_ind
      int4 :: i
    
      mu_num=pic2d%para2d%mu_num  
      comm=pic2d%layout2d%collective%comm
      size=pic2d%layout2d%collective%size
      rank=pic2d%layout2d%collective%rank
      allocate(rcounts(0:size-1),rdispls(0:size-1))
      allocate(tpgy2dmutmp(mu_num))
      allocate(tpgy2d_head)
      do i=1,mu_num
        tpgy2dmutmp(i)%ptr=>tp_gy2dallmu_head(i)%ptr
      end do

      call prepare_recv_for_gatherv(numtot,rcounts,rdispls,numleft,numgr,size,comm,rank) 
      if(rank==0) then
        allocate(outarray(1:(numtot/numgr)*(numgr-1)))
      else
        allocate(outarray(0:0))
      end if

      tpgy2dtmp=>tpgy2d_head
      do while(associated(tpgy2dmutmp(muind)%ptr))
        if(.not.associated(tpgy2dmutmp(muind)%ptr%next)) then
           exit
        else
          tpgy2dtmp%coords(1:3)=tpgy2dmutmp(muind)%ptr%coords(1:3)        
          allocate(tpgy2dtmp%next)
          tpgy2dtmp=>tpgy2dtmp%next
          tpgy2dmutmp(muind)%ptr=>tpgy2dmutmp(muind)%ptr%next
        end if
      end do

     call tp_gather_coords_to_rootprocess_gy(outarray,tpgy2d_head,numleft,numtot,numgr,rcounts, &
                                               rdispls,pic2d)

      if(rank==0) then
         write(fileitem, *) iter_num,outarray(:)
      end if

      deallocate(rcounts,rdispls)
        deallocate(outarray)
    end subroutine para_write_orbit_file_2d_gy 

    subroutine para_write_orbit_file_2d_gy_allmu(tp_gy2dallmu_head,numleft,numgr,fileitem,muind, & 
               pic2d,iter_num)
      class(tp_gy2dallmu_node), dimension(:),pointer,intent(in) :: tp_gy2dallmu_head
      int4, intent(in) :: fileitem,numleft,numgr,iter_num
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      int4, intent(in) :: muind
      int4, dimension(:),pointer :: rcounts,rdispls 
      real8, dimension(:), pointer :: outarray
      int4 :: numtot,size,rank,comm,mu_num
      int4 :: i      
      class(tp_gy2dallmu_node), dimension(:), pointer :: tpgy2dmutmp
      class(tp_gy2d_node), pointer :: tp_gy2d_head, tpgy2dtmp
 
      comm=pic2d%layout2d%collective%comm
      size=pic2d%layout2d%collective%size
      rank=pic2d%layout2d%collective%rank
      mu_num=pic2d%para2d%mu_num
      allocate(tpgy2dmutmp(1:mu_num))
      do i=1,mu_num
        allocate(tpgy2dmutmp(i)%ptr)
        tpgy2dmutmp(i)%ptr=>tp_gy2dallmu_head(i)%ptr
      end do
      allocate(rcounts(0:size-1),rdispls(0:size-1))
      call prepare_recv_for_gatherv(numtot,rcounts,rdispls,numleft,numgr,size,comm,rank) 
      if(rank==0) then
        allocate(outarray(1:(numtot/numgr)*(numgr-1)))
      else
        allocate(outarray(0:0))
      end if
   
     allocate(tp_gy2d_head)
     tpgy2dtmp=>tp_gy2d_head
      do while(associated(tpgy2dtmp)) 
         if(.not.associated(tpgy2dtmp%next)) then
           exit
         else
           tpgy2dtmp%coords(1:3)=tpgy2dmutmp(muind)%ptr%coords(1:3)
           allocate(tpgy2dtmp%next)
           tpgy2dtmp=>tpgy2dtmp%next
           tpgy2dmutmp(muind)%ptr=>tpgy2dmutmp(muind)%ptr%next
         end if
      end do
     
      call tp_gather_coords_to_rootprocess_gy(outarray,tp_gy2d_head,numleft,numtot,numgr,rcounts, &
                                               rdispls,pic2d)

      if(rank==0) then
         write(fileitem, *) iter_num,outarray(:)
      end if

      deallocate(rcounts,rdispls)
        deallocate(outarray)
      deallocate(tp_gy2d_head)
    end subroutine para_write_orbit_file_2d_gy_allmu 



    subroutine tp_gather_coords_to_rootprocess(outarray,tp_ful2d_head,numleft,numtot,numgr,rcounts, &
                                               rdispls,pic2d)
      class(tp_ful2d_node),pointer, intent(in) :: tp_ful2d_head
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      int4,intent(in) :: numleft, numgr,numtot    
      int4, dimension(:), pointer,intent(in) :: rdispls,rcounts
      real8, dimension(:) :: outarray
      int4 :: root,sizeone,comm,rank
      real8, dimension(:), pointer :: rbuf, sbuf
      class(tp_ful2d_node), pointer :: tptmp 
      int4 :: h,i,ierr

      sizeone=pic2d%layout2d%collective%size
      comm=pic2d%layout2d%collective%comm
      rank=pic2d%layout2d%collective%rank
      allocate(rbuf(0:numtot-1))
      root=0 
        tptmp=>tp_ful2d_head 
    if(numleft==0) then
        allocate(sbuf(0:0))
        call mpi_gatherv(sbuf,numleft*numgr,mpi_double_precision,rbuf,rcounts, &
                       rdispls,mpi_double_precision,root,comm,ierr)     
      nullify(tptmp) 
    else    
        allocate(sbuf(0:numgr*numleft-1))
        h=0
        do while(associated(tptmp)) 
           if(.not.associated(tptmp%next)) then
             exit
           else
             sbuf(numgr*h:numgr*h+numgr-2)=tptmp%coords(1:4)
             sbuf(numgr*h+numgr-1)=real(tptmp%tp,8)
             h=h+1
             tptmp=>tptmp%next
           end if
        end do
        call mpi_gatherv(sbuf,numleft*numgr,mpi_double_precision,rbuf,rcounts, &
                       rdispls,mpi_double_precision,root,comm,ierr)
        nullify(tptmp)
    end if  

      if(rank==0) then
        do i=0,numtot/numgr-1
           h=NINT(rbuf(i*numgr+numgr-1))
           outarray((numgr-1)*(h-1)+1:(numgr-1)*(h-1)+numgr-1)=rbuf(i*numgr:i*numgr+numgr-2) 
        end do
      end if

      deallocate(rbuf,sbuf)
  end subroutine tp_gather_coords_to_rootprocess

  subroutine tp_gather_coords_to_rootprocess_gy(outarray,tp_gy2d_head,numleft,numtot,numgr,rcounts, &
                                               rdispls,pic2d)
      class(tp_gy2d_node),pointer, intent(in) :: tp_gy2d_head
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      int4,intent(in) :: numleft, numgr,numtot    
      int4, dimension(:), pointer,intent(in) :: rdispls,rcounts
      real8, dimension(:) :: outarray
      int4 :: root,size,comm,rank
      real8, dimension(:), pointer :: rbuf, sbuf
      class(tp_gy2d_node), pointer :: tptmp 
      int4 :: h,i,ierr

      size=pic2d%layout2d%collective%size
      comm=pic2d%layout2d%collective%comm
      rank=pic2d%layout2d%collective%rank
      allocate(rbuf(0:numtot-1))
      root=0 
        tptmp=>tp_gy2d_head 
    if(numleft==0) then
        allocate(sbuf(0:0))
        call mpi_gatherv(sbuf,numleft*numgr,mpi_double_precision,rbuf,rcounts, &
                       rdispls,mpi_double_precision,root,comm,ierr)     
     nullify(tptmp) 
    else    
        allocate(sbuf(0:numgr*numleft-1))
        h=0
        do while(associated(tptmp)) 
           if(.not.associated(tptmp%next)) then
             exit
           else
             sbuf(numgr*h:numgr*h+numgr-2)=tptmp%coords(1:3)
             sbuf(numgr*h+numgr-1)=real(tptmp%tp,8)
             h=h+1
             tptmp=>tptmp%next
           end if
        end do
        call mpi_gatherv(sbuf,numleft*numgr,mpi_double_precision,rbuf,rcounts, &
                       rdispls,mpi_double_precision,root,comm,ierr)
        nullify(tptmp)
    end if   
      if(rank==0) then
         do i=0,numtot/numgr-1
           h=NINT(rbuf(i*numgr+numgr-1))
           outarray((numgr-1)*(h-1)+1:(numgr-1)*(h-1)+numgr-1)=rbuf(i*numgr:i*numgr+numgr-2) 

        end do
      end if

      deallocate(rbuf,sbuf)
  end subroutine tp_gather_coords_to_rootprocess_gy


        

end module para_write_file
