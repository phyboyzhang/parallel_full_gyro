module para_gyroaverage_2d_one
#include "work_precision.h"
use cartesian_mesh, only: cartesian_mesh_1d
  
use gyroaverage_utilities, only: s_compute_shape_circle
  
use spline_module, only: s_splcoefper1d0old,&
                         s_compute_splines_coefs_matrix_per_1d, &
                         s_compute_splines_coefs_matrix_nat_1d_new, &
                         s_localize_per_per, &
                         s_contribution_spl, &
                         splcoefper1d0old, &
                         s_localize_nat_per_new, &
                         compute_D_spl2D_per_per_noblock, &
                         compute_D_spl2D_nat_per_noblock, &
                         s_localize_new      

use m_para_spline, only: para_compute_spl2d_point_per_per_weight 
use piclayout, only: root_precompute_data
use paradata_type, only: pic_para_total2d_base
use gyroaverage_2d_base, only: gyropoint_node, &
                              gyropoint, &
                              gyroweight
                       

use m_parautilities, only: gather_field_to_rootprocess_per_per
use paradata_utilities, only: compute_process_of_point_per_per, &
                              coordinate_pointoutbound_per, &
                              coordinates_pointoutbound_per_per, &
                              globalind_from_localind_2d, &
                              startind_of_process, &
                              point_location_for_cubic_spline_2d
                              

implicit none
include "mpif.h"

  public :: sort_quadraturepoint_among_process, &
            para_compute_gyroaverage_stiff_matrix, &
            para_compute_gyroaverage_mesh_field, &
            muarray_euler_maclaurin_choice
!            store_data_on_rootprocess,  &
!            precompute_doublegyroaverage_matrix


contains

!  function find_maxbf_box(pgyro_2d,pic2d)
!    class(para_gyroaverage_2d), intent(inout) :: paragyro2d
!    class(pic_para_ful2d_base), intent(in) :: pic2d
!    
!    
!    
!  end function

!  subroutine initialize_gyro_boxextend(paragyro,mu,pic2d,geometry,boundary)  
!    type(para_gyroaverage_2d_plan), pointer,intent(inout) :: paragyro
!    class(pic_para_ful2d_base), pointer,intent(in) :: pic2d
!    character(len=*), intent(in) :: geometry
!    int4 :: rank,index(4),dimsize(2),coords(2),length(2)
!    real8 :: rho
!    logical :: periods(0:1)
!    int4 :: reorder=1,numproc1,numproc2,commold,comm=2
!    int4 :: nb(2)
!    allocate(paragyro)
!    numproc1=pic2d%para2d%numproc1
!    numproc2=pic2d%para2d%numproc2
!    dimsize(1)=numproc1
!    dimsize(2)=numproc2
!    paragyro%mu=mu
!    rho=sqrt(2._f64*mu)
!    rank=pic2d%layout2d%collective%rank
!    commold=pic2d%layout2d%collective%comm
!    length(1)=(pic2d%para2d%gxmax(1)-pic2d%para2d%gxmin(1))
!    length(2)=(pic2d%para2d%gxmax(2)-pic2d%para2d%gxmin(2)) 
!    nb(1)=floor(rho/length(1))+1
!    nb(2)=floor(rho/length(2))+1
!    if(real(n-1,8)*length(1).lt.2._f64*pic2d%para2d%m_x1%delta_eta) then
!      nb(1)=nb(1)+1
!    end if
!    if(real(n-1,8)*length(2).lt.2._f64*pic2d%para2d%m_x2%delta_eta) then
!      nb(2)=nb(2)+1
!    end if
!        
!    call get_layout_2d_box_index(layout2d,rank,index)
!    
!    select case(geometry)
!    case("cartesian")
!      select case(boundary)
!        case("double_per")
!          periods(0)=.true.
!          periods(1)=.true.          
!          call MPI_Cart_create(commold,2,dimsize,periods,reorder,comm,ierr)
!          allocate(boxextends((2*nb(1)+1)*pic2d%para2d%m_x1%num_cells+1,(2*nb(2)+1)*pic2d%para2d%m_x2%num_cells+1))         
!           
!        case("nat_per")
!
!        case default
!          stop
!      end select
!
!    case("polar")
!
!    case default
!      stop 
!    end select
!
!  end subroutine

!  function para_compute_gyroaverage_2d_point_double_per(x1,mu,gyro,pic2d,geometry,boundary) result(field_value)
!       type(gyroaverage_2d_plan), pointer, intent(in) ::  gyro
!       class(pic_para_ful2d_base),intent(in) :: pic2d
!       real8, intent(in) :: mu
!       real8, dimension(2),intent(in) :: x1
!       character(len=*),intent(in), optional :: gyroaverage,boundary
!!       character(len=*),intent(in) :: boundary, geometry
!       real8 :: field_value
!       real8 :: rho
!       int4 :: i,j,k,row(4),rank,NC(2),numproc(2),rankpoint
!       real8 :: points(3,gyro%N_points) 
!       real8 :: x(2),eta_min(2),eta_max(2),eta_star(2)
!       int4 :: ii(2), Nc(2), flag, ind(2), ell_1, ell_2, num=0
!       real8 :: val(-1:2,-1:2)
!       type(gyropoint), pointer :: pointhead,curpoint
!
!       call s_compute_shape_circle(points,gyro%N_points)
!    
!       Nc(1) = m_x1%nodes
!       Nc(2) = m_x2%nodes
!       rho = sqrt(2.0d0*mu)
!       eta_min(1)=m_x1%eta_min
!       eta_min(2)=m_x2%eta_min
!       eta_max(1)=m_x1%eta_max
!       eta_max(2)=m_x2%eta_max
!       row=pic2d%para2d%row
!   !    rank=pic2d%layout2d%collective%rank
!       numproc(1)=pic2d%para2d%numproc1
!       numproc(2)=pic2d%para2d%numproc2
!   
!       field_value=0._f64 
!       if(present(gyroaverage)) then
!          if(abs(x1(1)-gxmin(1)).lt.rho.or. abs(x1(1)-gxmax(1)).lt.rho.or.&
!             abs(x1(2)-gxmin(2)).lt.rho.or. abs(x1(2)-gxmax(2)).lt.rho)    then
!             field_value=para_compute_spl2d_field_point_per_per(x,eta_min,eta_max,row,rank,NC,numproc, &
!               pic2d%field2d%epwg_rw,pic2d%field2d%epwg_re,pic2d%field2d%epwg_rn,pic2d%field2d%epwg_rs, pic2d%field2d%epwg_rsw,&
!               pic2d%field2d%rse,pic2d%field2d%rne,pic2d%field2d%rnw,pic2d%field2d%rnw,pic2d%layout2d%collective%comm )
!          else 
!             allocate(pointhead) 
!             curpoint => pointhead           
!             do k=1,gyro%N_points
!             flag=1
!             x(1) = x1(1)+rho*points(1,k)
!             x(2) = x1(2)+rho*points(2,k)
!             do i=1, numproc(1)*numproc(2)
!               !!! obtain the rank of the box where the point locates
!               if(pic2d%para2d%gxmin(i,1).le.x(1).and.pic2d%para2d%gxmax(i,1).ge.x(1).and.pic2d%para2d%gxmin(i,2).le.x(2) &
!                 .and.pic2d%para2d%gxmax(i,2).ge.x(2)) then
!                 rankpoint=i
!               else 
!                 print*, "# the kth point locates out of the simulated domain."
!                 stop
!               end if
!             end do 
!             curpoint%rankgrid=pic2d%layout2d%collective%rank
!             curpoint%rankpoint=rankpoint
!             curpoint%gridind=(/x1(1),x1(2)/)
!             curpoint%xpoint=(/x(1),x1(2)/)
!
!          call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,Nc,flag)
!    
!          if(flag==1) then
!             call s_contribution_spl(eta_star,val)
!             val=val*points(3,k)  
!             do ell_2=-1,2
!                ind(2)=modulo(ii(2)+ell_2,Nc(2))
!                  do ell_1=-1,2           
!                    ind(1)=modulo(ii(1)+ell_1,Nc(1))   ! this is change to
!                    pafter=pafter+val(ell_1,ell_2)*(gyro%ep_weight_gyro(ind(1)+1,ind(2)+1))
!                  enddo
!             enddo
!          end if           
!        enddo
!
!       end if
!
!   else 
!     
!         if(abs(x1(1)-m_x1%eta_min).lt.rho.or. abs(x1(1)-m_x1%eta_max).lt.rho.or.&
!!! drift kineitc is implemented at the  boundary 
!            abs(x1(2)-m_x2%eta_min).lt.rho.or. abs(x1(2)-m_x2%eta_max).lt.rho)    then
!           pafter=compute_spl2d_field_point_per_per(x1, &
!                       field_2d%ep_weight_init,  &
!                       (/m_x1%eta_min,m_x2%eta_min/),    &
!                       (/m_x1%eta_max,m_x2%eta_max/),    &
!                       (/m_x1%delta_eta,m_x2%delta_eta/),  &
!                       (/m_x1%nodes,m_x2%nodes/))
!         else
!
!            do k=1,gyro%N_points
!              flag=1
!              x(1) = x1(1)+rho*points(1,k)
!              x(2) = x1(2)+rho*points(2,k)
!              call s_localize_per_per(x,eta_min,eta_max,ii,eta_star,Nc,flag)
!    
!              if(flag==1) then
!                 call s_contribution_spl(eta_star,val)
!                 val=val*points(3,k)  
!                 do ell_2=-1,2
!                   ind(2)=modulo(ii(2)+ell_2,Nc(2))
!                   do ell_1=-1,2           
!                     ind(1)=modulo(ii(1)+ell_1,Nc(1))   ! this is change to
!                     pafter=pafter+val(ell_1,ell_2)*(field_2d%ep_weight_init(ind(1)+1,ind(2)+1))
!                   enddo
!                 enddo
!              end if
!            enddo
!
!        end if
!
!    end if
!    
! !   print*, "pafter=",pafter, val(0,0:1)
!  end function para_compute_gyroaverage_2d_point_double_per


!!!! Sort the quadraturepoint on the current process accooridng to the processes they locate
  subroutine sort_quadraturepoint_among_process(mu,rank,pointhead,num_p,N_points,pic2d)
     class(pic_para_total2d_base), intent(in) :: pic2d
!     character(len=*), intent(in) :: geometry,boundary
  !   real8, dimension(:), intent(in) :: mu
     real8, intent(in) :: mu
     int4, intent(in) :: N_points,rank
     class(gyropoint_node),dimension(:), pointer,intent(inout) :: pointhead
     class(gyropoint_node),dimension(:), pointer :: curpoint
     int4, dimension(:),pointer, intent(inout) :: num_p
     real(8), dimension(:,:), allocatable :: points
     real8 :: rho,x1(2),x(2)
     real8 :: gxmin(2),gxmax(2)
     int4 :: num_mu,i,j,k,h,m,rankpoint,size
     int4 :: NC(2),flag     

  !   num_mu=size(mu)
     Nc(1) = pic2d%para2d%m_x1%nodes
     Nc(2) = pic2d%para2d%m_x2%nodes
     size=pic2d%layout2d%collective%size  
     allocate(points(N_points,3))
     allocate(curpoint(0:size-1))
     rho=sqrt(2.0*mu)
     do i=0, size-1 
       if(.not.associated(pointhead(i)%ptr)) then
          print*, "pointhead(i)%ptr is not associated."
          stop
       end if
       curpoint(i)%ptr=>pointhead(i)%ptr
     end do

    
  call s_compute_shape_circle(points,N_points)
!  h=0
  select case(pic2d%para2d%geometry)
    case ("cartesian")  
    select case(pic2d%para2d%boundary)
     case("nat_per")    
  !   do m=1,num_mu
     do i=1, Nc(1)
       x1(1)=pic2d%para2d%gboxmin(rank,1)+real(i-1,8)*pic2d%para2d%m_x1%delta_eta
       do j=1, Nc(2) 
         x1(2)=pic2d%para2d%gboxmin(rank,2)+real(j-1,8)*pic2d%para2d%m_x2%delta_eta
             do k=1,N_points
               flag=1
               x(1) = x1(1)+rho*points(1,k)
               x(2) = x1(2)+rho*points(2,k)
               x(2)=coordinate_pointoutbound_per(x(2),gxmin(2),gxmax(2))
               !!! For the dimension where the natural boundary condition is adopted, the point doesn't change.       
               if(x(1).lt.gxmin(1).or.x(1).gt.gxmax(1)) then
                  x(1)=x1(1)
               end if

               rankpoint= compute_process_of_point_per_per(x,pic2d%para2d%numproc,pic2d%para2d%gxmin, &
                          pic2d%para2d%gxmax,pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
               num_p(rankpoint)=num_p(rankpoint)+1   
               if(.not.associated(curpoint(rankpoint)%ptr)) then
                  stop
               else 
               curpoint(rankpoint)%ptr%gridrank=rank
                curpoint(rankpoint)%ptr%gridind(i)=i
                curpoint(rankpoint)%ptr%gridind(j)=j
                curpoint(rankpoint)%ptr%xpoint=(/x(1),x(2)/)
                allocate(curpoint(rankpoint)%ptr%next)
                curpoint(rankpoint)%ptr=>curpoint(rankpoint)%ptr%next   
               end if 
            end do !!! end k 
!          end if

         end do  !!! end j
       end do  !!! end i
 !    end do  !!! end m

     case ("double_per")
 !    do m=1,num_mu
     do i=1, Nc(1)
         x1(1)=pic2d%para2d%gboxmin(rank,1)+real(i-1,8)*pic2d%para2d%m_x1%delta_eta 
         do j=1, Nc(2) 
            x1(2)=pic2d%para2d%gboxmin(rank,2)+real(j-1,8)*pic2d%para2d%m_x2%delta_eta
            do k=1,N_points
               flag=1
               x(1) = x1(1)+rho*points(1,k)
               x(2) = x1(2)+rho*points(2,k)
               call coordinates_pointoutbound_per_per(x,pic2d%para2d%gxmin,pic2d%para2d%gxmax)
               rankpoint=compute_process_of_point_per_per(x,pic2d%para2d%numproc, &
                   pic2d%para2d%gxmin,pic2d%para2d%gxmax,pic2d%para2d%gboxmin,pic2d%para2d%gboxmax)
               num_p(rankpoint)=num_p(rankpoint)+1 
               curpoint(rankpoint)%ptr%gridrank=rank
               curpoint(rankpoint)%ptr%gridind(1)=i
               curpoint(rankpoint)%ptr%gridind(2)=j
               curpoint(rankpoint)%ptr%xpoint(:)=(/x(1),x(2)/)

!if(rank==3) then
!h=h+1
!print*, "rankpoint=",rankpoint,curpoint(rankpoint)%ptr%xpoint, h 
!end if 
               allocate(curpoint(rankpoint)%ptr%next)
                curpoint(rankpoint)%ptr=>curpoint(rankpoint)%ptr%next   
            end do !!! end k 
        end do  !!! end j
      end do  !!! end i
 !    end do  !!! end m
     case default
       stop
     end select
   
    case ("polar")
      stop
    case default
      stop
    
    end select

!   do i=0, size-1
!     deallocate(curpoint(i)) 
!     nullify(curpoint(i)%ptr)
!   end do  
   
   end subroutine sort_quadraturepoint_among_process

   !!!! First send the quadrature points to the processes they locates, then send all weight to the root process, 
   !!!! for the solving of the stiff matrix
   subroutine  para_compute_gyroaverage_stiff_matrix(pointhead,num_p,mu,mu_num,N_points,pic2d,rootdata) 
    class(pic_para_total2d_base), pointer,intent(in) :: pic2d
    class(root_precompute_data),  pointer, intent(inout) :: rootdata
!    real8,dimension(:), intent(in) :: mu
    real8, intent(in) :: mu
    int4, intent(in) :: N_points,mu_num
    class(gyropoint_node),dimension(:), pointer,intent(inout) :: pointhead
    class(gyropoint_node),dimension(:), pointer :: curpoint
    int4, dimension(:), pointer,intent(inout) :: num_p
 !   int4, intent(in) :: commone   !! The original communication world
    int4, dimension(:), pointer :: rbufone  !! store the number of quadrature
 !!!! sent from all the processes
    int4, dimension(:), pointer :: rcounts,scounts,sdispls,rdispls
    real8,dimension(:), pointer :: sbuf,rbuf
    int4, dimension(:), pointer :: sdispls2nd,rdispls2nd
    real8, dimension(:),pointer :: sbuf2nd,rbuf2nd
    int4, dimension(:), pointer :: rbuf3rd
    int4 :: ierr,numout,numsend,numrecv,comm
    int4 :: nearind(2),myrank,root,size
    real8 :: weight(-1:2,-1:2),val(-1:2,-1:2)
    real8 :: x(2),eta_min(2),eta_max(2),eta_star(2), rho
    int4 :: ii(2),iii(2), Nc(2), flag, ind(2), ell_1,  &
            ell_2,row(4),numproc(2),gridrank
    int4 :: contrindex(4,2)  ! store the global indexes of the comtribution points
    int4 :: stifmatind(2),globind(2),startind(2),nx1,nx2,gridind(2),startgridind(2)
    int4 :: num, rankpoint,localcoords(2)
    int4 :: i,j,h,l,m    

    root=0
    myrank=pic2d%layout2d%collective%rank
    size=pic2d%para2d%numproc(1)*pic2d%para2d%numproc(2)
    comm=pic2d%layout2d%collective%comm
    allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))  
    allocate(rbuf3rd(0:size-1),rbufone(0:size-1)) 
    allocate(curpoint(0:size-1))
    do i=0,size-1
      if(.not.associated(pointhead(i)%ptr)) then
        stop
      else
         curpoint(i)%ptr=>pointhead(i)%ptr 
      end if
    end do

!    call sort_quadraturepoint_among_process(mu,myrank,pointhead,num_p,N_points,pic2d)
!if(myrank==3) then
!do i=0,size-1
!print*,"i=",i, "pointhead(i)%ptr=",pointhead(i)%ptr%gridrank,  &
!pointhead(i)%ptr%gridind(1), pointhead(i)%ptr%gridind(2)
!enddo
!end if 


   call mpi_alltoall(num_p,1,mpi_integer,rbufone,1,mpi_integer,comm,ierr)
    numsend=0
    numrecv=0
    do i=0,size-1
      numsend=numsend+num_p(i)  
      numrecv=numrecv+rbufone(i)
      scounts(i)=num_p(i)*5
      rcounts(i)=rbufone(i)*5
      if(i==0) then
        sdispls(0)=0
        rdispls(0)=0
      else
        sdispls(i)=sdispls(i-1)+scounts(i-1)
        rdispls(i)=rdispls(i-1)+rcounts(i-1)
      end if 
    end do    
    allocate(sbuf(0:numsend*5-1))
    allocate(rbuf(0:numrecv*5-1))
    h=0
!    m=0
    do i=0, size-1
       do while(associated(curpoint(i)%ptr).and.associated(curpoint(i)%ptr%next))

!if(myrank==3) then 
!m=m+1      
!print*, "i=",i,curpoint(i)%ptr%xpoint, curpoint(i)%ptr%gridind(1), &
!            curpoint(i)%ptr%gridind(2),"m=",m
!end if
          sbuf(5*h)=real(curpoint(i)%ptr%gridrank,8)
          sbuf(5*h+1:5*h+2)=curpoint(i)%ptr%xpoint
          sbuf(5*h+3)=real(curpoint(i)%ptr%gridind(1),8)
          sbuf(5*h+4)=real(curpoint(i)%ptr%gridind(2),8)
          curpoint(i)%ptr=>curpoint(i)%ptr%next
          h=h+1
!          if(.not.associated(curpoint(i)%ptr%next)) then  !!! This is the last element
!             goto 50
!          end if  
        end do
50    end do

   call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcounts,rdispls,mpi_double_precision,comm,ierr) 

   !!!! now prepare the sendbuf and recvbuf for mpi_gatherv operation     
   !!!! The 33 elements are: for each contribution point, the first is the globalind(1), then 16 pairs of (val(i,j),globalind(2))
   !!!! The globalind(1) is for the grid points which the gyroaverage is computed for. The globalind(2) is for the quadraturepoints.
   numsend=0
   do i=0,size-1
      numsend=numsend+rbufone(i)
   end do 
   allocate(sbuf2nd(0:numsend*33-1), stat=ierr)

   root=0   
   call mpi_gather(numsend,1,mpi_integer,rbuf3rd,1,mpi_integer,root,comm,ierr)

   numrecv=0
   do i=0,size-1
      numrecv=numrecv+rbuf3rd(i)
      rcounts(i)=rbuf3rd(i)*33
      if(i==0) then
        rdispls(i)=0
      else
        rdispls(i)=rdispls(i-1)+rcounts(i-1)
      end if
   end do
   allocate(rbuf2nd(0:33*numrecv-1),stat=ierr)
   numout=numrecv

   !!!!! 14 elements: ptrank(which is the current rank), nearind(2), gridind(2), 9 weight(-1:1,-1:1)
!!                             O
!!                             *
!!                             O
!!                             *
!!                   ******O***O***O**O***          *
!!                             * nearind(2)
!!                             O
!!                             *


       Nc(1) = pic2d%para2d%m_x1%nodes
       Nc(2) = pic2d%para2d%m_x2%nodes
       rho = sqrt(2.0d0*mu)
       eta_min(1)=pic2d%para2d%m_x1%eta_min
       eta_min(2)=pic2d%para2d%m_x2%eta_min
       eta_max(1)=pic2d%para2d%m_x1%eta_max
       eta_max(2)=pic2d%para2d%m_x2%eta_max
       row=pic2d%para2d%row
       numproc(1)=pic2d%para2d%numproc(1)
       numproc(2)=pic2d%para2d%numproc(2)  

!!! compute the weight of the contribution points and store them in the global FEM matrix
  
       startind=startind_of_process(myrank,numproc,pic2d%layout2d)
   
        h=0
          do i=0,size-1
            num=0
            do while(num.lt.rbufone(i))
              x=rbuf(5*h+1:5*h+2)

        !!! compute the globalindex of the grid points where the gyroaveage happens   
              gridind(1)=NINT(rbuf(5*h+3))
              gridind(2)=NINT(rbuf(5*h+4))
              gridrank=NINT(rbuf(5*h))
              globind=globalind_from_localind_2d(gridind,numproc,gridrank,pic2d%layout2d,pic2d%para2d%boundary)
              stifmatind(1)=(globind(2)-1)*pic2d%layout2d%global_sz1+globind(1)

        !!!!! here -1 should be paid attention
              sbuf2nd(33*h)=real(stifmatind(1),8)

        !!! compute the  globalindex of contribution points
        !!! obtain the local index of the nearest points which is ii here. 
              select case(pic2d%para2d%geometry)
                case("cartesian")
                select case(pic2d%para2d%boundary)
                  case("double_per")
                    call s_localize_new(x,eta_min,eta_max,ii,eta_star,Nc-(/1,1/),flag)  
                  case("nat_per")
                    call s_localize_new(x,eta_min,eta_max,ii,eta_star,Nc-(/1,1/),flag)
                  case default
                    stop
                  end select

              case("polar")
                call s_localize_new(x,eta_min,eta_max,ii,eta_star,Nc,flag)
              case default
                stop
            end select

            call s_contribution_spl(eta_star,val)
            val=val/real(N_points,8)
            do m=-1,2
              iii(1)=ii(1)+m
              do l=-1,2
                iii(2)=ii(2)+l
                call point_location_for_cubic_spline_2d(rankpoint,localcoords,myrank,numproc,  &
                        iii,pic2d%para2d,pic2d%layout2d)
                globind=globalind_from_localind_2d(localcoords,numproc,rankpoint,pic2d%layout2d,pic2d%para2d%boundary)
                stifmatind(2)=(globind(2)-1)*pic2d%layout2d%global_sz1+globind(1)
                sbuf2nd(33*h+4*(m+1)+l+2)=val(m,l)           
                sbuf2nd(33*h+4*(m+1)+l+2+16)=real(stifmatind(2),8)
              end do
            end do
            num=num+1
            h=h+1
          end do
        end do
        root=0
        call MPI_GATHERV(sbuf2nd,33*numsend,MPI_DOUBLE_PRECISION,rbuf2nd,rcounts,rdispls,MPI_DOUBLE_PRECISION,root,comm,ierr) 

!!!! Now, get the global stiff matrix on the root process.    

        if(myrank==root) then
 !     allocate(rootdata%ACONTRI(pic2d%layout2d%global_sz2*pic2d%layout2d%global_sz1, &
 !              pic2d%layout2d%global_sz2*pic2d%layout2d%global_sz1),stat=ierr)
        rootdata%ACONTRI=0._F64
        do h=0, numout-1
        do m=-1,2
            do l=-1,2
          !!!! Tell whether the the point which h denotes for is counted or not to avoid the repeated counting of points at 
          !!!! the boundary of each box.
               IF(rootdata%ACONTRI(NINT(rbuf2nd(33*h)),NINT(rbuf2nd(33*h+4*(m+1)+l+2+16))).ne.0.or. &
                  rootdata%ACONTRI(NINT(rbuf2nd(33*h)),NINT(rbuf2nd(33*h+4*(m+1+1)+l+1+2+16))).ne.0.or. &
                  rootdata%ACONTRI(NINT(rbuf2nd(33*h)),NINT(rbuf2nd(33*h+4*(m+1+2)+l+2+2+16))).ne.0  ) then
                  goto 100 
               else

                  rootdata%ACONTRI(NINT(rbuf2nd(33*h)),NINT(rbuf2nd(33*h+4*(m+1)+l+2+16))) &
                  =rootdata%ACONTRI(NINT(rbuf2nd(33*h)),NINT(rbuf2nd(33*h+4*(m+1)+l+2+16))) &
                  +rbuf(33*h+4*(m+1)+l+2)
               end if
            end do
         end do
100    end do
    end if  
       
  end subroutine para_compute_gyroaverage_stiff_matrix

  !!! given the original field on the mesh,this subroutine is used to compute
  !the gyroaverage potential on the mesh
  subroutine  para_compute_gyroaverage_mesh_field(num_p,mu,mu_num,pic2d) 
        class(pic_para_total2d_base), pointer,intent(inout) :: pic2d
        !  real8,dimension(:), intent(in) :: mu
        real8, intent(in) :: mu
        class(gyropoint_node),dimension(:), pointer :: pointhead,curpoint
        int4, dimension(:), pointer,intent(inout) :: num_p
        int4, intent(in) :: mu_num
        int4, dimension(:), pointer :: rcounts,scounts,sdispls,rdispls,rcountsone, &
        scountstwo,rcountstwo
        real8,dimension(:), pointer :: sbuf,rbuf
        real8, dimension(:), pointer :: sbuf2nd,rbuf2nd
        int4 :: size,ierr,numout,comm
        int4 :: nearind(2), numproc(2)
        real8,dimension(:,:),pointer :: weight,val
        real8 :: x(2),eta_min(2),eta_max(2),eta_star(2)
        int4 :: ii(2), Nc(2), flag, ind(2), ell_1, ell_2,row(4)
        int4 :: contrindex(4,2)  ! store the global indexes of the comtribution points
        int4 :: globalind(2),startind(2),nx1,nx2,gridind(2),startgridind(2)
        real8 :: fieldvalue, rho
        int4 :: rank
        int4 :: l,i,j,h,m

        allocate(weight(-1:2,-1:2), val(-1:2,-1:2))

        rank=pic2d%layout2d%collective%rank
        size=pic2d%para2d%numproc(1)*pic2d%para2d%numproc(2)
        comm=pic2d%layout2d%collective%comm
        allocate(pointhead(0:size-1),curpoint(0:size-1))
        allocate(rcounts(0:size-1),scounts(0:size-1),sdispls(0:size-1),rdispls(0:size-1))  
        allocate(rcountsone(0:size-1),rcountstwo(0:size-1),scountstwo(0:size-1))
        do i=0,size-1
          allocate(pointhead(i)%ptr)
        curpoint(i)%ptr=>pointhead(i)%ptr
        end do

        call sort_quadraturepoint_among_process(mu,rank,pointhead,num_p,pic2d%para2d%N_points,pic2d)
        call mpi_alltoall(num_p,1,mpi_integer,rcounts,1,mpi_integer,comm,ierr)
        numout=0
        do i=0,size-1
          numout=numout+num_p(i)   
        end do    
        allocate(sbuf(0:numout*5-1))

        numout=0
        do i=0, size-1
          numout=numout+rcounts(i)
        end do
        allocate(rbuf(0:numout*5-1))

        do i=0,size-1
          scounts(i)=num_p(i)*5
          rcountsone(i)=rcounts(i)*5
          if(i==0) then
            sdispls(i)=0
            rdispls(i)=0
          else
            sdispls(i)=sdispls(i-1)+scounts(i-1)
            rdispls(i)=rdispls(i-1)+rcountsone(i-1)
          end if
        end do  

        h=0
        m=0
        do i=0, size-1
          do while(associated(curpoint(i)%ptr).and.associated(curpoint(i)%ptr%next))
        !if(rank==3) then  
        !m=m+1     
        !print*, "i=",i,curpoint(i)%ptr%xpoint, curpoint(i)%ptr%gridrank,curpoint(i)%ptr%gridind(1), &
        !            curpoint(i)%ptr%gridind(2),"m=",m
        !end if
            sbuf(5*h)=real(curpoint(i)%ptr%gridrank,8)
            sbuf(5*h+1:5*h+2)=curpoint(i)%ptr%xpoint
            sbuf(5*h+3)=real(curpoint(i)%ptr%gridind(1),8)
            sbuf(5*h+4)=real(curpoint(i)%ptr%gridind(2),8)
            curpoint(i)%ptr=>curpoint(i)%ptr%next
            h=h+1

          end do
        end do
        !print*, "rank=",rank, "sbuf=",sbuf 
        call mpi_alltoallv(sbuf,scounts,sdispls,mpi_double_precision,rbuf,rcountsone,rdispls, &
                        mpi_double_precision,comm,ierr)

        rho = sqrt(2.0d0*mu)
        row=pic2d%para2d%row
        numproc(1)=pic2d%para2d%numproc(1)
        numproc(2)=pic2d%para2d%numproc(2) 

        numout=0
        sdispls=0
        do i=0,size-1
        numout=numout+rcounts(i)
        scountstwo(i)=3*rcounts(i)
        if(i==0) then
          sdispls(0)=0
        else
          sdispls(i)=sdispls(i-1)+scountstwo(i-1)
        end if
        end do 
        allocate(sbuf2nd(0:numout*3-1), stat=ierr)

        numout=0
        rdispls=0
        do i=0,size-1
          numout=numout+num_p(i)
          rcountstwo(i)=3*num_p(i)
          if(i==0) then
            rdispls(i)=0
          else
            rdispls(i)=rdispls(i-1)+rcountstwo(i-1)
          end if
        end do
        allocate(rbuf2nd(0:3*numout-1),stat=ierr)

        h=0
        do i=0,size-1
          do j=1,rcounts(i)
            x=rbuf(5*h+1:5*h+2)
            call para_compute_spl2d_point_per_per_weight(weight,pic2d%para2d%m_x1,pic2d%para2d%m_x2,x,pic2d%para2d%row, &
                        pic2d%field2d%ep_weight,pic2d%field2d%epwg_w,pic2d%field2d%epwg_e,pic2d%field2d%epwg_n, &
                        pic2d%field2d%epwg_s, &
                        pic2d%field2d%epwg_sw,  pic2d%field2d%epwg_se,pic2d%field2d%epwg_ne,pic2d%field2d%epwg_nw )    
            call s_contribution_spl(x,val)
            fieldvalue=0._f64
            do m=-1,2
              do l=-1,2
                fieldvalue=fieldvalue+val(m,l)*weight(m,l)
              end do
            end do
            sbuf2nd(3*h)=fieldvalue
            sbuf2nd(3*h+1:3*h+2)=rbuf(5*h+3:5*h+4)  ! which is the lobal grid number
            h=h+1
          end do
        end do
        !print*, "rank=",rank, "sbuf2nd=",sbuf2nd
        call mpi_alltoallv(sbuf2nd,scountstwo,sdispls,mpi_double_precision,rbuf2nd,rcountstwo, &
                        rdispls,mpi_double_precision,comm,ierr) 
        !print*, "rank=",rank,"rbuf2d=",rbuf2nd 

        numout=0
        do i=1,size-1
          numout=numout+num_p(i)
        end do
        h=0
        do i=1,numout
          pic2d%field2d%epgyro(NINT(rbuf2nd(3*h+1)),NINT(rbuf2nd(3*h+2))) &
            =pic2d%field2d%epgyro(NINT(rbuf2nd(3*h+1)),NINT(rbuf2nd(3*h+2)))+rbuf2nd(3*h) 
          h=h+1
        end do

        !  deallocate(rcounts,scounts,sdispls,rdispls,rcountsone,scountstwo,rcountstwo, &
                        !             sbuf2nd,rbuf2nd,sbuf,rbuf)
        !
     end subroutine  para_compute_gyroaverage_mesh_field


   subroutine store_data_on_rootprocess(mu,mu_num,rank,rootdata,pic2d)
     class(root_precompute_data), pointer, intent(inout) :: rootdata
     class(pic_para_total2d_base), pointer, intent(in)   :: pic2d
  !   character(len=*), intent(in) :: geometry,boundary
     int4, intent(in) :: rank, mu_num
     real8, intent(in) :: mu
     class(gyropoint_node), dimension(:), pointer :: pointhead
     int4, dimension(:), pointer :: num_p
     int4 :: N_points
     int4 :: size
     int4 :: i,j
  
     N_points=pic2d%para2d%N_points
     size=pic2d%layout2d%collective%size
     select case(pic2d%para2d%geometry)
       case("cartesian")
         select case(pic2d%para2d%boundary)   
           case("double_per")
!             do j=1,mu_num
               allocate(num_p(0:size-1))
               allocate(pointhead(0:size-1)) 
               do i=0, size-1
                 allocate(pointhead(i)%ptr)
               end do 
        !!! compute the ACONTR matrix    
               call sort_quadraturepoint_among_process(mu,rank,pointhead,num_p,N_points,pic2d)
               call para_compute_gyroaverage_stiff_matrix(pointhead,num_p,mu,mu_num,N_points,pic2d,rootdata) 
!               rootdata%prematrix=rootdata%prematrix+matmul(rootdata%ACONTRI,rootdata%ASPL)*pic2d%para2d%mu_weights(j)
               deallocate(num_p)
               deallocate(pointhead)
!             enddo

            case("nat_per") 
!              call compute_D_spl2D_nat_per_noblock(pic2d%layout2d%global_sz1,pic2d%layout2d%global_sz2,rootdata%ASPL)        
        !!!! compute the ACONTR matrix    
              call sort_quadraturepoint_among_process(mu,rank,pointhead,num_p,N_points,pic2d)
              call para_compute_gyroaverage_stiff_matrix(pointhead,num_p,mu,mu_num,N_points,pic2d,rootdata) 

            case default
              stop
          end select

          case("polar")

          case default
            stop
        end select

   end subroutine store_data_on_rootprocess

   subroutine precompute_doublegyroaverage_matrix(rootdata,pic2d)
      class(root_precompute_data), pointer, intent(inout) :: rootdata
      class(pic_para_total2d_base), pointer, intent(in) :: pic2d
      real8, dimension(:), pointer :: density
      int4 :: rank,size,boxindex(4)
      int4 :: ierr, numdim,i,j,mu_num
      real8, dimension(:,:), pointer :: buf,buf1,buf2
      real8 :: mu
    
      int4 :: LDA, INFO, LWORK
      int4, dimension(:), pointer :: ipiv
      real8, dimension(:), pointer :: work

      numdim=pic2d%layout2d%global_sz1*pic2d%layout2d%global_sz2 
      rank=pic2d%layout2d%collective%rank
      size=pic2d%layout2d%collective%size
      mu_num=pic2d%para2d%mu_num
      boxindex(1)=pic2d%layout2d%boxes(rank)%i_min
      boxindex(2)=pic2d%layout2d%boxes(rank)%i_max
      boxindex(3)=pic2d%layout2d%boxes(rank)%j_min
      boxindex(4)=pic2d%layout2d%boxes(rank)%j_min
      allocate(density(numdim),stat=ierr)
      allocate(buf(boxindex(2)-boxindex(1)+1,boxindex(4)-boxindex(3)+1))

      buf= pic2d%field2d%den-pic2d%field2d%denequ    
      call gather_field_to_rootprocess_per_per(density,buf,rank, &
           size,boxindex,pic2d%para2d%numproc,pic2d%layout2d)

      if(rank==0) then
        allocate(buf1(numdim,numdim),buf2(numdim,numdim))
      end if
      do j=1,mu_num
         mu=pic2d%para2d%mu_nodes(j)   
         buf1=0.0 
         rootdata%ACONTRI=0.0  
         call store_data_on_rootprocess(mu,j,rank,rootdata,pic2d)   
         if(rank==0) then
           buf1=matmul(rootdata%ACONTRI,rootdata%ASPL)
           do i=1, numdim
             buf1(i,:)=buf1(i,:)*density(i)   ! ????
           end do
           buf1=matmul(rootdata%ASPL,buf1)
           buf2=buf2+matmul(rootdata%ACONTRI,buf1)*pic2d%para2d%mu_weights(j)         
         endif
      end do
      do i=1,numdim
        buf2(i,i)=buf2(i,i)+1.0_f64/pic2d%para2d%temp_i(i)+1.0_f64/pic2d%para2d%temp_e(i)
      end do
  
      LDA= numdim
      allocate(IPIV(NUMDIM))      
      allocate(work(NUMDIM))
      call dgetrf(numdim,numdim,buf2,LDA,IPIV,INFO)
      call dgetri(numdim,buf2,LDA,IPIV,WORK,LWORK,INFO) 
 
      rootdata%prematrix=buf2 
 
      deallocate(buf1,buf2,IPIV, work)
   end subroutine precompute_doublegyroaverage_matrix


     subroutine deri_forward_coef(p,d,coef)
        int4,intent(in) :: p,d  ! p is the order of the accuracy and d denotes the dth order of the derivative
        real8,dimension(p+d),intent(out) :: coef
        int4 :: i,j,k,imax,imin
        real8, dimension(:,:), allocatable :: buffer,buffer1,buffer2
        real8, dimension(:), allocatable :: work
        int4, dimension(:), allocatable :: ipiv
        real8 :: info
        imax=p+d-1
        imin=0
        allocate(buffer(imax+1,imax+1))
        allocate(buffer1(imax+1,imax+1))
        allocate(buffer2(imax+1,imax+1))
        allocate(work(imax+1),ipiv(imax+1))

        do i=1,imax+1
        do j=1,imax+1
          buffer(i,j)=real(j-1,8)**(i-1)
        end do
        end do

        buffer1=buffer
        buffer2=0.0

        call DGETRF(imax+1,imax+1,buffer,imax+1,ipiv,info)

        if(info.ne.0) then
          stop "matrix is numerically singulari"
        end if

        call DGETRI(imax+1,buffer,imax+1,ipiv,work,imax+1,info)

        if(info.ne.0) then
          stop "matrix inversion failed"
        end if

        do i=1,imax+1
          coef(i)=buffer(i,d+1)
        end do
        deallocate(buffer)
        deallocate(buffer1)
        deallocate(buffer2)

    end subroutine deri_forward_coef

    subroutine deri_central_coef(p,d,coef)
        int4,intent(in) :: p,d  ! p is the order of the accuracy and d denotes the dth order of the derivative
        real8,dimension(:),intent(inout) :: coef
        int4 :: i,j,k,imax,imin
        real8, dimension(:,:), allocatable :: buffer
        real8, dimension(:), allocatable :: work
        int4, dimension(:), allocatable :: ipiv
        int4 :: info, im
        imax=(p+d-1)/2
        imin=-(p+d-1)/2
        if(modulo(p+d,2)==0) then
        im=p+d-1
        else
        im=p+d
        end if
        allocate(buffer(im,im))
        allocate(work(im),ipiv(im))

        do i=1,im
          do j=1,im
            buffer(i,j)=real(imin+j-1,8)**(i-1)
          end do
        end do
        !  buffer1=buffer
        call DGETRF(im,im,buffer,im,ipiv,info)

        if(info.ne.0) then
          stop "matrix is numerically singulari"
        end if

        call DGETRI(im,buffer,im,ipiv,work,im,info)

        if(info.ne.0) then
          stop "matrix inversion failed"
        end if

        do i=1,im
          coef(i)=buffer(i,d+1)
        end do

    end subroutine deri_central_coef

    subroutine muarray_euler_maclaurin_choice(mu_bound,mu_num,mus,muweight,scheme)
        int4 :: pmax=8, d=6  ! d is the order of accuracy;pmax is the order of the derivative
        real8,intent(in) :: mu_bound
        int4,intent(in) :: mu_num,scheme
        real8, dimension(:), intent(inout) :: mus,muweight
        real8, dimension(:,:), allocatable :: coefmat
        real8, dimension(:), allocatable :: coef,bernumvec

        real8 :: binorm, bernum
        real8 :: integ
        real8 :: mul,mul1,mul2,mul3
        int4 :: i,j,k,l
        real8 :: B0
        real8 :: dvperp,vperp
        real8,dimension(:),allocatable ::trapeze(:),vptrapeze(:),g(:)

        allocate(coefmat(pmax/2,pmax+d))
        allocate(coef(pmax+d))
        allocate(bernumvec(pmax+4))
        allocate(trapeze(0:mu_num-1),vptrapeze(0:mu_num-1),g(0:mu_num-1))

        select case (scheme)

        case (1)   ! 1: forward finite difference
        do i=1,pmax/2
          coef=0.0
          call deri_forward_coef(d,i*2-1,coef)
          coefmat(i,:)=coef(:)
        end do
        case (2)  ! 2: central finite difference
        do i=1,pmax/2
          coef=0.0
          call deri_central_coef(d,i*2-1,coef)
          coefmat(i,:)=coef(:)
        end do
        case (3)  !3: michel's scheme
        do i=1,pmax/2
          coef=0.0
          call deri_central_coef(d,2*i,coef)
          coefmat(i,:)=coef(:)
        end do
        case default
        print*, "input correct scheme, scheme1=", scheme
        stop
        end select

!  print*, coefmat(2,:)

        !!Bernu number

        B0=1.0
        do i=1,pmax+4
          if(i==1)then
            bernumvec(1)=-0.5
          else
        ! compute the binormal factor
          mul1=1
          do j=1,i+1
            mul1=mul1*j
          end do
          bernum=0.0
          do k=0,i-1
            mul2=1
            mul3=1
            if(k==0) then
              mul2=1
            else
              do j=1,k
                mul2=mul2*j
              end do
            end if
            do j=1,i+1-k
              mul3=mul3*j
            end do
        ! compute the bernulli number
            binorm=real(mul1,8)/(real(mul2,8)*real(mul3,8))
            if(k==0) then
              bernum=bernum+binorm*B0
            else
              bernum=bernum+binorm*bernumvec(k)
            end if
          end do
          bernumvec(i)=-bernum/real(i+1,8)
        end if
        !           print*, "bernum",i,bernumvec(i)
      end do

      dvperp=sqrt(2.0*mu_bound)/real(mu_num-1,8)


      select case (scheme)

      case (1)
        !!! forward difference
        do i=0,mu_num-1
          if(i==0) then
            vptrapeze(i)=0.5
            do j=1,pmax/2
              vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,1)/real(2*j,8)
            end do
          else if(i.ge.1.and.i.le.(pmax+d-2)) then
            vptrapeze(i)=1.0
            do j=1,pmax/2
              if((2*j+d-2).ge.i) then
                vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,i+1)/real(2*j,8)
              end if
            end do
          else if(i==mu_num-1) then
            vptrapeze(mu_num-1)=1.0/2.0
          else
            vptrapeze(i)=1.0    !vptrapeze(i)+1.0
          end if
        end do

      case (2)
        !!!!central difference

        do i=0,mu_num-1
          if(i==0) then
            vptrapeze(i)=0.5
            do j=1,pmax/2
              vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*coefmat(j,(2*j+d+1)/2+1)/real(2*j,8)
!           print*, bernumvec(2*j), coefmat(j,(2*j+d+1)/2+1)
            end do
          else if(i.ge.1.and.i.le.(pmax+d-1)/2) then
            vptrapeze(i)=1.0
            do j=1,pmax/2
              if((2*j+d-2)/2.ge.i) then
                vptrapeze(i)=vptrapeze(i)+bernumvec(2*j)*(-coefmat(j,(2*j+d)/2-i)+ &
                        coefmat(j,(2*j+d)/2+i))/real(2*j,8)
              end if
            end do
          else if(i==mu_num-1) then
            vptrapeze(mu_num-1)=1.0/2.0
          else
            vptrapeze(i)=1.0
          end if
        end do
      case (3)
        !!!!Michel's scheme: central difference
        do i=0,mu_num-1
          if(i==0) then
            vptrapeze(i)=dvperp*bernumvec(2)/2.0
            do j=1,pmax/2
              vptrapeze(i)=vptrapeze(i)+dvperp*bernumvec(2*(j+1))*coefmat(j,(2*j+d-1)/2+1)/real(2*(j+1),8)
            end do
          else if(i.ge.1.and.i.le.(pmax+d-1)/2) then
            vptrapeze(i)=real(i,8)*dvperp
            do j=1,pmax/2
              if((2*j+d-1)/2.ge.i) then
                vptrapeze(i)=vptrapeze(i)+dvperp*bernumvec(2*j+2)*(coefmat(j,(2*j+d-1)/2+1-i)+ &
                        coefmat(j,(2*j+d-1)/2+1+i))/real(2*j+2,8)
              end if
            end do
          else if(i==mu_num-1) then
            vptrapeze(mu_num-1)=1.0/2.0*real(i,8)*dvperp
          else
            vptrapeze(i)=real(i,8)*dvperp
          end if
        end do

        case default
        print*, "input the correct scheme,scheme2=", scheme
        stop
        end select
!  print*,vptrapeze(:)

       select case (scheme)

         case (1)
           do i=1,mu_num
             vperp=real(i-1,8)*dvperp
             muweight(i)=vptrapeze(i-1)*vperp*dvperp
             mus(i)=vperp**2/2.0_f64
           end do

         case (2)
           do i=1,mu_num
             vperp=real(i-1,8)*dvperp
             muweight(i)=vptrapeze(i-1)*vperp*dvperp
             mus(i)=vperp**2/2.0_f64
           end do

         case (3)
           do i=1,mu_num
             vperp=real(i-1,8)*dvperp
             muweight(i)=vptrapeze(i-1)*dvperp
             mus(i)=vperp**2/2.0_f64
           end do
        
         case default
           print*, "input the correct scheme,scheme3=", scheme
           stop
         end select
       end subroutine muarray_euler_maclaurin_choice
 

end module para_gyroaverage_2d_one
