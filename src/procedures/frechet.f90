!=====================================================================
!=====================================================================
!
!			SUBROUTINE FRECHET
!
!> This subroutine calculates the frechet matrix (of partial derivatives)
!! and stored the number of rays passing through each node in synthe.fre
!! file
!!
!! The partial derivatives are the pathlength of the ray in the vincinity
!! of the node
!!
!! If direct problem, put the synthetic travel time in time1 file in
!! the same format as for the data file.
!!
!! It can add random gaussian noise if asked in PARAMETER.INP
!
! Calls: MAKEMAP, FORW(VEL3(INTMAP)
!                      TTMDER(GMATRIX)
!                      TTIME3(VEL3(INTMAP))
!                      SIMPLEX(SWITCH,SORT(SWITCH),TTIME2(VEL3(INTMAP)))
!                      RAYWEB(TTIME(VEL3(INTMAP)))
!                      )
! called by MAIN
!
!=====================================================================
!=====================================================================
    !INCLUDE 'INTERF/MOD_forw.f'
    !INCLUDE 'INTERF/MOD_makemap.f'

    SUBROUTINE FRECHET(vxnodes,vynodes,vznodes,bazin,rayparm,&
                       ieq,stc,ist,G,nrp,velco,invnod,signoisv,noisev)

      USE MOD_layer
      USE MOD_unit
      USE MOD_vdata
      USE MOD_iloc
      USE MOD_inv

      USE MOD_forw
      USE MOD_makemap

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of FRECHET
!=====================================================================
      integer,intent(in)                     :: noisev
      integer,intent(out)                    :: nrp,invnod
      integer,DIMENSION(:),pointer           :: ieq,ist

      real(kind=8),intent(in)                :: signoisv
      real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
      real(kind=8),DIMENSION(:),intent(in)   :: velco
      real(kind=8),DIMENSION(:),pointer      :: rayparm,bazin
      real(kind=8),DIMENSION(:,:),intent(inout):: G
      real(kind=8),DIMENSION(:,:),pointer    :: stc

!=====================================================================
! Declaration of the dummy arguments of FRECHET
!=====================================================================
      character(len=40)                     :: form1

      integer                               :: ier,nray,idum
      integer                               :: nr,i,j,k,m,ii,rline
      integer,DIMENSION(:),ALLOCATABLE      :: no_event
        
      real(kind=8)                          :: fstime
      real(kind=8),DIMENSION(:),ALLOCATABLE :: der_slo
      real(kind=8),DIMENSION(:),ALLOCATABLE :: tt_pred,tt_nois

!     character*3                         :: scount ! MP new     
!     character*20                        :: file_name  ! MP new

!=====================================================================
! Initialization of idum for random gaussian noise
!=====================================================================
      idum=-1
!=====================================================================
! Opening file of parameters of the Frechet matrix for inversion of
! teleseismic delay times INMAT
!=====================================================================
      open(inmat,file='synth.fre',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the SYNTHE.FRE file, logical &
              &unit ',inmat,'. Stooooop in FRECHET!'
         STOP
      end if

      open(raypath,file='raypath4gmt',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the RAYPATH4GMT file, logical &
              &unit ',raypath,'. Stooooop in FRECHET!'
         STOP
      end if

      open(time1,file='time.pred',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the TIME.PRED file, logical &
              &unit ',time1,'. Stooooop in FRECHET!'
         STOP
      end if

!      WRITE(scount,'(i3)') iiter
!      file_name  = "time.pred_iter_"//scount
!      open(unit = 444,file = file_name,status ='replace') 

      write(*,*)' '
      write(*,*)'3D RAY TRACING'
      write(*,*)' '

      write(inout,*)''
      write(inout,*)'3D RAY TRACING'
      write(inout,*)''
!=====================================================================
! Forward calculation of the 3D raytracing for all rays 1 -> nr
!=====================================================================
      ALLOCATE (tt_pred(n_data))
      ALLOCATE (tt_nois(n_data))
        
      write(*,*)'    number of rays to trace:',n_data
      write(*,*)''
      write(inout,*)'    number of rays to trace:',n_data
      write(inout,*)''

      rline=0
      do nr=1,n_data

         write(inout,*)'    tracing the ',nr,'-th ray ; Call FORW',rayparm(nr),bazin(nr)            
         CALL FORW(velco,nr,nrp,bazin,rayparm,ieq,stc,ist,G,&
              vxnodes,vynodes,vznodes,fstime)
         tt_pred(nr)=fstime
         tt_nois(nr)=fstime
!=====================================================================
! Counting how many lines there is in RAYPATH4GMT file for MAKEMAP
!=====================================================================
         rline=rline+nrp
      end do

      write(*,*)' '
      write(*,*)'END 3D RAYTRACING'
      write(*,*)' '

      write(inout,*)' '
      write(inout,*)'END 3D RAYTRACING'
      write(inout,*)' '
!=====================================================================
! If direct problem, put the synthetic travel time in time1 file in
! the same format as for the data file, and stop the program.
! It can add random gaussian noise if asked in PARAMETER.INP
!=====================================================================
      if(.not.INV) then
         form1='(i3,1x,i3,1x,f10.6,f8.3,f8.3,f8.3,f8.3)'
         write(time1,*)'Eq  sta    rayp      baz     tobs  tpred    &
              &tdiff  q'
         if(noisev.eq.1) then
            do nr=1,n_data
               CALL ADNOISE(tt_nois(nr),signoisv,idum)
            end do
         end if

         do nr=1,n_data
            write(time1,form1) ieq(nr),ist(nr),rayparm(nr),bazin(nr),&
                 tt_pred(nr),tt_nois(nr),(tt_nois(nr)-tt_pred(nr))
         end do

         write(*,*)''
         write(*,*)'End of forward calculation of travel times data'
         write(*,*)'Synthetic data stored in file TIME.PRED'
         write(*,*)'I''ve finished my job... BYE!'
         write(*,*)''
         write(inout,*)''
         write(inout,*)'End of forward calculation of travel times data'
         write(inout,*)'Synthetic data stored in file TIME.PRED'
         write(inout,*)'I''ve finished my job... BYE!'
         write(inout,*)''
         STOP
!=====================================================================
! If invers problem, put the synthetic travel time in time1 file just
! to check
!=====================================================================
      else
         do nr=1,n_data
            write(time1,'(f8.3)') tt_pred(nr)
!            write(444,'(f8.3)') tt_pred(nr)
         end do
!         close(444)
      end if
        
      close(raypath,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the RAYPATH4GMT file, logical &
              &unit ',raypath,'. Stooooop in FRECHET!'
         STOP
      end if

      close(time1,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the TIME.PRED file, logical &
              &unit ',time1,'. Stooooop in FRECHET!'
         STOP
      end if

!=====================================================================
! Making a map with the hits for each layer
!=====================================================================
      CALL MAKEMAP(vznodes,rline)
!=====================================================================
! Allocation of the dummy arrays no_event, der_slo
!=====================================================================
      ALLOCATE (no_event(n_data))
      ALLOCATE (der_slo(n_data))
!=====================================================================
! Starting loops over nodes
! m = number of the node
! nray = number of rays passing through the node
! invnod  = number of nodes with more than min_hit rays passing through it
! 
! in STECK, the loops are:
! do k = 1, nznode-1
!    do j = 2,nynode
!       do i = 2,nxnode
! but to fit with the parameter counter, it should be the way below...
!=====================================================================
      invnod=0
      nray = 0
      m=0
      do k = 1, nznode-1
         do j = 1,nynode
            do i = 1,nxnode
               m=m+1
               nray = 0
!=====================================================================
! Starting loops over events
!
! note which of rays 1 to n_ray hit the inversion-node (no_event)
! store the corresponding partial derivative in der_slo
!=====================================================================
               do nr=1,n_data
                  if(G(nr,m).ne.0) then
                     nray=nray+1
                     no_event(nray)=nr
                     der_slo(nray)=G(nr,m)
                  end if
               end do
!=====================================================================
! only hit nodes by at least min_hit are used
! INMAT is organized in blocks, with header line of:
!       - the number of node in the model (m)
!       - the number of rays passing through this node (nray)
!       - the new number of nodes to be inverted (invnod)
! Then for each block:
!       - the number of the event(no_event)
!       - the partial derivative(der_slo)
!=====================================================================
               if(nray.ge.min_hit) then
                  invnod=invnod+1
                  write(inmat,'(3(i5,2x))') m,nray,invnod
                  do ii=1, nray
                     write(inmat,'(i5,3x,F8.3)') no_event(ii),der_slo(ii)
                  end do
               end if
            end do
         end do
      end do

      write(inout,*)''
      write(inout,*)'Real number of velocity nodes to invert: ',invnod
!=====================================================================
! Deallocate the dummy array and closing the opened files
!=====================================================================
      DEALLOCATE (no_event)
      DEALLOCATE (der_slo)
      DEALLOCATE (tt_pred)
      DEALLOCATE (tt_nois)

      close(inmat,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the SYNTHE.FRE file, logical &
              &unit ',inmat,'. Stooooop in FRECHET!'
         STOP
      end if


    END SUBROUTINE FRECHET

!=====================================================================
!=====================================================================
!
!			SUBROUTINE MAKEMAP
!
!> This program looks for the points of each layer where
!! rays pass through. 
!Done by MJ
!
! Calls: none
! called by FRECHET
!
!=====================================================================
!=====================================================================
    SUBROUTINE MAKEMAP(vznodes,rline)

      USE MOD_unit
      USE MOD_layer
      USE MOD_vdata

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of MAKEMAP
!=====================================================================
      integer,intent(in)                   :: rline

      real(kind=8),DIMENSION(:),intent(in) :: vznodes
!=====================================================================
! Declaration of the dummy arguments of MAKEMAP
!=====================================================================
      integer                          :: ier,k,nr

      real(kind=8)                     :: x,y,z,px,py,pz
      real(kind=8)                     :: dlx,dly,dlz,ddz,l,l1
      real(kind=8)                     :: p1,p2,p3,x1,y1
!=====================================================================
! Opening some outputs for the subroutine.
! HITMAP is created here
! RAYPATHS4GMT was created in subroutine forw
!=====================================================================
      open(hitmap,file='hitmap',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file HITMAP, logical unit ',&
              hitmap,'. Stooooop in MAKEMAP!'
         STOP
      end if
        
      open(raypath,file='raypath4gmt',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file RAYPATH4GMT, logical &
              &unit ',raypath,'. Stooooop in MAKEMAP!'
         STOP
      end if

!=====================================================================
! Starting loops over the layers
!=====================================================================
      do k=2,nznode
           
         write(hitmap,*)' Layer ',k
!=====================================================================
! Starting loops over the points of a ray
!=====================================================================
         do nr=1,rline
            read(raypath,*) x,y,z
!=====================================================================
! If point belongs to the layer, writes it
!=====================================================================
            if(z.eq.vznodes(k))  then
               write(hitmap,*) x,y,z
!=====================================================================
! Else if point belongs to another layer, stores it
!=====================================================================
            else if(z.gt.vznodes(k)) then
               px=x
               py=y
               pz=z
!=====================================================================
! Else if point is between 2 layers, look for the coordinates of the 
! point on the layer:
!=====================================================================
            else if(z.lt.vznodes(k).and.pz.gt.vznodes(k)) then
               dlx=px-x
               dly=py-y
               dlz=pz-z
               ddz=pz-vznodes(k)
!=====================================================================
! calculate the distance between px,py,pz and x,y,z
!=====================================================================
               l=dsqrt(dlx**2+dly**2+dlz**2)
!=====================================================================
! determine the distance from point px,py,pz to the point in the layer
!=====================================================================
               l1=(ddz/dlz)*l
               
               x1=dlx*l1/l
               y1=dly*l1/l

               p1=px-x1
               p2=py-y1
               p3=pz-ddz

               px=x
               py=y
               pz=z

               write(hitmap,*) p1,p2,p3
            else
               px=x
               py=y
               pz=z
            end if

         end do
         REWIND(raypath)
      end do
!=====================================================================
! Closing the hitmap and raypath4gmt files
!=====================================================================
      close(hitmap,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the HITMAP file, logical unit ',&
              hitmap,'. Stooooop in MAKEMAP!'
         STOP
      end if
      close(raypath,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the RAYPATH4GMT file, logical &
              &unit ',raypath,'. Stooooop in MAKEMAP!'
         STOP
      end if

    END SUBROUTINE MAKEMAP
