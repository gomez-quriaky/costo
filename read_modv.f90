!=====================================================================
!=====================================================================
!
!			SUBROUTINE READ_MODV
!
!=====================================================================
!=====================================================================
!> input of velocity body information
!! reads in files velmod
!!
!! puts the velocity values and their variance in the
!! array __PAR__ and __VARPAR, VELS__ and __VVARI__
!! \param PAR is the velocity perturbation relative to initial velocity
!! \param VARPAR is the covariance
!! \param VELS velocity from the input model
!! \param VVARI variance velocity from the input model
!=====================================================================
!
! Calls none
! Called by MAIN
!=====================================================================

      SUBROUTINE READ_MODV(v0,par,varpar,nnod,nbod,npar,vels,vvari,vxnodes,&
                           vynodes,vznodes)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_inv
      USE MOD_delim

      IMPLICIT NONE
!=====================================================================
! Declaration of in/out/inout arguments of READ_MODV
!=====================================================================
      integer,intent(inout)                         :: nnod,npar
      integer,intent(in)                            :: nbod

      real(kind=8),DIMENSION(:),intent(inout)       :: vxnodes
      real(kind=8),DIMENSION(:),intent(inout)       :: vynodes
      real(kind=8),DIMENSION(:),intent(inout)       :: vznodes
      real(kind=8),DIMENSION(:,:,:,:),intent(inout) :: vels,vvari

      real(kind=8),DIMENSION(:),intent(inout)       :: par,varpar,v0
!=====================================================================
! Declaration of dummy arguments of READ_MODV
!=====================================================================
      integer                         :: i,j,k,ier
      integer                         :: temp1,temp2,temp3

      character(len=80)               :: dummy
!=====================================================================
! Velocities and their variances are read from vel.mod file
!=====================================================================
      write(*,*)''
      write(*,*)'READING THE VELOCITY MODEL FILE'
      write(inout,*)''
      write(inout,*)'READING THE VELOCITY MODEL FILE'

      open(invelm,file=velmod,status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*) 'Error in opening ',velmod,'file'
         write(*,*) 'Stooooop in READ_MODV!'
         STOP
      end if
!=================================================================
! Reading the velocity and the covariance from the file velmod
!=================================================================
      read(invelm,*) dummy
      read(invelm,*) temp1,temp2,temp3
      read(invelm,*) dummy
      read(invelm,*) (vxnodes(i),i=1,nxnode)
      read(invelm,*) (vynodes(i),i=1,nynode)
      read(invelm,*) (vznodes(i),i=1,nznode)
!=================================================================
! Check the coherence between node coordinates
!=================================================================
      do i=1,nxnode-1
         if(vxnodes(i).ge.vxnodes(i+1)) then
            write(*,*)''
            write(*,*)'Inconsistency in x-node coordinates:'
            write(*,'(''xnod1= '',f7.2,'' > xnod2= '',f7.2)') vxnodes(i),&
                 vxnodes(i+1)
            STOP 'in READ_MODV'
         end if
      end do
      do i=1,nynode-1
         if(vynodes(i).ge.vynodes(i+1)) then
            write(*,*)''
            write(*,*)'Inconsistency in y-node coordinates:'
            write(*,'(''ynod1= '',f7.2,'' > ynod2= '',f7.2)') vynodes(i),&
                 vynodes(i+1)
            STOP 'in READ_MODV'
         end if
      end do
      do i=1,nznode-1
         if(vznodes(i).ge.vznodes(i+1)) then
            write(*,*)''
            write(*,*)'Inconsistency in z-node coordinates:'
            write(*,'(''znod1= '',f7.2,'' > znod2= '',f7.2)') vznodes(i),&
                 vznodes(i+1)
            STOP 'in READ_MODV'
         end if
      end do

      nnod=nxnode*nynode*nznode

      do k=1,nznode
         read(invelm,*) dummy
         do j=nynode,1,-1
            read(invelm,*) (vels(i,j,k,1),i=1,nxnode)
         end do
         do j=nynode,1,-1
            read(invelm,*) (vvari(i,j,k,1),i=1,nxnode)
         end do
         v0(k) = vels(1,1,k,1)
      end do
!=====================================================================
! zlayerdepth = depth of the (last-1) layer + 5km of the velocity model
! vhalf = velocity of the last point of the 3D raytrac model
! these are used for the simplex subroutine, as the last layer is just
! here to take into account of the layer above for raytracing
!=====================================================================
      vhalf = vels(nxnode,nynode,nznode,1)
      zlayerdepth = vznodes(nznode-1)+5.0d0
!=================================================================
! Writting some outputs in parameter.out
!=================================================================

      write(inout,*)'    OUTPUT GRID INFO FOR 3D RAYTRACING'
      write(inout,*)'       Number of xnodes (EW direction): ',nxnode
      write(inout,*)'       Number of ynodes (NS direction): ',nynode
      write(inout,*)'       Number of znodes (layers depth): ',nznode
      write(inout,*)'       Total number of nodes..........: ',nnod
      write(inout,*)' '
      write(inout,*)'       Boundings in EW direction, from ',vxnodes(1),&
                    ' to ',vxnodes(nxnode)
      write(inout,*)'       Boundings in NS direction, from ',vynodes(1),&
                    ' to ',vynodes(nynode)
      write(inout,*)'       Boundings in depth from 0 to ',vznodes(nznode)
      write(inout,*)' '
      write(inout,*)'       Velocity and Depth for starting the 3D &
                    &raytracing:'
      write(inout,'(8x,''vhalf= '',f6.3,'' km/s, zlayerdepth= '',f8.2,&
           &'' km'')') vhalf,zlayerdepth
      write(inout,*)' '

!=====================================================================
! Set ndimn for the dimensionality of the velocity model.
! This parameter is important for subroutine Vel3
! ndimn = 1 for nxnode=3 and nynode=3
! ndimn = 2 for nxnode=3 and nynode/=3
! ndimn = 3 for nxnode/=3 and nynode=3
! ndimn = 4 for nxnode/=3 and nynode/=3
!=====================================================================
      if(nxnode.eq.3.and.nynode.eq.3) ndimn=1
      if(nxnode.eq.3.and.nynode.ne.3) ndimn=2
      if(nxnode.ne.3.and.nynode.eq.3) ndimn=3
      if(nxnode.ne.3.and.nynode.ne.3) ndimn=4
      write(inout,*)'    Dimension of the model, NDIMN= ',ndimn

!=====================================================================
! Store the velocities and covariances in par and varpar arrays
! at the indice npar, which is 0 if density not inverted, NBOD if
! inverted
! PAR is the velocity perturbation relative to initial velocity.
! Thus, it is 0 here
!=====================================================================
      do k=1,nznode
         do j=1,nynode
            do i=1,nxnode
               npar=npar+1
               par(npar)= vels(i,j,k,1) - vels(i,j,k,1)
               varpar(npar)=1./(vvari(i,j,k,1)*vvari(i,j,k,1))
            end do
         end do
      end do
      ibegin(2)=iend(1)+1
      iend(2)=iend(1)+nnod
      write(inout,*)'    Velocity model stored in arrays from i=',&
           ibegin(2),' to i=',iend(2)
      write(*,*)'    Velocity model stored in arrays from i=',&
           ibegin(2),' to i=',iend(2)
        
!=================================================================
! Closing the velocity model file vel.mod
!=================================================================
      close(invelm,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file ',velmod,'logical &
                   &unit ',invelm,'. Stooooop in READ_MODV!'
         STOP
      end if




      END SUBROUTINE READ_MODV
