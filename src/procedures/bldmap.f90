!=====================================================================
!=====================================================================
!   SUBROUTINE BLDMAP
!=====================================================================
!=====================================================================
!> set up an array which is used to point to node indices, for any x,y,z
!! for the 3D raytracing velocity model
! 
! Called from main, useful for INTMAP subroutine
! calls none
!=====================================================================

    SUBROUTINE BLDMAP(vxnodes,vynodes,vznodes)

      USE MOD_unit
      USE MOD_vdata
      USE MOD_layer
      USE MOD_iloc
        
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of BLDMAP
!=====================================================================
      real(kind=8),DIMENSION(:),intent(in)    :: vxnodes,vynodes,vznodes
!=====================================================================
! Declaration of the dummy arguments of BLDMAP
!=====================================================================
      integer                   :: i,ix,ix1,iy,iy1,iz,iz1
        
      real(kind=8)              :: xnow,ynow,znow
!=====================================================================
! Determination of the array dimension
!
! ixmax = distance in km between first and last nodes in x direction
! iymax = id. in y direction
! izmax = id. in z direction
!=====================================================================
      xl=1.d0 - vxnodes(1)
      ixmax=int(vxnodes(nxnode)+xl)

      yl=1.d0 - vynodes(1)
      iymax=int(vynodes(nynode)+yl)

      zl=1.d0 - vznodes(1)
      izmax=int(vznodes(nznode)+zl)

      write(inout,*)' '
      write(inout,*)'BUILDING MAP FROM THE VELOCITY MODEL'
      write(inout,*)' '
      write(inout,*)'    ARRAY SIZES IN KM (from vel.mod): '
      write(inout,*)'       -in X (EW) direction: ',ixmax-1
      write(inout,*)'       -in Y (NS) direction: ',iymax-1
      write(inout,*)'       -in Z (depth) direction: ',izmax-1
      write(inout,*)' '
!=====================================================================
! Allocation of the array ixloc
!=====================================================================
      ALLOCATE (ixloc(ixmax))
      ALLOCATE (iyloc(iymax))
      ALLOCATE (izloc(izmax))
!=====================================================================
! Storage of the array ixloc
!=====================================================================
      ix=1
      do i=1,ixmax
         ix1=ix+1
         xnow=dble(i)-xl
         if (xnow.ge.vxnodes(ix1)) ix=ix1
         ixloc(i)=ix
      end do
!=====================================================================
! Storage of the array iyloc
!=====================================================================
      iy=1
      do i=1,iymax
         iy1=iy+1
         ynow=dble(i)-yl
         if (ynow.ge.vynodes(iy1)) iy=iy1
         iyloc(i)=iy
      end do
!=====================================================================
! depth determinations according to velocity model, not density (as
! previously done by MJ)
!=====================================================================

      iz=1
      do i=1,izmax
         iz1=iz+1
         znow=dble(i)-zl
         if (znow.ge.vznodes(iz1)) iz=iz1
         izloc(i)=iz
      end do

    END SUBROUTINE BLDMAP


!=====================================================================
!=====================================================================
!   SUBROUTINE VCOEFS
!=====================================================================
!=====================================================================
!> This routine computes the coefficients for the linear velocity
!! interpolation used in VEL3.
!!
!! It has been added by W. PROTHERO to increase the efficiency
!! of the travel time computations
! 
! Called by MAIN
! calls none
!=====================================================================
    SUBROUTINE VCOEFS(velco,iiter,vels,vxnodes,vynodes,vznodes)

       
      USE MOD_unit
      USE MOD_vdata
      USE MOD_iloc
      USE MOD_layer
       
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of VCOEFS
!=====================================================================
      integer,intent(in)                          :: iiter

      real(kind=8),DIMENSION(:),intent(inout)     :: velco
      real(kind=8),DIMENSION(:,:,:,:),intent(in)  :: vels
      real(kind=8),DIMENSION(:),intent(in)        :: vxnodes,vynodes,vznodes

!=====================================================================
! Declaration of the dummy arguments of VCOEFS
!=====================================================================

      integer                         :: i,j,k,INDEY

      real(kind=8)                    :: V111,V112,V121,V211,V212,&
                                          V221,V222,V122
      real(kind=8)                    :: X1,X2,Y1,Y2,Z1,Z2,dxdydz,CY

!=====================================================================
! Rule for getting the index of velco:
! INDEX=ICOEF+8*(IP-1+nxnode*(JP-1)+nxnode*nynode*(KP-1)
!=====================================================================

      do i=1,nxnode-1
         do j=1,nynode-1
            do k=1,nznode-1
               INDEY=8*(i-1+nxnode*(j-1)+nxnode*nynode*(k-1))
               dxdydz=(vxnodes(i+1)-vxnodes(i))*(vynodes(j+1)-&
                    vynodes(j))*(vznodes(k+1)-vznodes(k))
               if (nxnode.le.1) then
                  dxdydz=(vynodes(j+1)-vynodes(j))*(vznodes(k+1)-vznodes(k))
               end if
               V111=vels(i,j,k,iiter)
               V121=vels(i,j+1,k,iiter)
               V211=vels(i+1,j,k,iiter)
               V112=vels(i,j,k+1,iiter)
               V122=vels(i,j+1,k+1,iiter)
               V212=vels(i+1,j,k+1,iiter)
               V221=vels(i+1,j+1,k,iiter)
               V222=vels(i+1,j+1,k+1,iiter)
                
               X1=vxnodes(i)
               X2=vxnodes(i+1)
               Y1=vynodes(j)
               Y2=vynodes(j+1)
               Z1=vznodes(k)
               Z2=vznodes(k+1)

               velco(INDEY+1)=(-V111+V211+V121-V221+V112-V212-V122+V222)&
                    /dxdydz

               velco(INDEY+2)=(Z2*(V111-V211-V121+V221)&
                    +Z1*(-V112+V212+V122-V222))/dxdydz

               CY=X2*(V111-V121-V112+V122)+X1*(-V211+V221+V212-V222)
               velco(INDEY+3)=CY/dxdydz

               CY=Y2*(V111-V211-V112+V212)+Y1*(-V121+V221+V122-V222)
               velco(INDEY+4)=CY/dxdydz

               CY=Y2*Z2*(-V111+V211)+Y1*Z2*(V121-V221)+Y2*Z1*(V112-V212)&
                    +Y1*Z1*(-V122+V222)
               velco(INDEY+5)=CY/dxdydz

               CY=X2*Z2*(-V111+V121)+X1*Z2*(V211-V221)+X2*Z1*(V112-V122)&
                    +X1*Z1*(-V212+V222)
               velco(INDEY+6)=CY/dxdydz

               CY=X2*Y2*(-V111+V112)+X1*Y2*(V211-V212)+X2*Y1*(V121-V122)&
                    +X1*Y1*(-V221+V222)
               velco(INDEY+7)=CY/dxdydz

               CY=X2*Y2*Z2*V111-X1*Y2*Z2*V211-X2*Y1*Z2*V121-X2*Y2*Z1*V112&
                    +X1*Y1*Z2*V221+X1*Y2*Z1*V212+X2*Y1*Z1*V122-X1*Y1*Z1*V222
               velco(INDEY+8)=CY/dxdydz

            end do
         end do
      end do

    END SUBROUTINE VCOEFS

