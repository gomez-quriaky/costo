!=====================================================================
!=====================================================================
!
!			SUBROUTINE TTMDER
!
! Subroutine to calculate the partial derivatives with
! respect to the velocity model along the raypath.
! The corresponding partial derivatives are distributed on
! a (variable) grid read in from the file vel.mod
!
!
! Called by FORW
! Calls GMATRIX
!       VEL3 (INTMAP)
!
!=====================================================================
!=====================================================================
          INCLUDE 'INTERF/MOD_gmatrix.f'
          INCLUDE 'INTERF/MOD_vel3.f'

          SUBROUTINE TTMDER(G,nr,nrp,rp,velco,vxnodes,vynodes,vznodes)

          USE MOD_unit
          USE MOD_layer
          USE MOD_vdata
          USE MOD_iloc
          
          USE MOD_vel3
          USE MOD_gmatrix

          IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of TTMDER
!=====================================================================
          integer,intent(in)                       :: nr,nrp

          real(kind=8),DIMENSION(:,:),intent(in)   :: rp
          real(kind=8),DIMENSION(:),intent(in)     :: velco
          real(kind=8),DIMENSION(:,:),intent(inout):: G
          real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
!=====================================================================
! Declaration of the dummy arguments of TTMDER
!=====================================================================
          character(len=72)       :: forma1,forma2,forma3

          integer                 :: nrp1,i,i1,ip,jp,kp
          integer                 :: nseg

          real(kind=8)            :: tt1,v,dt
          real(kind=8)            :: txp,typ,tzp
          real(kind=8)            :: rx,ry,rz,dx1,dy1,dz1,sl,fnsegi
          real(kind=8)            :: ssl,dxs1,dys1,dzs1,xp,yp,zp

          real(kind=8),DIMENSION(:,:,:),ALLOCATABLE  :: dtm
!=====================================================================
! Initialization of some variables and allocation of some arrays
!=====================================================================
          ALLOCATE (dtm(nxnode,nynode,nznode))
          tt1 = 0.
          dtm(:,:,:) = 0.d0
!=====================================================================
! Main loop over segments comprising the ray path
!=====================================================================
          nrp1 = nrp - 1
          
          do i=1,nrp1
             i1=i+1
             rx = rp(1,i)
             ry = rp(2,i)
             rz = rp(3,i)

             dx1 = rp(1,i1)
             dy1 = rp(2,i1)
             dz1 = rp(3,i1)      
             dx1 = dx1 - rx
             dy1 = dy1 - ry
             dz1 = dz1 - rz
!=====================================================================
! Decide on the number of sub-segments and compute their lengths
! ssl = Sub Segment Length
!=====================================================================
             sl = dsqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)

             nseg = int (sl/stepl) + 1
             fnsegi = 1.0/float(nseg)
             
             ssl = sl*fnsegi
             
             dxs1 = dx1*fnsegi
             dys1 = dy1*fnsegi
             dzs1 = dz1*fnsegi
             xp = rx - 0.5d0*dxs1
             yp = ry - 0.5d0*dys1
             zp = rz - 0.5d0*dzs1
             xp = xp + dxs1
             yp = yp + dys1
             zp = zp + dzs1

!=====================================================================
! calculate the partial derivative at this point
!=====================================================================
             CALL VEL3(xp,yp,zp,v,velco,ip,jp,kp)
             dt = ssl/v
             tt1 = tt1 + dt
!=====================================================================
! find indices for point xp,yp,zp in the variable grid
! and check if node(ip,jp,kp) is outside the model
!=====================================================================

             forma1='(''raypath'',i6,'' too close to edge, xp, x_node, &
                    &ip: '',f8.2,2i4)'
             forma2='(''raypath'',i6,'' too close to edge, yp, y_node, &
                    &jp: '',f8.2,2i4)'
             forma3='(''raypath'',i6,'' too close to edge, zp, z_node, &
                    &kp: '', f8.2,2i4)'
             if(ip.gt.nxnode .or. ip.lt.1) then
                write(inout,forma1) nr,xp,nxnode,ip
             end if
             if(jp.gt.nynode .or. jp.lt.1) then
                write(inout,forma2) nr,yp,nynode,jp
             end if
             if(kp.gt.nznode .or. kp.lt.1) then
                write(inout,forma3) nr,zp,nznode,kp
             end if
!=====================================================================
! Check that the nearest node is in the velocity model and that it is
! not one of the edge nodes.  This can happen if any raypath approaches
! closer than halfway to the edge
!=====================================================================
             txp = (vxnodes(ip) + vxnodes(ip + 1))/2.d0
             typ = (vynodes(jp) + vynodes(jp + 1))/2.d0
             tzp = (vznodes(kp) + vznodes(kp + 1))/2.d0
             if (xp .ge. txp) ip = ip + 1
             if (yp .ge. typ) jp = jp + 1
             if (zp .ge. tzp) kp = kp + 1

!=====================================================================
! Record the partial derivative associated with this node 
! and count the number of velocity nodes 'hit'.
!=====================================================================
             dtm(ip,jp,kp) = dtm(ip,jp,kp) + ssl
          end do

          CALL GMATRIX(G,nr,dtm)


          DEALLOCATE (dtm)

          END SUBROUTINE TTMDER

!=====================================================================
!=====================================================================
!         SUBROUTINE GMATRIX
!
! Build the "gmatrix" which contains the partial derivatives
! for every raypath.  Note that the boundary nodes are not 
! inverted for which is why the loops range over the values
! that they do.
! params sent: nr, dtm
! params returned: G (the matrix)
!
! called by  TTMDER_VAR
! Calls none
!=====================================================================
!=====================================================================
          SUBROUTINE GMATRIX(G,nr,dtm)

          USE MOD_layer

          IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of GMATRIX
!=====================================================================
          integer,intent(in)                       :: nr

          real(kind=8),DIMENSION(:,:,:),intent(in) :: dtm
          real(kind=8),DIMENSION(:,:),intent(inout):: G

!=====================================================================
! Declaration of the dummy arguments of GMATRIX
!=====================================================================
          integer            :: i,j,k,m
!=====================================================================
! Loop over the velocity nodes and store the corresponding
! partial derivative in the matrix G
! in STECK, the loops are:
! do k = 1, nznode-1
!             do j = 2,nynode-1
!                do i = 2, nxnode-1
!=====================================================================
          m = 0
          do k = 1, nznode-1
             do j = 1,nynode
                do i = 1, nxnode
                   m = m + 1
                   G(nr,m) = dtm(i,j,k)
                end do
             end do
          end do
	
          END SUBROUTINE GMATRIX
