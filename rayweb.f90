!=====================================================================
!=====================================================================
!
!			SUBROUTINE RAYWEB
!
! This routine just calculates the travel time for a straight 
! ray from the incident wave to the station.
!
! Called by FORW
! Calls TTIME(VEL3(INTMAP))
!
!=====================================================================
!=====================================================================
      INCLUDE 'INTERF/MOD_ttime.f'

      SUBROUTINE RAYWEB(xe,ye,ze,xs,ys,zs,fstime,rp,nrp,velco,nr)

      USE MOD_unit
      USE MOD_size
      USE MOD_layer
      USE MOD_vdata
      USE MOD_iloc

      USE MOD_ttime

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of RAYWEB
!=====================================================================
      integer,intent(out)                :: nrp
      integer,intent(in)                 :: nr
      real(kind=8),intent(in)            :: xs,ys,zs,xe,ye,ze
      real(kind=8),intent(out)           :: fstime

      real(kind=8),DIMENSION(:,:),intent(inout) :: rp
      real(kind=8),DIMENSION(:),intent(in)      :: velco

!=====================================================================
! Declaration of the dummy arguments of RAYWEB
!=====================================================================
      integer                    :: i
      integer                    :: ns1,nd,ns,npt,ic,n1,n2,n3
      real(kind=8)               :: a2log
      real(kind=8)               :: delx,dely,delz,sep
      real(kind=8)               :: sn1,xstep,ystep,zstep
      real(kind=8),DIMENSION(:),ALLOCATABLE :: trpth1

!=====================================================================
! Initialization of some arrays
!=====================================================================
      fstime = 0.0
      ns1 = 0
      ic = 0

      ALLOCATE (trpth1(3*iarsize))
!=====================================================================
! Compute source receiver distance
!=====================================================================
      a2log = 1.d0 / dlog10(2.d0)
      delx = xs-xe
      dely = ys-ye
      delz = zs-ze
      sep = dsqrt(delx*delx+dely*dely+delz*delz)
!=====================================================================
! Determine the number of path divisions for the ray
!=====================================================================
      nd = 1 + int(a2log * dlog10(sep/scale1))
      if (nd.gt.MAXSEG) then
         nd = MAXSEG
      end if
!=====================================================================
! Determine the number of path segments
!=====================================================================
      ns = 2**nd
!=====================================================================
! number of points on the path
!=====================================================================
      npt = ns + 1
      
      if (sep.le.scale1) then
         nd = 0
         ns = 1
         npt = 2
      end if
!=====================================================================
! Check wether indices will exceed array sizes
!=====================================================================
      if (npt .gt. iarsize) then
         write (*,*)''
         write (*,*)'Array size exceeded: npt (',npt,') greater than &
               &iarsize (',iarsize,'). Stooooop in RAYWEB !'
         STOP
      end if
!=====================================================================
!Determine points along straight-line path
!=====================================================================
      sn1   = 1.d0/ ns
      xstep = delx*sn1
      ystep = dely*sn1
      zstep = delz*sn1
      
      ns1=ns+1
      do i=1,ns1
         ic=ic+1
         trpth1(ic) = xe+xstep*(i-1)
         ic=ic+1
         trpth1(ic) = ye+ystep*(i-1)
         ic=ic+1
         trpth1(ic) = ze+zstep*(i-1)
      end do
!=====================================================================
! calculate the straight line travel time
!=====================================================================
      CALL TTIME(npt,trpth1,fstime,velco)
!=====================================================================
! transfer the raypath from a single to a three dimensional array
!=====================================================================
      do i=1,npt
         n3=3*i  
         n1=n3-2
         n2=n1+1
         rp(1,i)=trpth1(n1)
         rp(2,i)=trpth1(n2)
         rp(3,i)=trpth1(n3)
      end do
  
      nrp = npt
      DEALLOCATE (trpth1)
      END SUBROUTINE RAYWEB


!=====================================================================
!=====================================================================
!
!			SUBROUTINE TTIME
!
! travel time along path via trapezoidal rule integration
!
! Called by RAYWEB
! Calls VEL3(INTMAP)
!
!=====================================================================
!=====================================================================

      SUBROUTINE TTIME(npt,path,tt,velco)

      USE MOD_size
      USE MOD_iloc
      USE MOD_layer

      USE MOD_vel3

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of TTIME
!=====================================================================
      integer,intent(in)                     :: npt

      real(kind=8),intent(out)               :: tt
      real(kind=8),DIMENSION(:),intent(in)   :: path
      real(kind=8),DIMENSION(:),intent(in)   :: velco
!=====================================================================
! Declaration of the dummy arguments of TTIME
!=====================================================================
      integer                         :: j,i,k,l,i1,ip,jp,kp

      real(kind=8)                    :: dt,x,y,z,dx2,dy2,dz2,psep,v
      REAL(kind=8),DIMENSION(iarsize) :: v1

!=====================================================================
! Initialization of some variables
!=====================================================================
      tt = 0.0d0
      dt = 0.0d0
!=====================================================================
! compute the travel time
!=====================================================================
      j = 1
      if(size(path).lt.npt) then
         write(*,*)'Error in RAYWEB subroutine, dimension of array &
              & exceeds allowed one (npt= ',npt,'). Stooooop in TTIME!'
         STOP
      end if
      do i=1,npt
         k = j + 1
         l = j + 2
         x = path(j)
         y = path(k)
         z = path(l)
         j = j + 3
         CALL VEL3(x,y,z,v,velco,ip,jp,kp)
         v1(i) = v
      end do

      j = 1
      do i=2,npt
         k = j + 1
         l = j + 2
         i1 = i - 1
         dx2 = path(j) - path(j+3)
         dy2 = path(k) - path(k+3)
         dz2 = path(l) - path(l+3)
         j = j + 3
         psep = dsqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)
         if (v1(i) .ne. 0.d0 .and. v1(i1).ne. 0.d0 ) then
            dt = psep*(1.d0/v1(i) + 1.d0/v1(i1))
         end if
         tt = tt + dt
      end do
      tt = 0.5d0*tt

      END SUBROUTINE TTIME

