!=====================================================================
!=====================================================================
!
!			SUBROUTINE VEL3
!
!
! MODIFIED BY W. PROTHERO TO USE VELOCITY INTERPOLATION COEFFICIENTS
! COMPUTED IN SUBROUTINE VCOEFS. THIS ALLOWS THE VELOCITIES TO BE
! COMPUTED MUCH FASTER, SINCE THE INTERPOLATION COEFFICIENTS
! DON'T HAVE TO BE RECOMPUTED FOR EACH POINT.
! params sent:  x,y,z (location), velco (vel coeffs)
! params returned:  v (velocity at x, y, z)
!
! Called from FORW, TTIME, SIMPLEX
! Calls INTMAP
!
!=====================================================================
!=====================================================================
      INCLUDE 'INTERF/MOD_intmap.f'

      SUBROUTINE VEL3(x,y,z,v,velco,ip,jp,kp)

      USE MOD_size
      USE MOD_iloc
      USE MOD_layer

      USE MOD_intmap

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of VEL3
!=====================================================================
      integer,intent(out)                   :: ip,jp,kp

      real(kind=8),intent(in)               :: x,y,z
      real(kind=8),intent(out)              :: v
      real(kind=8),DIMENSION(:),intent(in)  :: velco
!=====================================================================
! Declaration of the dummy arguments of VEL3
!=====================================================================
      integer                  :: indey

      real(kind=8)             :: xy,xyz,xz,yz
      
!=====================================================================
! Get the location of the point x,y,z
!=====================================================================
      CALL INTMAP(x,y,z,ip,jp,kp)
!=====================================================================
! Compute index of ip,jp,kp coefficient in velco array
!=====================================================================
      indey=8*(ip-1+nxnode*(jp-1)+nxnode*nynode*(kp-1))
      if (indey .lt. 0) then
         write(*,*) 'Indey < 0 in Vel3 subroutine, indey,ip,jp,kp,x,y,z'
         write(*,*) indey,ip,jp,kp,x,y,z
      end if
      
      select case (ndimn)
         case (4)
            xy = x*y
            xyz = xy*z
            xz = x*z
            yz = y*z
            v = xyz*velco(1+indey)+xy*velco(2+indey)+yz*velco(3+indey)&
                +xz*velco(4+indey)+x*velco(5+indey)+y*velco(6+indey)&
                +z*velco(7+indey)+velco(8+indey)
! the velocity dependence is only on z
         case (1) 
            v=z*velco(7+indey)+velco(8+indey)
! the velocity dependence is only on z and y
         case (2)
            v=z*velco(7+indey)+velco(8+indey)+y*velco(6+indey)&
              +y*z*velco(3+indey)
! the velocity dependence is only on z and x
         case (3)
            v=z*velco(7+indey)+velco(8+indey)+x*velco(5+indey)&
              +x*z*velco(4+indey)
            
      end select

      END SUBROUTINE VEL3

!=====================================================================
!=====================================================================
!
!			SUBROUTINE INTMAP
!
! MODIFIED BY W. PROTHERO SO A SINGLE CALL CAN GET THE INDICES
!
! Called by VEL3
! Calls none
!=====================================================================
!=====================================================================

      SUBROUTINE INTMAP(x,y,z,ip,jp,kp)

      USE MOD_size
      USE MOD_layer
      USE MOD_iloc

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of TTIME
!=====================================================================
      integer,intent(out)                   :: ip,jp,kp

      real(kind=8),intent(in)               :: x,y,z
!=====================================================================
! Declaration of the dummy arguments of TTIME
!=====================================================================
      ip=INT(x+xl)
      ip=ixloc(ip)

      jp=INT(yl+y)
      jp=iyloc(jp)

      kp=INT(zl+z)
      if(kp.gt.izmax) then
         write(*,*)' Funky in intmap, kp > izmax'
         write(*,*)'x,y,z,kp,izloc =',x,y,z,kp,izloc(kp)
         write(*,*) 'Stooooop in INTMAP!'
         STOP
      else if (kp .lt. 1) then
         write(*,*) ' Funky kp < 1 , x,y,z,kp,izloc =',x,y,z,kp,izloc(1)
         kp = izloc(1)
!         STOP
      else
         kp=izloc(kp)
      end if

      END SUBROUTINE INTMAP

