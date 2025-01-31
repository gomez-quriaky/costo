!=====================================================================
!=====================================================================
!   SUBROUTINE STATCOORD
!=====================================================================
!=====================================================================
!
! This SUBROUTINE converts latitude and longitude coordinates of
! the stations into kilometric ones.
!
! From Bruce R. Julian, USGS Menlo Park, CA
!
! Called by RDATA
! Calls: DELAZ(CVRTOP)
!=====================================================================
    !INCLUDE 'INTERF/MOD_delaz.f'

    SUBROUTINE STATCOORD(stc,tha)

      USE MOD_vdata
      
      USE MOD_delaz

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of STATCOORD
!=====================================================================
      real(kind=8),intent(inout)            :: tha
      real(kind=8),DIMENSION(:,:),pointer   :: stc
!=====================================================================
! Declaration of the dummy arguments of STATCOORD
!=====================================================================
      integer                  :: n,i

      real(kind=8)             :: mlon,mlat,phi0,solat
      real(kind=8)             :: st0,ct0,sint,cost
      real(kind=8)             :: x,y,tn
      real(kind=8)             :: pio180
      real(kind=8),PARAMETER   :: c1 = 0.993305242609d0

      pio180 = dble(4.*atan(1.0)/180.)
!=====================================================================
! Convert coordinates to radians and correct latitude for ellipticity
!=====================================================================
      orlat = orlat*pio180
      orlon = orlon*pio180

      orlat = datan2(c1*dsin(orlat),dcos(orlat))
!=====================================================================
! Set up reference point
!=====================================================================
      st0 = dcos(orlat)
      ct0 = dsin(orlat)
      phi0 = orlon
      solat = orlat
!=====================================================================
! THA is the rotating angle counterclockwise from north to put the
! stations into the same axis of the coordinates
! I think we just do not need it
!=====================================================================
      tha = pio180*tha
      sint = dsin(tha)
      cost = dcos(tha)
!=====================================================================
! Input station list and compute relative distances (buffer through
! Character string to save input for iterative-model output)
!=====================================================================
      n=size(stc,1)
      if(n.ne.nstat) then
         write(*,*)''
         write(*,*)'Inconsistency in STATCOORD subroutine with the &
              &number of stations'
         STOP
      end if

      do i=1,nstat
!=====================================================================
! Convert degrees into radians
!=====================================================================
         stc(i,1) = stc(i,1)*pio180
         stc(i,2) = stc(i,2)*pio180
!=====================================================================
! Correct latitude for ellipticity
!=====================================================================
         stc(i,2) = datan2(c1*dsin(stc(i,2)),dcos(stc(i,2)))
!=====================================================================
! Compute x,y relative distances to orlat,orlon
!=====================================================================
         CALL delaz (stc(i,2),stc(i,1),x,y,ct0,st0,phi0,solat)

         stc(i,1) = x
         stc(i,2) = y

         if (tha .ne. 0.) then
            tn = cost*stc(i,2) - sint*stc(i,1)
            stc(i,1) = cost*stc(i,1) + sint*stc(i,2)
            stc(i,2) = tn
         end if

      end do

    END SUBROUTINE STATCOORD

!=====================================================================
!=====================================================================
!   SUBROUTINE DELAZ
!=====================================================================
!=====================================================================
!
! This SUBROUTINE calculates the geocentric distance and azimuths
!
! Called by STATCOORD
! Calls: CVRTOP
!=====================================================================

    SUBROUTINE DELAZ(lat,lon,x,y,ct0,st0,phi0,solat)

      USE MOD_cvrtop
      
      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of DELAZ
!=====================================================================
      real(kind=8),intent(in)      :: lat,lon,ct0,st0,phi0,solat
      real(kind=8),intent(out)     :: x,y
!=====================================================================
! Declaration of the dummy arguments of DELAZ
!=====================================================================
      real(kind=8),PARAMETER       :: erad = 6378.137
      real(kind=8)                 :: twopi,flat,pi
      real(kind=8)                 :: ct1,st1,delta,az0,az1
      real(kind=8)                 :: sdlon,cdlon,sdelt,cdelt
      real(kind=8)                 :: lambda,temp,radius,colat
!=====================================================================
! FLAT: Earth flattening constant (Chovitz, 1981, EOS, 62, 65-67)
! ERAD: Earth equatorial radius (Chovitz, 1981, EOS, 62, 65-67)
!=====================================================================
      pi = dble(4.*atan(1.0))
      twopi= 2.*pi
      flat = 1./298.257

      ct1 = dsin(lat)
      st1 = dcos(lat)

      if ((ct1 - ct0).eq.(0.) .and. (lon - phi0).eq.(0.)) then
        delta = 0.
        az0 = 0.
        az1 = 0.
      else
        sdlon = dsin(lon - phi0)
        cdlon = dcos(lon - phi0)
        cdelt = st0*st1*cdlon + ct0*ct1
!=====================================================================
! Calculate the polar coordinates sdelt,az0 from the cartesian one
!=====================================================================
        CALL cvrtop (st0*ct1 - st1*ct0*cdlon,st1*sdlon,sdelt,az0)

        delta = datan2(sdelt,cdelt)

        CALL cvrtop (st1*ct0 - st0*ct1*cdlon,-sdlon*st0,sdelt,az1)

        if (az0 .lt. 0.) az0 = az0 + twopi
        if (az1 .lt. 0.) az1 = az1 + twopi
      end if

      colat = pi/2. - (lat + solat)/2.

      lambda=1.d0 - flat
      lambda=lambda*lambda
      lambda=flat*(2.d0 - flat)/lambda
      temp=dcos(colat)
      temp=temp*temp
      temp=1.d0 + lambda*temp
      radius=erad/dsqrt(temp)

      y = radius*delta*dcos(az0)
      x = radius*delta*dsin(az0)


    END SUBROUTINE DELAZ

!=====================================================================
!=====================================================================
!   SUBROUTINE CVRTOP
!=====================================================================
!=====================================================================
!
! This SUBROUTINE converts from rectangular to polar coordinates
!
! From Bruce Julian, USGS Menlo Park, CA
!
! Called by DELAZ
! Calls: none
!=====================================================================

    SUBROUTINE CVRTOP(x,y,r,theta)
      
      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of CVRTOP
!=====================================================================
      real(kind=8),intent(in)      :: x,y
      real(kind=8),intent(out)     :: r,theta


      r = DSQRT(x*x + y*y)
      theta = datan2(y,x)

    END SUBROUTINE CVRTOP
