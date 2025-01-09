!=====================================================================
! Module for the common parts of velocity data and model information
!
! nstat = number of stations
! n_data = number of rays
! neq = number of earthquakes
!
! ioff = station geologic correction 1=yes, 0=no
!
! h = number of harmonics for raytracing
! i3d = 3D raytracing 1=yes, 0=no
! ioutext = extented output 1=yes, 0=no
! scale1 = step length along raypath
! stepl = step length for velocity model partial derivatives
! amp,ar,cf = amplitude, amplitude ratio, cutoff for harmonics
! nsweep = maximum number of sweep to do for the simplex algorithm
! tmin = minimum time difference between two deformed rays to stop
!        the simplex algorithm
! min_hit = minimum rays to pass through a node to invert it
!
! orlon,orlat = longitude and latitude for the center of the network
!=====================================================================
	module MOD_vdata
	  integer          :: nstat,n_data,neq
	  integer          :: ioff
	  integer          :: min_hit,nsweep
          integer          :: h,i3d,ioutext
	  real(kind=8)     :: orlon,orlat
          real(kind=8)     :: scale1,stepl,cf,ar,amp,tmin
	end module MOD_vdata
