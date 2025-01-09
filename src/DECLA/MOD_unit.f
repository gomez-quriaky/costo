!=====================================================================
! Logical units common to a lot of subroutines
! inpar = "parameter.inp" input parameter file
! inout = "parameter.out" output parameter control file
! indens = input density and variance model file
! invelm = input velocity and variance model file
! inddat = gravity data file
! invdat = delay times data file
! vpara = "raytrac.param" parameters for the velocity 3D ray tracing
! hitmap = "hitmap" file with the x,y location of rays for each interface
! raypath = "raypaths4gmt" file with the x,y,z location of rays
! inmat = "synthe.fre" file
! time1 = "time.pred" file of predicted travel times after forw subroutine
! blocnod = "bloc.nod" file with info about constrained density blocks
! grav1 = "grav.pred" file of predicted gravity anomaly after calgra
! resolden = "resol.den" file of resolution diagonal term of gravity
! resolvel = "resol.vel" file of resolution diagonal terms of velocity
! covarivel = "covari.vel" file of covariance diagonal terms of velocity
! covariden = "covari.den" file of covariance diagonal terms of density
! velout = "velocity.res" file of resulting velocity contrast
! denout = "density.res" file of resulting density
! correl = "correl.out" file of velocity/density correlation for each layer
! bcoeff = "bcoeff.res" file summing up the B-ceoff values for each iterations
!
! densmod = name of the density model
! velmod = name of the velocity model
!=====================================================================
	module MOD_unit
	   integer,PARAMETER  :: inpar=110
	   integer,PARAMETER  :: outmarie=50
	   integer,PARAMETER  :: inout=2
	   integer,PARAMETER  :: indens=3
	   integer,PARAMETER  :: invelm=4
	   integer,PARAMETER  :: inddat=120
	   integer,PARAMETER  :: invdat=8 
	   integer,PARAMETER  :: vpara=9
	   integer,PARAMETER  :: hitmap=12
 	   integer,PARAMETER  :: raypath=140
 	   integer,PARAMETER  :: inmat=150
	   integer,PARAMETER  :: time1=160
	   integer,PARAMETER  :: blocnod=170
	   integer,PARAMETER  :: grav1=180
	   integer,PARAMETER  :: resolden=182
	   integer,PARAMETER  :: resolvel=183
	   integer,PARAMETER  :: velout=184
	   integer,PARAMETER  :: denout=185
	   integer,PARAMETER  :: correl=186
	   integer,PARAMETER  :: bcoeff=187
	   integer,PARAMETER  :: covarivel=188
	   integer,PARAMETER  :: covariden=189
	   integer,PARAMETER  :: ftg1=200

	   character(len=20)  :: densmod,velmod
	end module MOD_unit
!=====================================================================
! initialization of some constants
! MAXPAR = Maximum number of parameters to be inverted for.
!          At this moment, should be 3
! MAXBOD = Maximum number of body to invert (comes from dens.mod,velm.mod)
! MAXSEG = Maximum of segments for one part of the ray
!
! iarsize = maximum number of points on raypaths
!           Sorry, can't get through this constant
!=====================================================================
	module MOD_size
!          integer,PARAMETER :: MAXPAR=3
!          integer,PARAMETER :: MAXBOD=4000
          integer,PARAMETER :: MAXSEG=7
	  integer,PARAMETER :: iarsize=300
	end module MOD_size
