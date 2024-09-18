!=====================================================================
!
! nxbloc,nybloc : maximum number of blocks for density model
! nxnode,nynode,nznode : maximum number of nodes for velocity model
! nlayer        : exact number of layers for BOTH model
! vhalf : velocity of the last node for 3D raytracing
! zlayerdepth : depth of the nlayer+5km for 3D raytracing
! ndimn : dimension of the model according to number of nodes
! 
!=====================================================================
	module MOD_layer
	  integer       :: nlayer,ndimn
	  integer       :: nxbloc,nybloc
	  integer       :: nxnode,nynode,nznode
	  real(kind=8)  :: vhalf,zlayerdepth
	end module MOD_layer

