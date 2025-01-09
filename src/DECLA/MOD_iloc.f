!=====================================================================
! Variables for location of points in arrays (BLDMAP,INTMAP...)
!=====================================================================
	module MOD_iloc
	 integer                         :: ixmax,iymax,izmax
	 integer,DIMENSION(:),pointer    :: ixloc,iyloc,izloc
	 real(kind=8)                    :: xl,yl,zl
	end module MOD_iloc

