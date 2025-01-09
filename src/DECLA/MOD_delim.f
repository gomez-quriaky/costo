!=====================================================================
! Delimiters of data blocks in array storage
! ibegin = beginning of the parameters
! iend   = ending of the parameters
! jbegin = beginning of the data
! jend   = end of the data
!=====================================================================
	module MOD_delim
	  integer,DIMENSION(4) :: jbegin,jend
	  integer,DIMENSION(3) :: ibegin,iend
	end module MOD_delim
