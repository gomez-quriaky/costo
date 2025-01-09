!=====================================================================
! Module for explicite interface of subroutine statcoord
!=====================================================================
	module MOD_statcoord
	 interface
          SUBROUTINE STATCOORD(stc,tha)
	   real(kind=8),intent(inout)            :: tha
	   real(kind=8),DIMENSION(:,:),pointer   :: stc
	  END SUBROUTINE STATCOORD
	 end interface
	end module MOD_statcoord
