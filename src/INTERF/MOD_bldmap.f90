!=====================================================================
! Module for explicite interface for the subroutine BDLMAP
!=====================================================================
	module MOD_bldmap
	 interface
	  SUBROUTINE BLDMAP(vxnodes,vynodes,vznodes)
	   real (kind=8),DIMENSION(:),intent(in) :: vxnodes,vynodes
	   real (kind=8),DIMENSION(:),intent(in) :: vznodes
	  END SUBROUTINE BLDMAP
	 end interface
	end module MOD_bldmap

