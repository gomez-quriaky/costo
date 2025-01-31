!=====================================================================
! Module interface for MAKEMAP subroutine call
!=====================================================================
	module MOD_makemap
	 interface
	   SUBROUTINE MAKEMAP(vznodes,rline)
	    integer,intent(in)                   :: rline
            real(kind=8),DIMENSION(:),intent(in) :: vznodes
	  END SUBROUTINE MAKEMAP
	 end interface
	end module MOD_makemap
