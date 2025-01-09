!=====================================================================
! Module for explicite interface of subroutine cvrtop
!=====================================================================
	module MOD_cvrtop
	 interface
          SUBROUTINE CVRTOP(x,y,r,theta)
	   real(kind=8),intent(in)      :: x,y
           real(kind=8),intent(out)     :: r,theta
	  END SUBROUTINE CVRTOP
	 end interface
	end module MOD_cvrtop
