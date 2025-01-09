!=====================================================================
! Module for explicite interface of subroutine INTMAP
!=====================================================================
	module MOD_intmap
	 interface
	  SUBROUTINE INTMAP(x,y,z,ip,jp,kp)
	   integer,intent(out)                   :: ip,jp,kp
	   real(kind=8),intent(in)               :: x,y,z
	  END SUBROUTINE INTMAP
	 end interface
	end module MOD_intmap
