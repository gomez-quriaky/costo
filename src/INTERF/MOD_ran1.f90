!=====================================================================
! Module for explicite interface of subroutine ran1
!=====================================================================
	module MOD_ran1
	 interface
          SUBROUTINE RAN1(idum,ran)
	   integer,intent(inout)    :: idum
	   real(kind=8),intent(out) :: ran
	  END SUBROUTINE RAN1
	 end interface
	end module MOD_ran1
