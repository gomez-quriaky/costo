!=====================================================================
! Module for explicite interface of subroutine timecal
!=====================================================================
	module MOD_timecal
	 interface
          SUBROUTINE TIMECAL(date1,date2)
	   integer,DIMENSION(8),intent(inout)  :: date1,date2
	  END SUBROUTINE TIMECAL
	 end interface
	end module MOD_timecal
