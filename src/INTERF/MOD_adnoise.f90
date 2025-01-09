!=====================================================================
! Module for explicite interface of subroutine adnoise
!=====================================================================
	module MOD_adnoise
	 interface
          SUBROUTINE ADNOISE(syndata,signois,idum)
	   integer,intent(inout)       :: idum
           real(kind=8),intent(in)     :: signois
           real(kind=8),intent(inout)  :: syndata
	  END SUBROUTINE ADNOISE
	 end interface
	end module MOD_adnoise
