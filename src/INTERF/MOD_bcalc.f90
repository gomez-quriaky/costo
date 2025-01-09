!=====================================================================
! Module for explicite interface of subroutine BCALC
!=====================================================================
	module MOD_bcalc
	 interface
	  SUBROUTINE BCALC(nbod,ilay,par,ibove,iiter,thrhold)
	   integer,intent(in)                      :: nbod,iiter
           integer,DIMENSION(:),intent(in)         :: ilay,ibove
	   real(kind=8),intent(in)                 :: thrhold
           real(kind=8),DIMENSION(:),intent(inout) :: par
	  END SUBROUTINE BCALC
	 end interface
	end module MOD_bcalc
