!=====================================================================
! Module for explicite interface of subroutine TTIME
!=====================================================================
	module MOD_ttime
	 interface
	  SUBROUTINE TTIME(npt,path,tt,velco)
	   integer,intent(in)                     :: npt
	   real(kind=8),intent(out)               :: tt
           real(kind=8),DIMENSION(:),intent(in)   :: path
           real(kind=8),DIMENSION(:),intent(in)   :: velco
	  END SUBROUTINE TTIME
	 end interface
	end module MOD_ttime
