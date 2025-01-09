!=====================================================================
! Module for explicite interface of subroutine DENVEL
!=====================================================================
	module MOD_denvel
	 interface
	  SUBROUTINE DENVEL(bderi,h1,ilay,ibove,par,ddvr,nbod)
	   integer(kind=8),intent(in)                         :: nbod
           integer(kind=8),DIMENSION(:),intent(in)            :: ilay,ibove
           real(kind=8),DIMENSION(:),intent(in)       :: par,ddvr
           real(kind=8),DIMENSION(:),intent(inout)    :: h1
           real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
	  END SUBROUTINE DENVEL
	 end interface
	end module MOD_denvel
