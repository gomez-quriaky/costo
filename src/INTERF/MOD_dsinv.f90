!=====================================================================
! Module for explicite interface of subroutine DSINV
!=====================================================================
	module MOD_dsinv
	 interface
	  SUBROUTINE DSINV(bderi,npar,toler,ier)
	   integer,intent(in)                         :: npar
	   integer,intent(inout)                      :: ier
	   real(kind=8),intent(in)                    :: toler
	   real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
	  END SUBROUTINE DSINV
	 end interface
	end module MOD_dsinv
