!=====================================================================
! Module for explicite interface of subroutine INVERMAT
!=====================================================================
	module MOD_invermat
	 interface
	  SUBROUTINE INVERMAT(bderi,npar,iiter,regul,lambda)
	   integer,intent(in)                          :: npar,iiter,regul
	   real(kind=8),intent(in)		       :: lambda
	   real(kind=8),DIMENSION(:,:),intent(inout)   :: bderi
	  END SUBROUTINE INVERMAT
	 end interface
	end module MOD_invermat
