!======================================================================
! Module for explicite interface of subroutine READBCO
!======================================================================
	module MOD_readbco
	 interface
	  SUBROUTINE READBCO(binit,varb,par,varpar,npar)
	   integer,intent(inout)                      :: npar
           real(kind=8),DIMENSION(:),intent(in)       :: binit,varb
           real(kind=8),DIMENSION(:),intent(inout)    :: par,varpar
	  END SUBROUTINE READBCO
	 end interface
	end module MOD_readbco
