!=====================================================================
! Module for explicite interface of subroutine SWITCH
!=====================================================================
	module MOD_switch
	 interface
	  SUBROUTINE SWITCH(a,b,ttb,ttw,nrp)
	   integer,intent(in)                           :: nrp
	   real(kind=8),intent(inout)                   :: ttb,ttw
	   real(kind=8),DIMENSION(:,:),intent(inout)    :: a,b
	  END SUBROUTINE SWITCH
	 end interface
	end module MOD_switch
