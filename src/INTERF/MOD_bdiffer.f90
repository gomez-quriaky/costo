!=====================================================================
! Module for explicite interface of subroutine BDIFFER
!=====================================================================
	module MOD_bdiffer
	 interface

	  SUBROUTINE BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,&
	  		     velco,v0,nbod,dtot,ibove)
	   integer,intent(in)                        :: iiter,nbod
           integer,DIMENSION(:),intent(in)           :: ilay,ibove
           real(kind=8),DIMENSION(:),intent(in)      :: par,par0,varpar
           real(kind=8),DIMENSION(:),intent(in)      :: velco,v0,varb
           real(kind=8),DIMENSION(:,:),intent(in)    :: xb,yb,zb
	   real(kind=8),DIMENSION(:,:),intent(inout) :: dtot
	  END SUBROUTINE BDIFFER

	 end interface
	end module MOD_bdiffer
