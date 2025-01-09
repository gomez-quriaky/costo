!=====================================================================
! Module for explicite interface of subroutine perturb
!=====================================================================
	module MOD_perturb
	 interface
          SUBROUTINE PERTURB(bderi,h1,npar,drmax,ilay,par,par0,vels,iiter,&
			     ibove,xb,yb,zb)
	   integer,intent(in)                            :: npar,iiter
           integer,DIMENSION(:),intent(in)               :: ilay,ibove
           real(kind=8),DIMENSION(:),intent(in)          :: h1,drmax
	   real(kind=8),DIMENSION(:),intent(inout)       :: par,par0
           real(kind=8),DIMENSION(:,:),intent(in)        :: bderi,xb,yb,zb
           real(kind=8),DIMENSION(:,:,:,:),intent(inout) :: vels
	  END SUBROUTINE PERTURB
	 end interface
	end module MOD_perturb
