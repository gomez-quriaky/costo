!=====================================================================
! Module for explicite interface of subroutine DMFSD
!=====================================================================
	module MOD_dmfsd
	 interface
	  SUBROUTINE DMFSD(mat,npar,toler,ier)
	   integer,intent(in)                         :: npar
	   integer,intent(inout)                      :: ier
	   real(kind=8),intent(in)                    :: toler
	   real(kind=8),DIMENSION(:),intent(inout)    :: mat
	  END SUBROUTINE DMFSD
	 end interface
	end module MOD_dmfsd
