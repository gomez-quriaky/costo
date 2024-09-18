!=====================================================================
! Module for explicite interface of subroutine resol
!=====================================================================
	module MOD_resol
	 interface
          SUBROUTINE RESOL(bderi,cderi,npar,xb,yb,ilay,vxnodes,vynodes,ibove)
	   integer,intent(in)                       :: npar
           integer,DIMENSION(:),intent(in)          :: ilay,ibove
           real(kind=8),DIMENSION(:),intent(in)     :: vxnodes,vynodes
           real(kind=8),DIMENSION(:,:),intent(in)   :: bderi,cderi,xb,yb
	  END SUBROUTINE RESOL
	 end interface
	end module MOD_resol
