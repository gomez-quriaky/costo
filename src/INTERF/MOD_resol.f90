!=====================================================================
! Module for explicite interface of subroutine resol
!=====================================================================
	module MOD_resol
	 interface
          SUBROUTINE RESOL(bderi,aderi,npar,ndat,xb,yb,ilay,vxnodes,vynodes,ibove,punvar)
	   integer,intent(in)                       :: npar,ndat
           integer,DIMENSION(:),intent(in)          :: ilay,ibove
           real(kind=8),DIMENSION(:),intent(in)     :: vxnodes,vynodes
           real(kind=8),DIMENSION(:,:),intent(in)   :: bderi,aderi,xb,yb
	   real(kind=8),DIMENSION(:),pointer        :: punvar
	  END SUBROUTINE RESOL
	 end interface
	end module MOD_resol
