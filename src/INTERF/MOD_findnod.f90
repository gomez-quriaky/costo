!=====================================================================
! Module for explicite interface of subroutine findnod
!=====================================================================
	module MOD_findnod
	 interface
          SUBROUTINE FINDNOD(xb,yb,zb,vxnodes,vynodes,vznodes,nbod,ibove)
	   integer,intent(in)                     :: nbod
           integer,DIMENSION(:),intent(in)        :: ibove
           real(kind=8),DIMENSION(:,:),intent(in) :: xb,yb,zb
	   real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
	  END SUBROUTINE FINDNOD
	 end interface
	end module MOD_findnod
