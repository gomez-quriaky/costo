!=====================================================================
! Module for explicit interface of subroutine TTMDER
!=====================================================================
	module MOD_ttmder
	 interface
	  SUBROUTINE TTMDER(G,nr,nrp,rp,velco,vxnodes,vynodes,vznodes)
	   integer,intent(in)                       :: nr,nrp
           real(kind=8),DIMENSION(:,:),intent(in)   :: rp
           real(kind=8),DIMENSION(:),intent(in)     :: velco
           real(kind=8),DIMENSION(:,:),intent(inout):: G
           real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
	  END SUBROUTINE TTMDER
	 end interface
	end module MOD_ttmder

