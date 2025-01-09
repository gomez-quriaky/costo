!=====================================================================
! Module for explicite interface of subroutine FORW
!=====================================================================
	module MOD_forw
	 interface
	  SUBROUTINE FORW(velco,nr,nrp,bazin,rayparm,ieq,stc,ist,G,&
                          vxnodes,vynodes,vznodes,fstime)
	   integer,intent(in)                     :: nr
	   integer,intent(out)                    :: nrp
           integer,DIMENSION(:),pointer           :: ieq,ist
	   real(kind=8),intent(out)               :: fstime
           real(kind=8),DIMENSION(:),intent(in)   :: velco
           real(kind=8),DIMENSION(:),intent(in) :: vxnodes,vynodes,vznodes
           real(kind=8),DIMENSION(:),pointer      :: rayparm,bazin
           real(kind=8),DIMENSION(:,:),intent(inout):: G
           real(kind=8),DIMENSION(:,:),pointer    :: stc
	  END SUBROUTINE FORW
	 end interface
	end module MOD_forw
