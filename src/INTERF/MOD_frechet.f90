!=====================================================================
! Module for explicite interface of subroutine FRECHET
!=====================================================================
	module MOD_frechet
	 interface
	  SUBROUTINE FRECHET(vxnodes,vynodes,vznodes,bazin,rayparm,&
	                    ieq,stc,ist,G,nrp,velco,invnod,signoisv,noisev)
	   integer,intent(in)                     :: noisev
	   integer,intent(out)                    :: nrp,invnod
           integer,DIMENSION(:),pointer           :: ieq,ist
	   real(kind=8),intent(in)                :: signoisv
           real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
           real(kind=8),DIMENSION(:),intent(in)   :: velco
           real(kind=8),DIMENSION(:),pointer      :: rayparm,bazin
           real(kind=8),DIMENSION(:,:),intent(inout):: G
           real(kind=8),DIMENSION(:,:),pointer    :: stc
	  END SUBROUTINE FRECHET
	 end interface
	end module MOD_frechet
