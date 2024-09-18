!=====================================================================
! Module for explicite interface of subroutine VSMOOTH
!=====================================================================
	module MOD_vsmooth
	 interface
	  SUBROUTINE VSMOOTH(smooth,ivside,jvside,par,ibove,vxnodes,vynodes,&
			     vsum,vcont,bderi,h1,ismooth)
           integer,intent(out)                               :: vcont
	   integer,OPTIONAL,intent(in)                       :: ismooth
           integer,DIMENSION(:),intent(in)                   :: ivside,ibove
           integer,DIMENSION(:,:),intent(in)                 :: jvside
           real(kind=8),intent(out)                          :: vsum
           real(kind=8),DIMENSION(:),intent(in)              :: smooth,par
           real(kind=8),DIMENSION(:),intent(inout)           :: h1
	   real(kind=8),DIMENSION(:),intent(in)              :: vxnodes,vynodes
           real(kind=8),DIMENSION(:,:),intent(inout)         :: bderi
	  END SUBROUTINE VSMOOTH
	 end interface
	end module MOD_vsmooth
