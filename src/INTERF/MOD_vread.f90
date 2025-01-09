!=====================================================================
! Module for explicite interface of subroutine vread
!=====================================================================
	module MOD_vread
	 interface
	  SUBROUTINE VREAD(fdata,ifil,rayparm,bazin,weight,ttt,npts,ist,&
			  ieq,stc,ptvar2,modvar,rtvar)
	   character(len=80),intent(in)          :: fdata
           integer,intent(in)                    :: ifil
	   integer,DIMENSION(:),intent(in)       :: modvar
	   integer,DIMENSION(:),intent(inout)    :: npts
	   integer,DIMENSION(:),pointer          :: ieq,ist
	   real(kind=8),DIMENSION(:),intent(in)  :: rtvar
           real(kind=8),DIMENSION(:),pointer     :: rayparm,bazin,weight,ptvar2
           real(kind=8),DIMENSION(:,:),pointer   :: ttt,stc
	  end SUBROUTINE VREAD
	 end interface
	end module MOD_vread
