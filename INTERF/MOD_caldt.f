!=====================================================================
! Module for explicite interface of subroutine caldt
!=====================================================================
	module MOD_caldt
	 interface
          SUBROUTINE CALDT(ibove,invnod,caldata,vels,ieq,aderi,par)
	   integer,intent(in)                         :: invnod
           integer,DIMENSION(:),intent(inout)         :: ibove
	   integer,DIMENSION(:),pointer               :: ieq
	   real(kind=8),DIMENSION(:),intent(in)       :: par
	   real(kind=8),DIMENSION(:),intent(inout)    :: caldata
	   real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
	   real(kind=8),DIMENSION(:,:,:,:),intent(in) :: vels
	  END SUBROUTINE CALDT
	 end interface
	end module MOD_caldt
