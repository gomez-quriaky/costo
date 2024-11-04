!=====================================================================
! Module for explicite interface of subroutine caldt
!=====================================================================
	module MOD_caldt
	 interface
          SUBROUTINE CALDT(ibove,invnod,caldata,vels,ieq,aderi,par,&
						val_ad_s,columnAD,rowAD)
	   			integer,intent(in)                         :: invnod
           		integer,DIMENSION(:),intent(inout)         :: ibove
	   		integer,DIMENSION(:),pointer               :: ieq
	   		real(kind=8),DIMENSION(:),intent(in)       :: par
	   real(kind=8),DIMENSION(:),intent(inout)    :: caldata
	   real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
	   real(kind=8),DIMENSION(:,:,:,:),intent(in) :: vels
	   real(kind=8),DIMENSION(:),intent(inout) 		:: val_ad_s
	   integer, DIMENSION(:),intent(inout)   :: columnAD
	   integer, DIMENSION(:), intent(inout)  :: rowAD

	  END SUBROUTINE CALDT
	 end interface
	end module MOD_caldt
