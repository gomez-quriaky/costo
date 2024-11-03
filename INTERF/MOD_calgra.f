!=====================================================================
! Module for explicite interface of subroutine calgra
!=====================================================================
	module MOD_calgra
	 interface
          SUBROUTINE CALGRA(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,&
                            yb,zb,noised,signoisd,val_ad_s,columnAD,rowAD)
	   integer,intent(in)                         :: iiter,nbod,noised
	   real(kind=8),intent(in)                    :: signoisd
           real(kind=8),DIMENSION(:),intent(in)       :: par
           real(kind=8),DIMENSION(:),intent(inout)    :: caldata
           real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
           real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb
           real(kind=8),DIMENSION(:),pointer          :: FX,FY,FZ
           real(kind=8), DIMENSION(:), intent(inout)  :: val_ad_s
           integer, DIMENSION(:), intent(inout)        :: columnAD
           integer, DIMENSION(:), intent(inout)         ::rowAD
	  END SUBROUTINE CALGRA
	 end interface
	end module MOD_calgra
