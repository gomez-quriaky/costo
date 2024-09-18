!=====================================================================
! Module for explicite interface of subroutine calgra
!=====================================================================
	module MOD_calgra
	 interface
          SUBROUTINE CALGRA(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,&
                            yb,zb,noised,signoisd)
	   integer,intent(in)                         :: iiter,nbod,noised
	   real(kind=8),intent(in)                    :: signoisd
           real(kind=8),DIMENSION(:),intent(in)       :: par
           real(kind=8),DIMENSION(:),intent(inout)    :: caldata
           real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
           real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb
           real(kind=8),DIMENSION(:),pointer          :: FX,FY,FZ
	  END SUBROUTINE CALGRA
	 end interface
	end module MOD_calgra
