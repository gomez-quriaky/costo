	module MOD_calgradio
	 interface

	  SUBROUTINE CALGRADIO(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
	  							noised,signoisd,rtvar,ncomp)

      integer,intent(in)                         :: iiter,nbod,noised
      real(kind=8),intent(in)                    :: signoisd      
      real(kind=8),DIMENSION(:),intent(in)       :: par      
      real(kind=8),DIMENSION(:),intent(inout)    :: caldata
      real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
      real(kind=8),DIMENSION(:),pointer          :: FX,FY,FZ
      real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb

		real(kind=8),DIMENSION(:),intent(in)    :: rtvar
		integer,intent(inout)                   :: ncomp
		
	  END SUBROUTINE CALGRADIO

	 end interface
	end module MOD_calgradio