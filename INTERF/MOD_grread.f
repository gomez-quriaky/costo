!=====================================================================
! MP18 Module for explicite interface of GRREAD subroutine
!=====================================================================
	module MOD_grread
	 interface
    	SUBROUTINE GRREAD(fdata,ifil,modvar,rtvar,npts,grFX,grFY,grFZ,&
                        XX,XY,XZ,YY,YZ,ZZ,ptvar3,ncomp)

      character(len=80),intent(in)          :: fdata

      integer,intent(in)                      :: ifil
      integer,DIMENSION(:),intent(inout)      :: npts
      integer,intent(inout)                   :: ncomp
      integer,DIMENSION(:),intent(in)         :: modvar
      real(kind=8),DIMENSION(:),intent(in)    :: rtvar
      real(kind=8),DIMENSION(:),pointer       :: ptvar3

      real(kind=8),DIMENSION(:),pointer       :: grFX,grFY,grFZ
      real(kind=8),DIMENSION(:),pointer       :: XX,XY,XZ,YY,YZ,ZZ

	  end subroutine GRREAD
	 end interface
	end module MOD_grread