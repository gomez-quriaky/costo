!=====================================================================
! Module of explicite interface for data_irr subroutine
!=====================================================================
	module MOD_datairr
	 interface
	  SUBROUTINE DATA_IRR(fdata,ifil,modvar,rtvar,&
                     npts,ptvar1,gdata,gFX,gFY,gFZ,DVAL,AMASK)
	   character(len=80),intent(in)          :: fdata

           integer,intent(in)                    :: ifil
           integer,DIMENSION(:),intent(inout)    :: npts
           integer,DIMENSION(:),intent(in)       :: modvar

           real(kind=8),DIMENSION(:),pointer     :: ptvar1
           real(kind=8),DIMENSION(:),pointer     :: gdata,gFX,gFY,gFZ
           real(kind=8),DIMENSION(:),intent(inout) :: DVAL,AMASK
           real(kind=8),DIMENSION(:),intent(in)  :: rtvar
	  end SUBROUTINE DATA_IRR
	 end interface
	end module MOD_datairr
