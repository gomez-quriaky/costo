!=====================================================================
! Module interface for RDATA subroutine call
!=====================================================================
	module MOD_rdata
	 interface
	   SUBROUTINE RDATA(nfil,codata,dfile,modvar,rtvar,&
                        FX,FY,FZ,DVAL,AMASK,usedata,npts,rayparm,bazin,&
			weight,ndat,ist,stc,ieq,ttt,punvar,varv,varg,vargr,ncomp)

	    character(len=80),DIMENSION(:),intent(in) :: dfile

	    integer,intent(in)                  :: nfil
        integer,DIMENSION(:),intent(in)     :: codata,modvar
	    integer,DIMENSION(:),pointer        :: ieq,ist
        integer,DIMENSION(nfil),intent(out) :: npts
	    integer,intent(out)                 :: ndat

	    real(kind=8),intent(out)                :: varv,varg
        real(kind=8),DIMENSION(6),intent(out)   :: vargr

        real(kind=8),DIMENSION(:),intent(in)    :: rtvar
        real(kind=8),DIMENSION(:),pointer       :: FX,FY,FZ,punvar
        real(kind=8),DIMENSION(:),intent(inout) :: DVAL,AMASK
	    real(kind=8),DIMENSION(:),pointer       :: rayparm,bazin,weight
        real(kind=8),DIMENSION(:),pointer       :: usedata
	    real(kind=8),DIMENSION(:,:),pointer     :: stc,ttt
	    
	    integer,intent(inout)                   :: ncomp

	   END SUBROUTINE RDATA
	 end interface
	end module MOD_rdata
