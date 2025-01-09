!=====================================================================
! Module for explicite interface of subroutine output
!=====================================================================
	module MOD_output
	 interface
          SUBROUTINE OUTPUT(iter,ddtot,grdtot,vdtot,bdtot,velvar,rovar,par,ibove,xb,yb,&
			    zb,vxnodes,vynodes,vznodes,vels,varv,varg,vargr,epsi,epsil,&
			    rmsg,rmst,rmsgr,errclt,errclgr,errclg,rtvar,ncomp)

	    integer,intent(in)                         :: iter
	    integer,DIMENSION(:),intent(in)            :: ibove

	    real(kind=8),intent(in)                    :: varv,varg,epsi,epsil
        real(kind=8),DIMENSION(6),intent(in)	   :: vargr	   
      	real(kind=8),DIMENSION(:),intent(in)       :: par
	    real(kind=8),DIMENSION(:),intent(in)       :: rmsg,rmst,errclt,errclg
	    real(kind=8),DIMENSION(:),intent(in)       :: vxnodes,vynodes,vznodes
      	real(kind=8),DIMENSION(:,:),intent(in)     :: ddtot,vdtot,bdtot
	    real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb,rovar,velvar
	    real(kind=8),DIMENSION(:,:,:,:),intent(in) :: vels
        real(kind=8),DIMENSION(:),intent(inout)    :: grdtot
		real(kind=8),DIMENSION(:,:),intent(inout)        :: rmsgr,errclgr

      real(kind=8),DIMENSION(:),intent(in)    :: rtvar         ! MP36
      integer,intent(inout)                   :: ncomp         ! MP36

	  END SUBROUTINE OUTPUT
	 end interface
	end module MOD_output
