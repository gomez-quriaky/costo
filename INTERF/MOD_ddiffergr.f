	module MOD_ddiffergr
	 interface

	  SUBROUTINE DDIFFERGR(usedata,caldata,rmsgr,iiter,par,par0,punvar,rovar,&
                     	varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
            			diff,dtot,grdtot,h1,bderi,npts,errclgr,rtvar,ncomp)

	    integer,intent(in)                               :: iiter
        real(kind=8),DIMENSION(:),pointer                :: usedata,caldata
        real(kind=8),DIMENSION(:,:),intent(inout)        :: rmsgr,errclgr
        real(kind=8),DIMENSION(:),intent(in)             :: par,par0        
		real(kind=8),DIMENSION(:,:),intent(inout)        :: rovar
        real(kind=8),DIMENSION(:),intent(inout)          :: diff,grdtot
		real(kind=8),DIMENSION(:),pointer                :: punvar

		! pour la partie modèle de densité (utilisé uniquement si pas d'inversion gravi avec la gradio)
		integer :: dcont
		integer,intent(in)                               :: ismooth
		integer,DIMENSION(:),intent(in)                  :: iside
		integer,DIMENSION(:,:),intent(in)                :: jside
		real(kind=8)                                     :: difs,dsum
		real(kind=8),DIMENSION(:),intent(in)             :: smooth
		real(kind=8),DIMENSION(:,:),intent(in)           :: xb,yb,zb
		real(kind=8),DIMENSION(:),intent(inout)          :: h1
		real(kind=8),DIMENSION(:,:),intent(inout)        :: bderi,dtot
		real(kind=8),DIMENSION(:),intent(in)             :: varpar

      	integer,DIMENSION(:),intent(in)      :: npts

		real(kind=8),DIMENSION(:),intent(in)    :: rtvar
		integer,intent(inout)                   :: ncomp

	  END SUBROUTINE DDIFFERGR

	 end interface
	end module MOD_ddiffergr

