!=====================================================================
! Module for explicite interface of subroutine DSMOOTH
!=====================================================================
	module MOD_dsmooth
	 interface
	  SUBROUTINE DSMOOTH(smooth,iside,jside,xb,yb,par,dsum,&
			     dcont,bderi,h1,ismooth)
           integer,intent(out)                               :: dcont
	   integer,OPTIONAL,intent(in)                       :: ismooth
           integer,DIMENSION(:),intent(in)                   :: iside
           integer,DIMENSION(:,:),intent(in)                 :: jside
           real(kind=8),intent(out)                          :: dsum
           real(kind=8),DIMENSION(:,:),intent(in)            :: xb,yb
           real(kind=8),DIMENSION(:),intent(in)              :: smooth,par
           real(kind=8),DIMENSION(:),intent(inout)           :: h1
           real(kind=8),DIMENSION(:,:),intent(inout)         :: bderi
	  END SUBROUTINE DSMOOTH
	 end interface
	end module MOD_dsmooth
