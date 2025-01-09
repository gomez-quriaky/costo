!=====================================================================
! Module for explicite interface of subroutines DDIFFER.
!=====================================================================
	module MOD_ddiffer
	 interface

	  SUBROUTINE DDIFFER(usedata,caldata,rmsg,iiter,par,par0,punvar,rovar,&
                           varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
			   diff,dtot,h1,bderi,errclg)
	   integer,intent(in)                               :: iiter,ismooth
	   integer,DIMENSION(:),intent(in)                  :: iside
           integer,DIMENSION(:,:),intent(in)                :: jside
           real(kind=8),DIMENSION(:),intent(inout)          :: rmsg,errclg
           real(kind=8),DIMENSION(:),intent(in)             :: par,par0,varpar
           real(kind=8),DIMENSION(:),intent(in)             :: smooth
	   real(kind=8),DIMENSION(:),intent(inout)          :: h1
	   real(kind=8),DIMENSION(:,:),intent(inout)        :: bderi
           real(kind=8),DIMENSION(:,:),intent(in)           :: xb,yb,zb
	   real(kind=8),DIMENSION(:),intent(inout)          :: diff
	   real(kind=8),DIMENSION(:,:),intent(inout)        :: dtot,rovar
           real(kind=8),DIMENSION(:),pointer                :: punvar
           real(kind=8),DIMENSION(:),pointer                :: usedata,caldata
	  END SUBROUTINE DDIFFER

	 end interface
	end module MOD_ddiffer
