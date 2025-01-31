!=====================================================================
! Module for explicite interface of subroutine VDIFFER
!=====================================================================
	module MOD_vdiffer
	 interface

	  SUBROUTINE VDIFFER(usedata,caldata,rmst,iiter,par,par0,punvar,&
                          varpar,vxnodes,vynodes,ismooth,smooth,ivside,&
                          jvside,ibove,diff,dtot,velvar,h1,bderi,errclt)
	   integer,intent(in)                       :: iiter,ismooth
	   integer,DIMENSION(:),intent(in)          :: ivside,ibove
	   integer,DIMENSION(:,:),intent(in)        :: jvside
           real(kind=8),DIMENSION(:),intent(inout)  :: rmst,errclt
           real(kind=8),DIMENSION(:),intent(in)     :: par,par0,varpar
           real(kind=8),DIMENSION(:),intent(in)     :: smooth
	   real(kind=8),DIMENSION(:),intent(inout)  :: h1
           real(kind=8),DIMENSION(:),intent(in)     :: vxnodes,vynodes
	   real(kind=8),DIMENSION(:),intent(inout)  :: diff
	   real(kind=8),DIMENSION(:,:),intent(inout):: bderi
	   real(kind=8),DIMENSION(:,:),intent(inout):: dtot,velvar
           real(kind=8),DIMENSION(:),pointer        :: punvar
           real(kind=8),DIMENSION(:),pointer        :: usedata,caldata
	  END SUBROUTINE VDIFFER

	 end interface
	end module MOD_vdiffer
