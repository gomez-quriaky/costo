!=====================================================================
! Module for explicite interface of subroutine BLDMAT
!=====================================================================
	module MOD_bldmat

     !USE MKL_SPBLAS
	 interface
	  SUBROUTINE BLDMAT(iiter,aderi,bderi,punvar,varpar,h1,diff,npar,&
                        npar1,ndat,smooth,iside,jside,ivside,jvside,&
                        xb,yb,vxnodes,vynodes,par,ibove,ismooth,ilay,&
                        ddvr,nbod, val_ad_s, rowAD, columnAD, nnz)
                        !, Aderi_crs)
	   integer,intent(in)                         :: iiter,npar,npar1,ndat
           integer,intent(in)                         :: ismooth,nbod
           integer,DIMENSION(:),intent(in)            :: iside,ivside,ibove,ilay
           integer,DIMENSION(:,:),intent(in)          :: jside,jvside
           real(kind=8),DIMENSION(:),intent(in)       :: diff,smooth,par,varpar
           real(kind=8),DIMENSION(:),intent(in)       :: ddvr
           real(kind=8),DIMENSION(:),intent(in)       :: vxnodes,vynodes
           real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb
	       real(kind=8),DIMENSION(:,:),intent(in)     :: aderi
           real(kind=8),DIMENSION(:),intent(inout)    :: h1
           real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
           real(kind=8),DIMENSION(:),pointer          :: punvar
           real(kind=8), DIMENSION(:), intent(in)     :: val_ad_s
           integer, DIMENSION(:), intent(in)          :: rowAD
           integer, DIMENSION(:), intent(in)          :: columnAD
           integer, intent(in)                        :: nnz
          !type(sparse_matrix_t)                      :: Aderi_crs
	  END SUBROUTINE BLDMAT
	 end interface
	end module MOD_bldmat
