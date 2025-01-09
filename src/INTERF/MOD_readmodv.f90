!======================================================================
! Module for explicite interface of subroutine READ_MODV
!======================================================================
	module MOD_readmodv
	 interface
	  SUBROUTINE READ_MODV(v0,par,varpar,nnod,nbod,npar,vels,vvari,vxnodes,&
                               vynodes,vznodes)
	   integer,intent(inout)                         :: nnod,npar
	   integer,intent(in)                            :: nbod
           real(kind=8),DIMENSION(:),intent(inout)       :: vxnodes
           real(kind=8),DIMENSION(:),intent(inout)       :: vynodes
           real(kind=8),DIMENSION(:),intent(inout)       :: vznodes
           real(kind=8),DIMENSION(:,:,:,:),intent(inout) :: vels,vvari
           real(kind=8),DIMENSION(:),intent(inout)       :: par,varpar,v0
	  END SUBROUTINE READ_MODV
	 end interface
	end module MOD_readmodv
