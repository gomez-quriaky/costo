!=====================================================================
! Module for explicite interface of subroutine vcoefs
!=====================================================================
	module MOD_vcoefs
	 interface
          SUBROUTINE VCOEFS(velco,iiter,vels,vxnodes,vynodes,vznodes)
	   integer,intent(in)                          :: iiter
	   real(kind=8),DIMENSION(:),intent(out)       :: velco
	   real(kind=8),DIMENSION(:,:,:,:),intent(in)  :: vels
	   real(kind=8),DIMENSION(:),intent(in)        :: vxnodes,vynodes,vznodes
	  END SUBROUTINE VCOEFS
	 end interface
	end module MOD_vcoefs
