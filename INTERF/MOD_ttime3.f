!=====================================================================
! Module for explicite interface of subroutine TTIME3
!=====================================================================
	module MOD_ttime3
	 interface
	  SUBROUTINE TTIME3(nrp,pth,tt,dpath,velco,vznodes)
	   integer,intent(in)                     :: nrp
	   real(kind=8),intent(out)               :: tt,dpath
           real(kind=8),DIMENSION(:,:),intent(in) :: pth
           real(kind=8),DIMENSION(:),intent(in)   :: velco,vznodes
	  END SUBROUTINE TTIME3
	 end interface
	end module MOD_ttime3
