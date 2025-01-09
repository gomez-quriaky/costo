!=====================================================================
! Module for explicite interface of subroutine TTIME2
!=====================================================================
	module MOD_ttime2
	 interface
	  SUBROUTINE TTIME2(nrp,pth,tt,velco)
	   integer,intent(in)                     :: nrp
	   real(kind=8),intent(out)               :: tt
           real(kind=8),DIMENSION(:,:),intent(in) :: pth
           real(kind=8),DIMENSION(:),intent(in)   :: velco
	  END SUBROUTINE TTIME2
	 end interface
	end module MOD_ttime2
