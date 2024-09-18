!=====================================================================
! Module for explicite interface of subroutine SIMPLEX
!=====================================================================
	module MOD_simplex
	 interface
	  SUBROUTINE SIMPLEX(nrp,fstime,rp,velco)
	   integer,intent(in)                        :: nrp
	   real(kind=8),intent(inout)                :: fstime
	   real(kind=8),DIMENSION(:,:),intent(inout) :: rp
	   real(kind=8),DIMENSION(:),intent(in)      :: velco
	  END SUBROUTINE SIMPLEX
	 end interface
	end module MOD_simplex
