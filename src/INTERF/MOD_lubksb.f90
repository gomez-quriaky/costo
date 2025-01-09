!=====================================================================
! Module for explicite interface of subroutine LUBKSB
!=====================================================================
	module MOD_lubksb
	 interface
	  SUBROUTINE LUBKSB(a,indx,b)
	   integer,DIMENSION(:),intent(in)         :: indx
           real(kind=8),DIMENSION(:,:),intent(in)  :: a
           real(kind=8),DIMENSION(:),intent(inout) :: b
	  END SUBROUTINE LUBKSB
	 end interface
	end module MOD_lubksb
