!=====================================================================
! Module for explicite interface of subroutine SORT
!=====================================================================
	module MOD_sort
	 interface
	  SUBROUTINE SORT(a,b,c,tt1,tt2,tt3,nrp)
	   integer,intent(in)                        :: nrp
	   real(kind=8),intent(inout)                :: tt1,tt2,tt3
	   real(kind=8),DIMENSION(:,:),intent(inout) :: a,b,c
	  END SUBROUTINE SORT
	 end interface
	end module MOD_sort
