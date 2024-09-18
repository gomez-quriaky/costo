!=====================================================================
! Module for explicite interface of subroutine LUDCMP
!=====================================================================
	module MOD_ludcmp
	 interface
	  SUBROUTINE LUDCMP(a,indx,d)
	   integer,DIMENSION(:),intent(out)              :: indx
           real(kind=8),DIMENSION(:,:),intent(inout)     :: a
           real(kind=8),intent(out)                      :: d
	  END SUBROUTINE LUDCMP
	 end interface
	end module MOD_ludcmp
