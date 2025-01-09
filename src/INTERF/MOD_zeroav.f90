!=====================================================================
! Module for explicite interface of subroutine zeroav
!=====================================================================
	module MOD_zeroav
	 interface
          SUBROUTINE ZEROAV(x_block,y_block,par,drmax,ibove)
	   integer,DIMENSION(:),intent(in)         :: x_block,y_block
	   integer,DIMENSION(:),intent(in)         :: ibove
           real(kind=8),DIMENSION(:),intent(in)    :: drmax
           real(kind=8),DIMENSION(:),intent(inout) :: par
	  END SUBROUTINE ZEROAV
	 end interface
	end module MOD_zeroav
