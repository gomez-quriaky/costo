!=====================================================================
! Module for explicite interface of subroutine SUMTA2
!=====================================================================
	module MOD_sumta2
	 interface
	  SUBROUTINE SUMTA2(x,y1,y2,isign,sum)
	   integer,intent(in)                 :: isign
           real(kind=8),intent(in)            :: y1,y2
           real(kind=8),intent(inout)         :: sum,x
	  END SUBROUTINE SUMTA2
	 end interface
	end module MOD_sumta2
