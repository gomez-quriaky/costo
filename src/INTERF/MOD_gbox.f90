!=====================================================================
! Module for explicite interface of subroutine GBOX
!=====================================================================
	module MOD_gbox
	 interface
	  SUBROUTINE GBOX(gFX,gFY,gFZ,xb1,yb1,zb1,xb2,yb2,zb2,rhoinit,&
			  g,sum)
	   real(kind=8),intent(in)              :: rhoinit
           real(kind=8),intent(in)              :: gFX,gFY,gFZ
	   real(kind=8),intent(in)              :: xb1,yb1,zb1,xb2,yb2,zb2
           real(kind=8),intent(out)             :: g,sum
	  END SUBROUTINE GBOX
	 end interface
	end module MOD_gbox
