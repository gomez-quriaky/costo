!=====================================================================
! Module for explicite interface of subroutine BLOCKMEAN
!=====================================================================
	module MOD_blockmean
	 interface
	  SUBROUTINE BLOCKMEAN(xb1,xb2,yb1,yb2,zb1,zb2,Vm,velco)
	   real(kind=8),intent(out)              :: Vm
           real(kind=8),intent(in)               :: xb1,xb2,yb1,yb2,zb1,zb2
           real(kind=8),DIMENSION(:),intent(in)  :: velco
	  END SUBROUTINE BLOCKMEAN
	 end interface
	end module MOD_blockmean
