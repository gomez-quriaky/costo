!=====================================================================
! Module for explicite interface of subroutine RAYWEB
!=====================================================================
	module MOD_rayweb
	 interface
	  SUBROUTINE RAYWEB(xe,ye,ze,xs,ys,zs,fstime,rp,nrp,velco,nr)
	   integer,intent(out)                :: nrp
      	   integer,intent(in)                 :: nr
      	   real(kind=8),intent(in)            :: xs,ys,zs,xe,ye,ze
      	   real(kind=8),intent(out)           :: fstime
      	   real(kind=8),DIMENSION(:,:),intent(inout) :: rp
      	   real(kind=8),DIMENSION(:),intent(in)    :: velco
	  END SUBROUTINE RAYWEB
	 end interface
	end module MOD_rayweb
