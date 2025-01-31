!=====================================================================
! Module for explicite interface of subroutine delaz
!=====================================================================
	module MOD_delaz
	 interface
          SUBROUTINE DELAZ(lat,lon,x,y,ct0,st0,phi0,solat)
	   real(kind=8),intent(in)      :: lat,lon,ct0,st0,phi0,solat
           real(kind=8),intent(out)     :: x,y
	  END SUBROUTINE DELAZ
	 end interface
	end module MOD_delaz
