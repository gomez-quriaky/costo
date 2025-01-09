!=====================================================================
! Module for explicite interface of subroutine VEL3
!=====================================================================
	module MOD_vel3
	 interface
	  SUBROUTINE VEL3(x,y,z,v,velco,ip,jp,kp)
	   integer,intent(out)                   :: ip,jp,kp
           real(kind=8),intent(in)               :: x,y,z
           real(kind=8),intent(out)              :: v
           real(kind=8),DIMENSION(:),intent(in)  :: velco
	  END SUBROUTINE VEL3
	 end interface
	end module MOD_vel3
