!=====================================================================
! Module for explicite interface of subroutine GMATRIX
!=====================================================================
	module MOD_gmatrix
	 interface
	  SUBROUTINE GMATRIX(G,nr,dtm)
	   integer,intent(in)                       :: nr
           real(kind=8),DIMENSION(:,:,:),intent(in) :: dtm
           real(kind=8),DIMENSION(:,:),intent(inout):: G
	  END SUBROUTINE GMATRIX
	 end interface
	end module MOD_gmatrix
