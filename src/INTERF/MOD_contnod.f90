!=====================================================================
! Module for explicite interface of subroutine CONTNOD
!=====================================================================
	module MOD_contnod
	 interface
	  SUBROUTINE CONTNOD(ibove,ivside,jvside)
           integer,DIMENSION(:),intent(in)        :: ibove
           integer,DIMENSION(:),intent(inout)     :: ivside
           integer,DIMENSION(:,:),intent(inout)   :: jvside
	  END SUBROUTINE CONTNOD
	 end interface
	end module MOD_contnod
