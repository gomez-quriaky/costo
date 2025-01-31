!=====================================================================
! Module for explicite interface of subroutine CONTBOD
!=====================================================================
	module MOD_contbod
	 interface
	  SUBROUTINE CONTBOD(xb,yb,nbod,ilay,iside,jside)
	   integer,intent(in)                     :: nbod
           integer,DIMENSION(:),intent(in)        :: ilay
           integer,DIMENSION(:),intent(inout)     :: iside
           integer,DIMENSION(:,:),intent(inout)   :: jside
           real(kind=8),DIMENSION(:,:),intent(in) :: xb,yb
	  END SUBROUTINE CONTBOD
	 end interface
	end module MOD_contbod
