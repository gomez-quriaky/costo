!=====================================================================
! Module for explicite interface of subroutine DGEMV
!=====================================================================
	module MOD_dgemv
	 interface
	  SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
	   character(len=1),intent(in)            :: TRANS
           integer,intent(in)                     :: INCX, INCY, LDA, M, N
           real(kind=8),intent(in)                :: ALPHA, BETA
           real(kind=8),DIMENSION(:),intent(in)   :: X
           real(kind=8),DIMENSION(:),intent(out)  :: Y
           real(kind=8),DIMENSION(:,:),intent(in) :: A
	  END SUBROUTINE DGEMV
	 end interface
	end module MOD_dgemv
