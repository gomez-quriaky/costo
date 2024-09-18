!=====================================================================
! Module for explicite interface of subroutine DGEMM
!=====================================================================
	module MOD_dgemm
	 interface
	  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
	   character(len=1),intent(in)               :: TRANSA, TRANSB
           integer ,intent(in)                       :: M, N, K, LDA, LDB, LDC
           real(kind=8),intent(in)                   :: ALPHA, BETA
           real(kind=8),DIMENSION(:,:),intent(in)    :: A, B
           real(kind=8),DIMENSION(:,:),intent(inout) :: C
	  END SUBROUTINE DGEMM
	 end interface
	end module MOD_dgemm
