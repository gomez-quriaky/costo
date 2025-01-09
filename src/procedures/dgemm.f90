!=====================================================================
!=====================================================================
!
! File of the DOUBLE PRECISION Level-3 BLAS.
! ==========================================
!
!     SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
!
!     See:
!
!        Dongarra J. J.,   Du Croz J. J.,   Duff I.  and   Hammarling S.
!        A set of  Level 3  Basic Linear Algebra Subprograms.  Technical
!        Memorandum No.88 (Revision 1), Mathematics and Computer Science
!        Division,  Argonne National Laboratory, 9700 South Cass Avenue,
!        Argonne, Illinois 60439.
!
!        This subroutine was initially in F77, and rewritten in F90
!        by C. Tiberi (10/05/01).
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
! Call by BLDMAT
! Calls LSAME (contains function)
!
!=====================================================================
!=====================================================================

    SUBROUTINE DGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DGEMM
!=====================================================================

      character(len=1),intent(in)               :: TRANSA, TRANSB

      integer ,intent(in)                       :: M, N, K, LDA, LDB, LDC

      real(kind=8),intent(in)                   :: ALPHA, BETA
      real(kind=8),DIMENSION(:,:),intent(in)    :: A, B
      real(kind=8),DIMENSION(:,:),intent(inout) :: C

!=====================================================================
! Declaration of the dummy arguments of DGEMM
!=====================================================================
      logical         :: NOTA, NOTB

      integer         :: I, INFO, J, L, NCOLA, NROWA, NROWB

      real(kind=8)    :: TEMP
      real(kind=8),PARAMETER :: ONE = 1.0D+0
      real(kind=8),PARAMETER :: ZERO = 0.0D+0 
!=====================================================================
! Executable Statements
!=====================================================================

!=====================================================================
! Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
! transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
! and  columns of  A  and the  number of  rows  of  B  respectively.
!=====================================================================

      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      if( NOTA ) then
         NROWA = M
         NCOLA = K
      else
         NROWA = K
         NCOLA = M
      end if
      if( NOTB ) then
         NROWB = K
      else
         NROWB = N
      end if
!=====================================================================
! Test the input parameters.
!=====================================================================
      INFO = 0
      if((.not.NOTA).AND.(.not.LSAME(TRANSA,'C')).AND.(.not.LSAME(TRANSA,'T'))) then
         INFO = 1
      else if((.not.NOTB).AND.(.not.LSAME(TRANSB,'C')).AND.(.not.LSAME(TRANSB,'T'))) then
         INFO = 2
      else if( M  .lt.0) then
         INFO = 3
      else if( N  .lt.0) then
         INFO = 4
      else if( K  .lt.0) then
         INFO = 5
      else if( LDA.lt.MAX( 1, NROWA ) ) then
         INFO = 8
      else if( LDB.lt.MAX( 1, NROWB ) ) then
         INFO = 10
      else if( LDC.lt.MAX( 1, M     ) ) then
         INFO = 13
      end if
      if( INFO.ne.0 ) then
         write(*,'(''On entry DGEMM, parameter no'',i3,'' had illegal value'')') INFO
         STOP 'in DGEMM'
      end if

!=====================================================================
! Quick return if possible.
!=====================================================================
      if((M.eq.0).OR.(N.eq.0).OR.(((ALPHA.eq.ZERO).OR.(K.eq.0)).AND.(BETA.eq.ONE))) then
         RETURN
      end if
!=====================================================================
! And if  alpha.eq.zero.
!=====================================================================
      if( ALPHA.eq.ZERO ) then
         if( BETA.eq.ZERO ) then
            do J = 1, N
               do I = 1, M
                  C( I, J ) = ZERO
               end do
            end do
         else
            do J = 1, N
               do I = 1, M
                  C( I, J ) = BETA*C( I, J )
               end do
            end do
         end if
         RETURN
      end if
!=====================================================================
! Start the operations.
!=====================================================================
      if( NOTB ) then
         if( NOTA ) then
!=====================================================================
! Form  C := alpha*A*B + beta*C.
!=====================================================================
            do J = 1, N
               if( BETA.eq.ZERO ) then
                  do I = 1, M
                     C( I, J ) = ZERO
                  end do
               else if( BETA.ne.ONE ) then
                  do I = 1, M
                     C( I, J ) = BETA*C( I, J )
                  end do
               end if
               do L = 1, K
                  if( B( L, J ).ne.ZERO ) then
                     TEMP = ALPHA*B( L, J )
                     do I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     end do
                  end if
               end do
            end do
         else
!=====================================================================
! Form  C := alpha*A'*B + beta*C
!=====================================================================
            do J = 1, N

               do I = 1, M
                  TEMP = ZERO
                  do L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
                  end do
                  if( BETA.eq.ZERO ) then
                     C( I, J ) = ALPHA*TEMP
                  else
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  end if
               end do

            end do
         end if
      else
         if( NOTA ) then
!=====================================================================
! Form  C := alpha*A*B' + beta*C
!=====================================================================
            do J = 1, N
               if( BETA.eq.ZERO ) then
                  do I = 1, M
                     C( I, J ) = ZERO
                  end do
               else if( BETA.ne.ONE ) then
                  do I = 1, M
                     C( I, J ) = BETA*C( I, J )
                  end do
               end if
               do L = 1, K
                  if( B( J, L ).ne.ZERO ) then
                     TEMP = ALPHA*B( J, L )
                     do I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     end do
                  end if
               end do
            end do
         else
!=====================================================================
! Form  C := alpha*A'*B' + beta*C
!=====================================================================
            do J = 1, N
               do I = 1, M
                  TEMP = ZERO
                  do L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
                  end do
                  if( BETA.EQ.ZERO ) then
                     C( I, J ) = ALPHA*TEMP
                  else
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  end if
               end do
            end do
         end if
      end if

    CONTAINS
!=====================================================================
!=====================================================================
!
! LOGICAL FUNCTION LSAME ( CA, CB )
!
! Purpose
!  =======
!
!  LSAME  tests if CA is the same letter as CB regardless of case.
!  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
!  CA is either the same as CB or the equivalent lower case letter.
!
!  N.B. This version of the routine is only correct for ASCII code.
!       Installers must modify the routine for other character-codes.
!
!       For EBCDIC systems the constant IOFF must be changed to -64.
!       For CDC systems using 6-12 bit representations, the system-
!       specific code in comments must be activated.
!
!  Parameters
!  ==========
!
!  CA     - CHARACTER*1
!  CB     - CHARACTER*1
!           On entry, CA and CB specify characters to be compared.
!           Unchanged on exit.
!
!  Auxiliary routine for Level 2 Blas.
!
!  -- Written on 20-July-1986
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, Nag Central Office.
!
!=====================================================================
!=====================================================================
      LOGICAL FUNCTION LSAME ( CA, CB )

        IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of LSAME
!=====================================================================

        character(len=1),intent(in)  :: CA, CB

!=====================================================================
! Declaration of the dummy arguments of LSAME
!=====================================================================
        integer,PARAMETER            :: ioff=32

        LSAME = CA .eq. CB
!=====================================================================
! Now test for equivalence
!=====================================================================
        if ( .not.LSAME ) then
           LSAME = ICHAR(CA) - IOFF .eq. ICHAR(CB)
        end if

      END FUNCTION LSAME

!=====================================================================
! End of DGEMM .
!=====================================================================
      
    END SUBROUTINE DGEMM
