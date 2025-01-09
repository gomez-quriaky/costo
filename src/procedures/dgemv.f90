!=====================================================================
!=====================================================================
!
!     File of the DOUBLE PRECISION  Level-2 BLAS.
!     ===========================================
!
!     SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
!
!     See:
!
!        Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
!        An  extended  set of Fortran  Basic Linear Algebra Subprograms.
!
!        Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
!        and  Computer Science  Division,  Argonne  National Laboratory,
!        9700 South Cass Avenue, Argonne, Illinois 60439, US.
!
!        Or
!
!        NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
!        Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
!        OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
!        Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.
!
! Purpose
! =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
! Parameters
! ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
! Call by BLDMAT, PERTURB
! Calls LSAME (contains function)
!=====================================================================
!=====================================================================

    SUBROUTINE DGEMV (TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of DGEMV
!=====================================================================
      character(len=1),intent(in)            :: TRANS

      integer,intent(in)                     :: INCX, INCY, LDA, M, N

      real(kind=8),intent(in)                :: ALPHA, BETA
      real(kind=8),DIMENSION(:),intent(in)   :: X
      real(kind=8),DIMENSION(:),intent(out)  :: Y
      real(kind=8),DIMENSION(:,:),intent(in) :: A
!=====================================================================
! Declaration of the dummy arguments of DGEMV
!=====================================================================
      integer                  :: I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY

      real(kind=8)             :: TEMP
      real(kind=8),PARAMETER   :: ONE = 1.0D+0
      real(kind=8),PARAMETER   :: ZERO = 0.0D+0
!=====================================================================
! Executable Statements
!=====================================================================

!=====================================================================
! Test the input parameters
!=====================================================================
      INFO = 0
      if ( .not.LSAME( TRANS, 'N' ).AND. .not.LSAME( TRANS, 'T' ).AND.&
           .not.LSAME( TRANS, 'C' )) then
         INFO = 1
      else if( M.lt.0 ) then
         INFO = 2
      else if( N.lt.0 ) then
         INFO = 3
      else if( LDA.lt.MAX( 1, M ) ) then
         INFO = 6
      else if( INCX.eq.0 ) then
         INFO = 8
      else if( INCY.eq.0 ) then
         INFO = 11
      end if
      if( INFO.ne.0 ) then
         write(*,'(''On entry DGEMM, parameter no'',i3,'' had illegal value'')') INFO
         STOP 'in DGEMM'
      end if
!=====================================================================
!     Quick return if possible.
!=====================================================================
      if((M.eq.0).OR.(N.eq.0).OR.((ALPHA.eq.ZERO).AND.(BETA.eq.ONE))) then
         RETURN
      end if
!=====================================================================
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!=====================================================================
      if( LSAME( TRANS, 'N' ) ) then
         LENX = N
         LENY = M
      else
         LENX = M
         LENY = N
      end if
      if( INCX.gt.0 ) then
         KX = 1
      else
         KX = 1 - ( LENX - 1 )*INCX
      end if
      if( INCY.gt.0 ) then
         KY = 1
      else
         KY = 1 - ( LENY - 1 )*INCY
      end if
!=====================================================================
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!=====================================================================
      if( BETA.ne.ONE ) then
         if( INCY.eq.1 ) then
            if( BETA.eq.ZERO ) then
               do I = 1, LENY
                  Y( I ) = ZERO
               end do
            else
               do I = 1, LENY
                  Y( I ) = BETA*Y( I )
               end do
            end if
         else
            IY = KY
            if( BETA.eq.ZERO ) then
               do I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
               end do
            else
               do I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
               end do
            end if
         end if
      end if
      if( ALPHA.eq.ZERO ) then
         RETURN
      end if
      if( LSAME( TRANS, 'N' ) ) then
!=====================================================================
!        Form  y := alpha*A*x + y.
!=====================================================================
         JX = KX
         if( INCY.eq.1 ) then
            do J = 1, N
               if( X( JX ).ne.ZERO ) then
                  TEMP = ALPHA*X( JX )
                  do I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
                  end do
               end if
               JX = JX + INCX
            end do
         else
            do J = 1, N
               if( X( JX ).ne.ZERO ) then
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  do I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
                  end do
               end if
               JX = JX + INCX
            end do
         end if
      else
!=====================================================================
!        Form  y := alpha*A'*x + y.
!=====================================================================
         JY = KY
         if( INCX.eq.1 ) then
            do J = 1, N
               TEMP = ZERO
               do I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
               end do
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
            end do
         else
            do J = 1, N
               TEMP = ZERO
               IX   = KX
               do I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
               end do
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
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
!     End of DGEMV .
!=====================================================================

    END SUBROUTINE DGEMV
