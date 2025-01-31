!=====================================================================
!=====================================================================
!
!			SUBROUTINE LUBKSB
!                    (LU BacK SuBsititution)
!
! Taken From Numerical Recipes in FROTRAN90.
!
! Solves the set of N equations A.X=B. Here the NxN matrix A is input,
! not as the original matrix A, but rather as its LU decomposition,
! determined by the routine LUDCMP. INDX is input as the permutation
! vector length N returned by LUDCMP. B is input as the right-hand-side
! vector B, also of length N, and returns with the solution vector X.
! A and INDX are not modified by this routine and can be left in place
! for successive calls with different right-hand-sides B.
! This routine takes into account the possibility that B will begin with
! many zero elements, so it is efficient for use in matrix inversion.
!
! Modification for two last lines of the original programm: the principal
! loop runs from i=1 to N, but the calculation implicates i+1, which
! produces an out of range error.
!
! Other bugs (22/05/01) concerning the last line : the division by
! a(i,i) should be done after the end if.
!
! Calls none
! called by INVMAT
!=====================================================================
!=====================================================================

    SUBROUTINE LUBKSB(a,indx,b)

      USE mynrutil, ONLY : assert_eq

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of LUBKSB
!=====================================================================
      integer,DIMENSION(:),intent(in)         :: indx

      real(kind=8),DIMENSION(:,:),intent(in)  :: a
      real(kind=8),DIMENSION(:),intent(inout) :: b
!=====================================================================
! Declaration of the dummy arguments of LUBKSB
!=====================================================================
      integer           :: i,n,ii,ll
      real(kind=8)      :: summ
!=====================================================================
! When ii is set to a positive value, it will become the index of the
! first nonvanishing element of B.
!=====================================================================
      n=assert_eq(size(a,1),size(a,2),size(indx),'LUBKSB')
      ii=0
!=====================================================================
! We now do the forward substitution, equation (2.3.6).
! The only new wrinkle is to unscramble the permutation as we go.
!=====================================================================
      do i=1,n
         ll=indx(i)
         summ=b(ll)
         b(ll)=b(i)

         if (ii /= 0) then
            summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
!=====================================================================
! else, if a nonzero element was encountered, from now on we will have
! to do the dot product above
!=====================================================================
         else if (summ /= 0.0) then
            ii=i
         end if

         b(i)=summ
      end do
!=====================================================================
! Now we do the back-substitution (equation 2.3.7)
!=====================================================================
      do i=n,1,-1
         if (i .lt. n) then
            b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))
         end if
         b(i) = b(i)/a(i,i)
      end do

    END SUBROUTINE LUBKSB
