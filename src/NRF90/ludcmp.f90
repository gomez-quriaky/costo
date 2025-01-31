!=====================================================================
!=====================================================================
!
!			SUBROUTINE LUDCMP
!
! Given a NxN input matrix A, this subroutine replaces it by the LU 
! decomposition of a rowwise permutation of itself. On output, A is
! arranged as in equation (2.3.14) of Numerical Recipes; INDX is an
! output vector of length N that records the row permutation effected
! by the partial pivoting.
! D is output as +/- 1 depending of wether the number of row interchanges
! was even or odd, respectively.
! This routine is used in combinaison with LUBKSB to solve linera equations
! or invert matrix
!
! Calls none
! called by INVMAT
!=====================================================================
!=====================================================================

      SUBROUTINE LUDCMP(a,indx,d)

      USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of LUDCMP
!=====================================================================
      integer,DIMENSION(:),intent(out)          :: indx

      real(kind=8),DIMENSION(:,:),intent(inout) :: a
      real(kind=8),intent(out)                  :: d
!=====================================================================
! Declaration of the dummy arguments of LUDCMP
!=====================================================================
      integer                           :: j,n,imax
      real(kind=8),DIMENSION(size(a,1)) :: vv
      real(kind=8),PARAMETER            :: TINY=1.0e-20_sp
!=====================================================================
! vv stores the implicit scaling of each row
! TINY is a small number
!=====================================================================
      n=assert_eq(size(a,1),size(a,2),size(indx),'LUDCMP')
      d=1.0
      vv=maxval(abs(a),dim=2)

      if (any(vv == 0.0)) call nrerror('singular matrix in LUDCMP')
!=====================================================================
! Save the scaling in vv
!=====================================================================
      vv=1.0_sp/vv

      do j=1,n
!=====================================================================
! Find the pivot row
!=====================================================================
         imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
!=====================================================================
! If we need to interchange rows, do so:
!    - change the parity of D
!    - interchange the scaling factor vv
!=====================================================================
         if (j /= imax) then
            call swap(a(imax,:),a(j,:))
            d=-d
            vv(imax)=vv(j)
         end if

         indx(j)=imax
!=====================================================================
! If the pivot element is zero, the matrix is singular (at least to
! the precision of the algorithm). For some applications on
! singular matrices, it is desirable to substitute TINY for zero
!=====================================================================
         if (a(j,j) == 0.0) a(j,j)=TINY
!=====================================================================
! Divide the pivot element
!=====================================================================
         a(j+1:n,j)=a(j+1:n,j)/a(j,j)
!=====================================================================
! Reduce the remaining submatrix
!=====================================================================
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))

      end do


      END SUBROUTINE LUDCMP
