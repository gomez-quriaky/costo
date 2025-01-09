!=====================================================================
!=====================================================================
!
!			SUBROUTINE READBCO
!
!=====================================================================
!=====================================================================
!> reads in file parameter.inp the initial values of B-coeff for
!! linear relation velocity/density
!!
!! puts the values in the array __PAR__ and __VARPAR__
!! \param par is values of B-coef
!! \param varpar is covariance of B-coef
!=====================================================================
!
! Calls none
! Called by MAIN
!=====================================================================

    SUBROUTINE READBCO(binit,varb,par,varpar,npar)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_inv
      USE MOD_delim

      IMPLICIT NONE
!=====================================================================
! Declaration of in/out/inout arguments of READBCO
!=====================================================================
      integer,intent(inout)                      :: npar

      real(kind=8),DIMENSION(:),intent(in)       :: binit,varb
      real(kind=8),DIMENSION(:),intent(inout)    :: par,varpar
!=====================================================================
! Declaration of dummy arguments of READBCO
!=====================================================================
      integer          :: i
!=====================================================================
! input of the B-value
!=====================================================================

      write(inout,*)' '
      write(inout,*)'READING B-COEFFICIENT MODEL'
      write(*,*)''
      write(*,*)'READING B-COEFFICIENT MODEL'

!=====================================================================
! Case of B inverted or density AND velocity are inverted
! In this case, we use the Bvalue indicated in parameter.inp file
!=====================================================================
      write(inout,'(5x,''B is inverted, reading the informa&
           &tion from parameter.inp file'')')
      write(inout,'(5x,i2,'' fixed coefficients read'')') nlayer

      do i=1,nlayer
         npar=npar+1
         par(npar)=binit(i)
         varpar(npar)=1./(varb(i)*varb(i))
         write(inout,'(7x,''Layer '',i2,'': B-coeff= '',f10.2,&
              &'' Covariance= '',f10.2)') i,par(npar),varpar(npar)
      end do
      
      ibegin(3)=iend(2)+1
      iend(3)=iend(2)+nlayer
      
      write(inout,'(5x,''B coeficients stored in arrays from i='',&
           &i6,'' to i='',i6)') ibegin(3),iend(3)
      write(*,'(5x,''B coeficients stored in arrays from i='',&
           &i6,'' to i='',i6)') ibegin(3),iend(3)
       


      END SUBROUTINE READBCO
