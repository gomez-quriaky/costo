!=====================================================================
!=====================================================================
!
!			SUBROUTINE INVERMAT
!
!> This subroutine inverts the matrix BDERI previously computed in
!! subroutine BLDMAT.
!!
!! WARNING: THIS INVERSION IS WORKING ONLY IF BDERI IS SYMMETRIC !!!!
!!
!! First, try a simple inversion of the matrix. If failure, then
!! a LU (Lower-Upper triangular part) decomposition is proceeded.
!! The output is the matrix BDERI.
!
! Calls : DSINV(DMFSD),LUDCMP,LUBKSB
! called by MAIN
!
!=====================================================================
!=====================================================================
    !INCLUDE 'INTERF/MOD_dsinv.f'
    !INCLUDE 'INTERF/MOD_ludcmp.f'
    !INCLUDE 'INTERF/MOD_lubksb.f'

    SUBROUTINE INVERMAT(bderi,npar,iiter,regul,lambda)

      USE MOD_unit

      USE MOD_dsinv

      !USE MOD_ludcmp
      !USE MOD_lubksb
      USE LAPACK95
      USE F95_PRECISION
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of INVERMAT
!=====================================================================
      integer,intent(in)                          :: npar,iiter,regul

      real(kind=8),intent(in)			               :: lambda
      real(kind=8),DIMENSION(:,:),intent(inout)   :: bderi
!=====================================================================
! Declaration of the dummy arguments of INVERMAT
!=====================================================================
      integer                                  :: ier,i,j,jj
      integer,DIMENSION(npar)                  :: indx

      real(kind=8)                             :: toler,d
!      real(kind=8),DIMENSION(npar,npar)        :: matemp,dcmp
      real(kind=8),DIMENSION(:,:),ALLOCATABLE  :: dcmp
      real(kind=8),DIMENSION(npar)             :: vv
      
      integer                                   :: status
!=====================================================================
! Initialization of some variables
! DCMP-matrix is first the identity matrix
!=====================================================================
      ier=0
      toler=1.0d-7
      status = 0
!      matemp(:,:) = bderi(:,:)

      write(*,*)''
      write(*,*)'INVERTING THE MATRIX'
      write(*,*)'    Iteration no ',iiter
      write(*,*)''

      write(inout,*)''
      write(inout,*)'INVERTING THE MATRIX'
      write(inout,*)'    Iteration no ',iiter
      write(inout,*)''
!=====================================================================
! Multiply by lambda the diagonal terms to be sure the matrix is able 
! to be inversed
!=====================================================================
      if(regul.eq.1.and.iiter.eq.1) then
	      do i=1,npar
 	         bderi(i,i)=bderi(i,i)*lambda
         end do
      end if
!=====================================================================
! Try a simple invertion of bderi matrix.
! If failure, proceed to the LU decomposition, bderi unchanged
!=====================================================================
      !CALL DSINV(bderi,npar,toler,ier)
      !call dpotrf( uplo, n, a, lda, info )
      
      CALL dpotrf('L', npar, bderi, npar , status)

      print *, 'result factorizaton Cholesky ', status
!=====================================================================
! If no problem encountered during DSINV, add the symmetrical part of 
! the matrix
!=====================================================================
      !if(ier.eq.0) then
      if(status.eq.0) then

         do i=2,npar
            jj=i-1
            do j=1,jj
               bderi(i,j)=bderi(j,i)
            end do
         end do

      else
!=====================================================================
! Start the Lower-Upper triangular matrix decomposition (Working even 
! for non-symmetric matrices)
! MATEMP matrix is then its LU decomposition 
! Bderi matrix has been unchanged by the DSINV subroutine...
!=====================================================================
         write(*,*)''
         write(*,*)'    Problem in inverting matrix, in INVERMAT ier= ',ier
         write(*,*)'    Have to proceed to LU decomposition..........'

         write(inout,*)''
         write(inout,*)'    Problem in inverting matrix, in INVERMAT ier= ',ier
         write(inout,*)'    Have to proceed to LU decomposition..........'

         ALLOCATE(dcmp(npar,npar))
         
!         CALL LUDCMP(matemp,indx,d)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !ADAPT to MKL SPARSE BLAS
! 	 CALL LUDCMP(bderi,indx,d)
!          write(*,*)''
!          write(*,*)'    LU decomposition done.'

!          dcmp(:,:)=0.d0
!          do i=1,npar
!             dcmp(i,i)=1.
!          end do

!          do j=1,npar
!             do i=1,npar
!                vv(i)=dcmp(i,j)
!             end do
! !            CALL LUBKSB(matemp,indx,vv)
! 	    CALL LUBKSB(bderi,indx,vv)
!             do i=1,npar
!                dcmp(i,j)=vv(i)
!             end do
!          end do
!          write(*,*)''
!          write(*,*)'    LU Back Substitution done.'

!          bderi(:,:)=dcmp(:,:)
	 DEALLOCATE(dcmp)

      end if

!=====================================================================
! Storage of the diagonal elements of BDERI through the iterations
! so that we can retreive the covariance of each parameters
!=====================================================================



    END SUBROUTINE INVERMAT

!=====================================================================
!=====================================================================
!
!			SUBROUTINE DSINV
!
! This subroutine inverts the matrix BDERI previously computed in
! subroutine BLDMAT
!
! Calls : DMFSD
! called by INVERMAT
!
!=====================================================================
!=====================================================================
    !INCLUDE 'INTERF/MOD_dmfsd.f'

   SUBROUTINE DSINV(bderi,npar,toler,ier)

      USE MOD_unit

      
      USE MOD_dmfsd

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DSINV
!=====================================================================
      integer,intent(in)                         :: npar
      integer,intent(inout)                      :: ier

      real(kind=8),intent(in)                    :: toler
      real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
!=====================================================================
! Declaration of the dummy arguments of DSINV
!=====================================================================
      integer                            :: ka,ia,ja,i,j,k,ipiv,min,la
      integer                            :: kend,lanf,ind,lhor,lver,ka1

      real(kind=8)                       :: din,work
      real(kind=8),DIMENSION(npar*npar)  :: mat
!=====================================================================
! Beginning of REAL BORING mathematica...
! if anyone understands...
! Is it the Cholesky decomposition?????
!=====================================================================
      mat=reshape(bderi,(/npar*npar/))
!=====================================================================
! First, it takes the upper triangular part of the MAT-matrix
! and store it in mat
!=====================================================================
      ka=0
      do ia=1,npar
         do ja=1,ia
            ka1=npar*(ia-1)+ja
            mat(ka+ja)=mat(ka1)
         end do
         ka=ka+ia
      end do

      CALL DMFSD(mat,npar,toler,ier)

!=====================================================================
! If ier < 0 simple inversion is impossible, proceed to a LU 
! decomposition
! If ier >= 0 , proceed to the inversion
!=====================================================================
      if(ier.lt.0) then
         write(*,*)''
         write(*,*)'    IER.lt.0, stop the DSINV subroutine and proceed to &
              &a LU decomposition.'

         write(inout,*)''
         write(inout,*)'    IER.lt.0, stop the DSINV subroutine and proceed &
              &to a LU decomposition.'
      else
         ipiv=npar*(npar+1)/2
         ind=ipiv
         do i=1,npar
            din=1./mat(ipiv)
            mat(ipiv)=din
            min=npar
            kend=i-1
            lanf=npar-kend
            
            if(kend.gt.0) then
               j=ind
               do k=1,kend
                  work=0.d0
                  min=min-1
                  lhor=ipiv
                  lver=j
                  do la=lanf,min
                     lver=lver+1
                     lhor=lhor+la
                     work=work+mat(lver)*mat(lhor)
                  end do
                  mat(j)=-work*din
                  j=j-min
               end do
            end if
            ipiv=ipiv-min
            ind=ind-1
         end do
      
         do i=1,npar
            ipiv=ipiv+i
            j=ipiv
            do k=i,npar
               work=0.d0
               lhor=j
               do la=k,npar
                  lver=lhor+k-i
                  work=work+mat(lhor)*mat(lver)
                  lhor=lhor+la
               end do
               mat(j)=work
               j=j+k
            end do
         end do

         ka=ipiv
         do ia=npar,1,-1
            ka=ka-ia
            do ja=ia,1,-1
               ka1=npar*(ia-1)+ja
               mat(ka1)=mat(ka+ja)
            end do
         end do

         bderi=reshape(mat,(/npar,npar/))

      end if

    END SUBROUTINE DSINV

!=====================================================================
!=====================================================================
!
!			SUBROUTINE DMFSD
!
! I don't know exactly what this routine does... May be Cholsky 
! decomposition ?
!
! Calls : none
! called by DSINV
!
!=====================================================================
!=====================================================================

    SUBROUTINE DMFSD(mat,npar,toler,ier)

      USE MOD_unit

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DMFSD
!=====================================================================
      integer,intent(in)                         :: npar
      integer,intent(inout)                      :: ier

      real(kind=8),intent(in)                    :: toler
      real(kind=8),DIMENSION(:),intent(inout)    :: mat
!=====================================================================
! Declaration of the dummy arguments of DMFSD
!=====================================================================
      integer               :: kpiv,k,ind,lend,i,la,lanf,lind

      real(kind=8)          :: tol,dsum,dpiv
!=====================================================================
! If npar=1, return an error message
!=====================================================================
      if((npar-1).lt.0) then
         ier=-1
         write(*,*)''
         write(*,*)'    Number of parameter is 1. IER=-1 from DMFSD'

         write(inout,*)''
         write(inout,*)'    Number of parameter is 1. IER=-1 from DMFSD'
      else
         ier=0
         kpiv=0
!=====================================================================
! Loop over the columns
!=====================================================================
kloop:   do k=1,npar
            kpiv=kpiv+k
            ind=kpiv
            lend=k-1
            tol=ABS(toler*mat(kpiv))
            
            do i=k,npar
               dsum=0.d0

               if(lend.ne.0) then
                  do la=1,lend
                     lanf=kpiv-la
                     lind=ind-la
                     dsum=dsum+mat(lanf)*mat(lind)
                  end do
               end if
               dsum=mat(ind)-dsum

               if((i-k).eq.0) then
                  
                  if((dsum-tol).le.0) then
                     if(dsum.le.0) then
                        ier=-1
                        write(*,*)'    Problems in DMFSD IER return =-1'
                        write(inout,*)'    Problems in DMFSD IER return =-1'
                        exit kloop
                     else
                        if(ier.le.0) then
                           ier=k-1
                           if(ier.ne.0) then
                              write(*,*)'    Problems in DMFSD IER = ',ier
                              write(inout,*)'    Problems in DMFSD IER = ',ier
                           end if
                        end if
                        dpiv=dsqrt(dsum)
                        mat(kpiv)=dpiv
                        dpiv=1./dpiv
                        ind=ind+i
                        cycle
                     end if
                  else
                     dpiv=dsqrt(dsum)
                     mat(kpiv)=dpiv
                     dpiv=1./dpiv
                     ind=ind+i
                     cycle
                  end if
               else
                  mat(ind)=dsum*dpiv
                  ind=ind+i
               end if

            end do

         end do kloop

      end if

    END SUBROUTINE DMFSD
