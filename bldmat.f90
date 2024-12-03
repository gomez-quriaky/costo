!=====================================================================
!=====================================================================
!
!			SUBROUTINE BLDMAT
!
!> subroutine calculates the matrix to invert in case of INV=.TRUE.
!!
!! The matrix is (At*Cd**(-1)*A + Cp**(-1) - Cs**(-1)*Ds - Cb**(-1)*Db)
!! with:
!!    - A the matrix of partial derivatives (ADERI) and At its transpose
!!    - Cd the covariance matrix of the data
!!    - Cp the covariance matrix of the parameters
!!    - Cb**(-1)*Db the linear relation term
!!    - Cs**(-1)*Ds the smoothing parameters
!!
!! It also computes the other term needed for the inversion H1:
!!
!! At*Cd**(-1)*(obs-calc) - Cs**(-1)*Smooth - Cb**(-1)*brel
!!
!! with:
!!    - (obs-calc) difference between observed and calculated data
!!    - Smooth smoothing constraint
!!    - brel linear relation constraint
!
!
! Calls: DSMOOTH,VSMOOTH,DENVEL,DGEMM,DGEMV
! called by MAIN
!
!=====================================================================
!=====================================================================
      INCLUDE 'INTERF/MOD_denvel.f'
      INCLUDE 'INTERF/MOD_dgemm.f'
      INCLUDE 'INTERF/MOD_dgemv.f'
      
      SUBROUTINE BLDMAT(iiter,aderi,bderi,punvar,varpar,h1,diff,npar,&
                        npar1,ndat,smooth,iside,jside,ivside,jvside,&
                        xb,yb,vxnodes,vynodes,par,ibove,ismooth,ilay,&
                        ddvr,nbod)
        
      USE MOD_delim
      USE MOD_unit
      USE MOD_inv

      USE MOD_dsmooth
      USE MOD_vsmooth
      USE MOD_denvel
      USE MOD_dgemm
      USE MOD_dgemv

      !USE BLAS95
      !USE F95_precision

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of BLDMAT
!=====================================================================
      integer,intent(in)                         :: iiter,npar,npar1,ndat
      integer,intent(in)                         :: ismooth,nbod
      integer,DIMENSION(:),intent(in)            :: iside,ivside,ibove,ilay
      integer,DIMENSION(:,:),intent(in)          :: jside,jvside

      real(kind=8),DIMENSION(:),intent(in)       :: diff,smooth,par,varpar
      real(kind=8),DIMENSION(:),intent(in)       :: ddvr
      real(kind=8),DIMENSION(:),intent(in)       :: vxnodes,vynodes
      real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb
      real(kind=8),DIMENSION(:,:),intent(in)     :: aderi
      real(kind=8),DIMENSION(:),intent(inout)    :: h1
      real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
      real(kind=8),DIMENSION(:),pointer          :: punvar
!=====================================================================
! Declaration of the dummy arguments of BLDMAT
!=====================================================================
      integer                                    :: i,ii,j,k,l,dcont,vcont
      integer,DIMENSION(8)                       :: time_miter1,time_miter2, time_multi
      integer,DIMENSION(8)                       :: t_h11,t_h12 

      real(kind=8)                               :: dsum,vsum
      real(kind=8)                               :: tracetot,vectrtot
      real(kind=8),DIMENSION(4)                  :: traces,vectrs
      real(kind=8),DIMENSION(3)                  :: itrace
      real(kind=8),DIMENSION(3,4)                :: vectra,trace
      real(kind=8),DIMENSION(npar,ndat)          :: matinter !to compute aderi*punvar
      real(kind=8),PARAMETER                     :: one = 1.0d0
      real(kind=8),PARAMETER                     :: zero = 0.0d0
!=====================================================================
! Initialization of some arrays:
! VECTRA = sum of the H1 terms for each parameter type (dens., vel...)
!          2: sum of At*Cd**(-1)*DIFF
!          3: sum of At*Cd**(-1)*DIFF + smoothing
!          4: sum of At*Cd**(-1)*DIFF + smoothing + linear relation
! TRACE  = trace of the BDERI array for each parameter type
!          1: trace of At*Cd**(-1)*A
!          2: trace of At*Cd**(-1)*A + Cp**(-1)
!          3: trace of At*Cd**(-1)*A + Cp**(-1) + Dp*Cp**(-1)
!          4: trace of At*Cd**(-1)*A + Cp**(-1) + Ds*Cs**(-1) + Db*Cb**(-1)
! ITRACE = number of diagonal elts of BDERI for each parameter type
! TRACETOT = trace of the final BDERI matrix/npar
! VECTRTOT = total sum of At*Cd**(-1)*DIFF + smoothing + linear relation/npar
!=====================================================================
      if(iiter.eq.1) then
         trace(:,:)=0.d0
         vectra(:,:)=0.d0
         itrace(:)=0.d0

         tracetot=0.d0
         vectrtot=0.d0
         traces(:)=0.d0
         vectrs(:)=0.d0
      end if

      h1(:)=0.d0

      write(inout,*)''
      write(inout,*)'COMPUTING THE MATRIX TO BE INVERTED (BLDMAT)'
      write(inout,*)'    Iteration no ',iiter

      write(*,*)''
      write(*,*)'COMPUTING THE MATRIX TO BE INVERTED (BLDMAT)'
      write(*,*)'    Iteration no ',iiter
      write(*,*)''
!=====================================================================
! Calculation of At*A (of dimension nparXnpar)
! and store the result in the array BDERI 
!=====================================================================
     ! write(inout,*)'    Computing the partial derivative second part'
     ! write(*,*)'    Computing the partial derivative second part'
!=====================================================================
! Time consumming with BLAS subroutine is much better than with MATMUL!!
!=====================================================================
     ! do i=1,npar
      !   do j=1,i
        !    do l=1,ndat
         !      bderi(i,j)=bderi(i,j)+aderi(l,i)*aderi(l,j)*punvar(l)
          !  end do
           ! if(j.ne.i) bderi(j,i)=bderi(i,j)
        ! end do
      ! end do
      CALL DATE_AND_TIME(VALUES=time_miter1)
      write(inout,*)'    Computing the partial derivative first part'
      write(*,*)'    Computing the partial derivative first part'
      do j=1,ndat
         do i=1,npar
            matinter(i,j)=aderi(j,i)*punvar(j)
         end do
      end do
      write(inout,*)'    ....OK'
      write(*,*)'    ....OK'
      CALL DATE_AND_TIME(VALUES=time_miter2)

      CALL TIMECAL(time_miter1,time_miter2)
      write(inout,*)'    Computing the partial derivative second part'
      write(*,*)'    Computing the partial derivative second part'

      CALL DGEMM ('N','N',npar,npar,ndat,one,matinter,npar,aderi,ndat,&
                  zero,bderi, npar)
      CALL DATE_AND_TIME(VALUES=time_multi)
      CALL TIMECAL(time_miter2,time_multi)
      write(inout,*)'    ....OK'
      write(*,*)'    ....OK'
!=====================================================================
! Calculation of At*Cd**(-1)*A. (really does At*A*Cd**(-1))
! Cd**(-1) was previously calculated in array PUNVAR (as 1./sigma**2)
! Initialization of H1 to 0
!=====================================================================

!=====================================================================
! Calculation of H1=At*Cd**(-1)*DIFF for each parameter
!=====================================================================
      write(inout,*)'    Computing the second end member'
      write(*,*)'    Computing the second end member'
      
      !do i=1,npar
       ! do j=1,ndat
        !  h1(i) = h1(i)+aderi(j,i)*punvar(j)*diff(j)
        !end do
     ! end do

      CALL DGEMV('N',npar,ndat,one,matinter,npar,diff,1,zero,h1,1)

      CALL DATE_AND_TIME(VALUES=t_h11)
      CALL TIMECAL(time_multi,t_h11)
      if(iiter.eq.1) then
         do i=1,npar
            j=1
            do while(iend(j).lt.i)
               j=j+1
            end do
            vectra(j,2)=vectra(j,2)+h1(i)
         end do
      end if
      write(inout,*)'    ....OK'
      write(*,*)'    ....OK'

!=====================================================================
! Calculation of At*Cd**(-1)*A + Cp**(-1)
! If (as supposed) all parameters are independent of each other, only
! diagonal elements of Cp**(-1) has to be added. Thus, we just have to
! add VARPAR elements.
!=====================================================================
      write(inout,*)'    Computing the parameter variance part'
      write(*,*)'    Computing the parameter variance part'
      do i=1,npar
	      if(iiter.eq.1) then
	         j=1
	         do while(iend(j).lt.i)
	            j=j+1
	         end do
	         trace(j,1)=trace(j,1)+bderi(i,i)
	         itrace(j)=itrace(j)+1
	      end if

      	bderi(i,i)=bderi(i,i)+varpar(i)

	      if(iiter.eq.1) then
	         j=1
	         do while(iend(j).lt.i)
               j=j+1
            end do
            trace(j,2)=trace(j,2)+bderi(i,i)
	      end if
      end do
      write(inout,*)'    ....OK'
      write(*,*)'    ....OK'
!=====================================================================
! If smoothing constraint is on:
! Calculation of At*Cd**(-1)*A + Cp**(-1) - Ds*Cs**(-1)
! i.e. we add smoothing terms on the diagonal terms of BDERI
!=====================================================================
      if(ismooth.ne.0) then
         write(inout,*)'    Computing the smoothing part'
         write(*,*)'    Computing the smoothing part'
         if(INVD.or.INVGR) then
            CALL DSMOOTH(smooth,iside,jside,xb,yb,par,dsum,&
                         dcont,bderi,h1,ismooth)
         end if
         if(INVV) then
            CALL VSMOOTH(smooth,ivside,jvside,par,ibove,vxnodes,&
                         vynodes,vsum,vcont,bderi,h1,ismooth)
         end if
         write(inout,*)'    ....OK'
         write(*,*)'    ....OK'
      end if

      if(iiter.eq.1) then
         do i=1,npar
            j=1
            do while (iend(j).lt.i)
               j=j+1
            end do
            trace(j,3)=trace(j,3)+bderi(i,i)
            vectra(j,3)=vectra(j,3)+h1(i)
         end do
      end if
!=====================================================================
! Calculation of At*Cd**(-1)*A + Cp**(-1) - Ds*Cs**(-1) - Db*Cb**(-1)
! Adding the matrix connecting densities and velocities
!=====================================================================
      if((INVD.and.INVV).or.(INVGR.and.INVV)) then
         if(ibegin(3).eq.0) then
            write(*,*)''
            write(*,*)'Inconsistency between INVV, INVD (or INVGR) and ibegin(3)'
            write(*,*)'ibegin(3) should not be 0!'
            STOP 'in BLDMAT'
         end if
         write(inout,*)'    Computing the linear relation part'
         write(*,*)'    Computing the linear relation part'
         CALL DENVEL(bderi,h1,ilay,ibove,par,ddvr,nbod) !MP27 aucune correction nécessaire
         write(inout,*)'    ....OK'
         write(*,*)'    ....OK'
      end if

      if(iiter.eq.1) then

         do i=1,npar
            j=1
            do while (iend(j).lt.i)
               j=j+1
            end do
            trace(j,4)=trace(j,4)+bderi(i,i)
            vectra(j,4)=vectra(j,4)+h1(i)
            tracetot=tracetot+bderi(i,i)
            vectrtot=vectrtot+h1(i)
         end do

         tracetot=tracetot/npar
         vectrtot=vectrtot/npar

         do i=1,3
            if(itrace(i).gt.0) then
               do j=4,1,-1
                  if(j.gt.1) then
                     trace(i,j)=trace(i,j)-trace(i,j-1)
                     traces(j)=traces(j)+trace(i,j)
                     trace(i,j)=trace(i,j)/itrace(i)
                     vectra(i,j)=vectra(i,j)-vectra(i,j-1)
                     vectrs(j)=vectrs(j)+vectra(i,j)
                     vectra(i,j)=vectra(i,j)/itrace(i)
                  else
                     traces(1)=traces(1)+trace(i,1)
                     trace(i,1)=trace(i,1)/itrace(i)
                  end if
               end do
            end if
         end do

         traces(:)=traces(:)/npar
         vectrs(:)=vectrs(:)/npar
!=====================================================================
! Writting some outputs
!=====================================================================
         write(inout,*)''
         write(inout,*)'AFTER FIRST ITERATION OF INVERSION:'
         write(inout,*)'    Partial traces: Data, Par, Smooth, B:'
         write(inout,'(4x,4E12.4)') traces(:)
         write(inout,*)'    Total traces, matrices and vectors:'
         write(inout,'(4x,2E12.4)') tracetot,vectrtot

         write(*,*)''
         write(*,*)'AFTER FIRST ITERATION OF INVERSION:'
         write(*,*)'    Partial traces: Data, Par, Smooth, B:'
         write(*,'(4x,4E12.4)') traces(:)
         write(*,*)'    Total traces, matrices and vectors:'
         write(*,'(4x,2E12.4)') tracetot,vectrtot

      end if


      END SUBROUTINE BLDMAT
