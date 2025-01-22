!=====================================================================
!=====================================================================
!
!			SUBROUTINE PERTURB
!
!> This subroutine calculates the perturbation of the parameters
!! from the inversion of the matrix (BDERI)
!!
!! Then, it checks if the perturbation is not outside the limits
!! given for the gravity in the parameter.inp file
!!
!! If yes, it chooses a factor between 0.1 and 1 that can reduce the
!! perturbation so that the parameter fits in the limits
!!
!! If the factor should be less than 0.1, it means that the perturbation
!! is really too big, and it fixes the parameter to its value without
!! the perturbation.
!
! 04/07/2001: Fixing a bug in line 165... Dealing with par(ibove(ii))
!             without testing if ibove(ii).ne.0
!
! Calls : none
! called by MAIN
!
!=====================================================================
!=====================================================================

    SUBROUTINE PERTURB(bderi,h1,npar,drmax,ilay,par,par0,vels,iiter,&
                       ibove,xb,yb,zb)

      USE MOD_unit
      USE MOD_layer
      USE MOD_inv
      USE MOD_delim

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of PERTURB
!=====================================================================
      integer,intent(in)                            :: npar,iiter
      integer,DIMENSION(:),intent(in)               :: ilay,ibove

      real(kind=8),DIMENSION(:),intent(in)          :: h1,drmax
      real(kind=8),DIMENSION(:),intent(inout)       :: par,par0
      real(kind=8),DIMENSION(:,:),intent(in)        :: bderi,xb,yb,zb
      real(kind=8),DIMENSION(:,:,:,:),intent(inout) :: vels
!=====================================================================
! Declaration of the dummy arguments of PERTURB
!=====================================================================
      integer                         :: i,j,k,ibloc,jfail,nro,nsl
      integer                         :: ivel,ier,ipar,ii,count

      real(kind=8)                    :: factor,newf,romean,rostdv
      real(kind=8)                    :: x,y,z
      real(kind=8)                    :: slmean,slstdv,slmin,slmax,pertu2
      real(kind=8)                    :: romin,romax,slvar,rovar
      real(kind=8),DIMENSION(npar)    :: pertu
!=====================================================================
! Calculation of the perturbation for all the parameters in PERTU
! BDERI is the inverted matrix, H1 is the second hand member of the
! equation.
! new_parameter = previous_parameter + BDERI*H1
!=====================================================================
      pertu(:)=0.d0
      ibloc=0
      factor=1.d0

      write(*,*)''
      write(*,*)'CALCULATION OF THE PARAMETER PERTURBATION'
      write(inout,*)''
      write(inout,*)'CALCULATION OF THE PARAMETER PERTURBATION'

      do i=1,npar
         pertu(i)= dot_product(bderi(i,1:npar),h1(1:npar))
      end do

!=====================================================================
! Writing the perturbation into a dummy file
! for Marie Lopez master seuential inversion with local tomography
! 19-05-2016
!=====================================================================
      open(outmarie,file='marie.out',status='unknown',iostat=ier)
      write(*,*)'ECRITURE DES PERTU'
      if(ier.ne.0) then
         write(*,*)''
         write(*,*)'Error in opening file: MARIE.OUT, logical unit ', &
              outmarie
         STOP 'in MAIN!'
      end if

      do i=1,npar
         write(outmarie,'(2f8.3)') pertu(i), par(i)
      end do
      close(outmarie)
!=====================================================================
! Checking if resulting densities are inside the given limits
! Allowed rho is between +/- DRMAX (read in parameter.inp file)
!=====================================================================
      if(INVD.or.INVGR) then
         count = 0
         write(*,*)'    Checking if resulting densities are inside &
              &the limits'
         write(inout,*)'    Checking if resulting densities are inside &
              &the limits'
!=====================================================================
! Proceeds to linear derivatives
!=====================================================================
         do j=ibegin(1),iend(1)
            ibloc=ibloc+1
            pertu2=DABS(par(j)+pertu(j))
            if(pertu2.gt.drmax(ilay(ibloc))) then
               newf=(drmax(ilay(ibloc))-DABS(par(j)))/DABS(pertu(j))
!=====================================================================
! If the factor by which the parameter increments have to be
! multiplied becomes less than 0.1 (10%), the corresponding parameter
! is fixed for this iteration (ie perturbation set to 0.).
! Because the reduction of value becomes too big.
!=====================================================================
               if(newf.lt.(1.d-1)) then
                  count=count+1
                  newf=1.d0
                  pertu(j)=0.d0
                  x=(xb(ibloc,1)+xb(ibloc,2))/2
                  y=(yb(ibloc,1)+yb(ibloc,2))/2
                  z=(zb(ibloc,1)+zb(ibloc,2))/2
                  write(inout,*)''
                  write(inout,'(5x,''Density parameter no.'',i5,''('',f6.1,&
                       &'','',f6.1,'','',f6.1,'') fixed'')') j,x,y,z
               end if
               if(newf.lt.factor) then
                  factor=newf
		            pertu(j)=pertu(j)*factor
                  jfail=j
               end if
            end if
         end do
         if(count.gt.0) then
            write(*,*)'    WARNING: some density parameters had to be fixed'
         end if
      end if
      
      if(factor.lt.1.d0) then
         write(inout,*)''
         write(inout,'(''    Differential parameter vector reduced &
              &by'',f10.4)') factor
         write(inout,'(''    Parameter no. '',i5,'' touches limits'')') jfail

         write(*,*)''
         write(*,'(''    Differential parameter vector reduced by'',&
              &f10.4)') factor
         write(*,'(''    Parameter no. '',i5,'' touches limits'')') jfail
      end if
!=====================================================================
! Calculation of the new parameter in PAR array, and storage of the
! old one in PAR0
! BUG: the line par(i)=par(i)+pertu(i)*factor is wrong here... because
!      factor will have the value of the last pass in the if loop...
!      It has to be 1 for the velocity
!      I suppressed the line, and added pertu(j)=pertu(j)*factor above
!=====================================================================
      do i=1,npar
         par0(i)=par(i)
	      par(i)=par(i)+pertu(i)
!         par(i)=par(i)+pertu(i)*factor
      end do
!=====================================================================
! Calculation of the standard deviation of parameters
!=====================================================================
      romean = 0.d0
      rostdv = 0.d0
      slmean = 0.d0
      slstdv = 0.d0
      nro = 0
      nsl = 0
      ivel = 0
      ii = 0
!=====================================================================
! First deal with the velocity
!=====================================================================

      if(INVV) then
         slmin = par(ibegin(2))
         slmax = par(ibegin(2))
         do i=ibegin(2),iend(2)
            ivel=ivel+1
            if(ibove(ivel).ne.0) then
               nsl = nsl +1
               slmean = slmean + par(ibove(ivel))
               if(slmin.gt.par(ibove(ivel))) slmin = par(ibove(ivel))
               if(slmax.lt.par(ibove(ivel))) slmax = par(ibove(ivel))
            end if
         end do

         if(nsl.gt.1) then
            slmean = slmean/nsl
            do i=ibegin(2),iend(2)
               ii=ii+1
               if(ibove(ii).ne.0) then
                  slstdv = slstdv + (par(ibove(ii))-slmean)**2
               end if
            end do
            slvar = slstdv/(nsl-1)
            slstdv = DSQRT(slvar)
            write(*,*)'    Writting velocity statistics in PARAMETER.OUT'
            write(*,*)'    Velocities:  Min.    Max.    Mean.    St.Dev.'
            write(*,'(10x,3(f10.3),f10.2)') slmin,slmax,slmean,slstdv
            write(inout,*)''
            write(inout,*)'    Velocities:  Min.    Max.    Mean.    St.Dev.'
            write(inout,'(10x,3(f10.3),f10.2)') slmin,slmax,slmean,slstdv
         end if
!=====================================================================
! Have to store the new perturbed velocity in VELS(i,j,k,iiter+1) array
! so that the raytracing model will be also modified for the next
! iteration
!=====================================================================
         write(*,*)''
         write(*,*)'    Storing the new velocity model in VELS array'
         write(inout,*)''
         write(inout,*)'    Storing the new velocity model in VELS array'
         ipar=ibegin(2)
         do k=1,nznode-1
            do j=1,nynode
               do i=1,nxnode
                  vels(i,j,k,iiter+1)=vels(i,j,k,1) + par(ipar)
!                  vels(i,j,k,iiter+1)=vels(i,j,k,iiter) + pertu(ipar)
                  ipar=ipar+1
               end do
            end do
         end do
!=====================================================================
! The last layer of the velocity model is just here for the 3D raytracing
! limits.
! So, we just have to set it to the previous value without perturbation
!=====================================================================
         vels(:,:,nznode,iiter+1)=vels(:,:,nznode,iiter)

      end if

!=====================================================================
! Second deal with the density
!=====================================================================
      if(INVD.or.INVGR) then
         romin = par(ibegin(1))
         romax = par(ibegin(1))
         do i=ibegin(1),iend(1)
            nro = nro + 1
            romean = romean + par(i)
            if(romin.gt.par(i)) romin = par(i)
            if(romax.lt.par(i)) romax = par(i)
         end do

         if(nro.gt.1) then
            romean = romean/nro
            do i=ibegin(1),iend(1)
               rostdv = rostdv + (par(i)-romean)**2
            end do
            rovar = rostdv/(nro-1)
            rostdv = DSQRT(rovar)
            write(*,*)'    Writting density statistics in PARAMETER.OUT'
            write(*,*)'    Densities:  Min.    Max.    Mean.    St.Dev.'
            write(*,'(10x,3(f10.3),f10.2)') romin,romax,romean,rostdv
            write(inout,*)''
            write(inout,*)'    Densities:  Min.    Max.    Mean.    St.Dev.'
            write(inout,'(10x,3(f10.3),f10.2)') romin,romax,romean,rostdv
         end if
      end if

!=====================================================================
! Finally deal with the density-velocity coeff
!=====================================================================
      if((INVD.and.INVV).or.(INVGR.and.INVV)) then
         write(*,*)''
         write(*,*)'RESULTS OF PARAMETERS FOR B-COEFF'
         write(*,*)'  PAR(i)      PAR0(i)  PERTURBATION(i)'
         do i=ibegin(3),iend(3)
            write(*,'(3(f10.4,1x))') par(i),par0(i),pertu(i)
         end do
         write(*,*)''

         write(inout,*)''
         write(inout,*)'RESULTS OF PARAMETERS FOR B-COEFF'
         write(inout,*)'  PAR(i)      PAR0(i)  PERTURBATION(i)'
         do i=ibegin(3),iend(3)
            write(inout,'(3(f10.4,1x))') par(i),par0(i),pertu(i)
         end do
         write(inout,*)''
      end if

    END SUBROUTINE PERTURB




