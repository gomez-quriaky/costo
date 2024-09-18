!=====================================================================
!=====================================================================
!   SUBROUTINE CALGRA
!=====================================================================
!=====================================================================
!> This subroutine computes in case of INVD=.True. :
!!
!!     - the frechet matrix (ADERI) for the first iteration
!!     - the calculated gravity anomaly CALDATA (in mGal) for
!!       iteration.gt.1 with new parameters PAR
!!     - in case of direct problem, CALDATA (in mGal) with the model
!!       geometry in vel.mod and density PAR
!!     - in case of direct problem, can add random gaussian noise
! 
! Called by MAIN
! calls GBOX,ADNOISE(RAN1)
!=====================================================================
      INCLUDE 'INTERF/MOD_gbox.f'
      INCLUDE 'INTERF/MOD_adnoise.f'

    SUBROUTINE CALGRA(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
                      noised,signoisd)

      USE MOD_unit
      USE MOD_delim
      USE MOD_inv

      USE MOD_gbox
      USE MOD_adnoise
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of CALGRA
!=====================================================================
      integer,intent(in)                         :: iiter,nbod,noised

      real(kind=8),intent(in)                    :: signoisd
      real(kind=8),DIMENSION(:),intent(in)       :: par
      real(kind=8),DIMENSION(:),intent(inout)    :: caldata
      real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
      real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb
      real(kind=8),DIMENSION(:),pointer          :: FX,FY,FZ
!=====================================================================
! Declaration of the dummy arguments of CALGRA
!=====================================================================
      integer                            :: i,j,idp,ndp,ibod,ipar
      integer                            :: idum,ier,ii

      real(kind=8)                       :: rhoinit,sum,moy
      real(kind=8)                       :: xb1,xb2,yb1,yb2,zb1,zb2
      real(kind=8)                       :: g,gFX,gFY,gFZ      
      real(kind=8),PARAMETER             :: gamma = 6.670d-11
      real(kind=8),PARAMETER             :: si2mg = 1.d5
      real(kind=8),PARAMETER             :: km2m  = 1.d3

      real(kind=8)                       :: resz,res2z      
!=====================================================================
! Initialization of caldata to 0. and idum for random noise
!=====================================================================
      caldata(jbegin(1):jend(1)) = 0.d0
      idum=-1

      write(*,*)''
      write(*,*)'FORWARD CALCULATION OF GRAVITY SYNTHETICS'
      write(inout,*)''
      write(inout,*)'FORWARD CALCULATION OF GRAVITY SYNTHETICS'
!=====================================================================
! If this is NOT the first iteration, aderi/=0, and calculate the
! caldata directly (Bouguer anomaly)
!=====================================================================
!      if(iiter.gt.1) then
!
!         do j=jbegin(1),jend(1)
!            do i=ibegin(1),iend(1)
!               caldata(j) = caldata(j) + par(i)*aderi(j,i)
!            end do
!         end do
!=====================================================================
! If this is the first iteration, calculate the gravity effect and
! ADERI array
!=====================================================================
!      else
!=====================================================================
! loop over the data points (idp)
!=====================================================================
         moy = 0.d0
         ndp = 0
         do idp = jbegin(1),jend(1)
            ndp = ndp+1
            gFX = FX(idp)
            gFY = FY(idp)
            gFZ = FZ(idp)
!=====================================================================
! loop over the bodies (ipar,ibod)
!=====================================================================
            do ibod = 1,nbod
               ipar = ibegin(1)+ibod-1
               rhoinit = par(ipar)

               xb1 = xb(ibod,1)
               xb2 = xb(ibod,2)
               yb1 = yb(ibod,1)
               yb2 = yb(ibod,2)
               zb1 = zb(ibod,1)
               zb2 = zb(ibod,2)
   
! on va remplacer GBOX par prism_gx, une des fonctions en C de tesseroid
   call tesseroid_gz(xb1,xb2,yb1,yb2,zb1,zb2,rhoinit,gFX,gFy,-gFZ,resz,res2z) 

               g=resz
               sum=res2z    
       
               caldata(idp) = caldata(idp) + g
!=====================================================================
! In case of direct problem, can add gaussian random noise to the data
!=====================================================================
               if((.not.INV) .and. (noised.eq.1)) then
                  CALL ADNOISE(caldata(idp),signoisd,idum)
               end if
!=====================================================================
! Calculation of the derivatives (ADERI array) in mGal/(g/cm**3)
!=====================================================================
               if(INV) then
                  if(rhoinit.ne.0.d0) then
                     aderi(idp,ipar) = g/rhoinit
                  else
                     aderi(idp,ipar) = sum*gamma*si2mg*km2m*1.d3
                  end if
               end if
!=====================================================================
! End of loop over bodies (ipar)
!=====================================================================
            end do
    
!=====================================================================
! End of loop over data point (idp) and calcul of the mean
!=====================================================================
            moy = moy + caldata(idp)
         end do
         moy = moy/ndp
!=====================================================================
! Remove the mean of the calculated data.
! As the average data input should be zero, the calculated ones should
! be normalized the same way.
!=====================================================================
         do idp=jbegin(1),jend(1)
            caldata(idp) = caldata(idp) - moy
         end do
!=====================================================================
! End of calculation if iiter .gt. 1
!=====================================================================
!      end if

!=====================================================================
! Store the synthetic gravity data in file GRAV.PRED
!=====================================================================
      open(grav1,file='grav.pred',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file GRAV.PRED, logical unit ',&
                    grav1,'.Stooooop in CALGRA!'
         STOP
      end if

      do i=jbegin(1),jend(1)
         write(grav1,'(4f15.3)') FX(i),FY(i),FZ(i),caldata(i)
      end do

      close(grav1,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file GRAV.PRED, logical unit ',&
                    grav1,'.Stooooop in CALGRA!'
         STOP
      end if

      write(*,*)'    Synthetic data stored in GRAV.PRED file'
      write(inout,*)'    Synthetic data stored in GRAV.PRED file'
!=====================================================================
! In case of direct problem, stop the program
!=====================================================================
      if(.not.INV) then
         write(*,*)''
         write(*,*)'End of forward calculation of gravity anomaly'
         write(*,*)'Synthetic data stored in file GRAV.PRED'
         write(*,*)'I''ve finished my job... BYE!'
         write(*,*)''
         write(*,*)'Number of data calculated: ',jend(1)-jbegin(1)+1
         write(*,*)''
         write(inout,*)''
         write(inout,*)'End of forward calculation of gravity anomaly'
         write(inout,*)'Synthetic data stored in file GRAV.PRED'
         write(inout,*)'I''ve finished my job... BYE!'
         write(inout,*)''
         write(inout,*)'Number of data calculated: ',jend(1)-jbegin(1)+1
         write(inout,*)''
         STOP 'THE FORWARD CALCULATION'
      end if

    END SUBROUTINE CALGRA