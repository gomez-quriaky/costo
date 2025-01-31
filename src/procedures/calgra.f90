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
      !INCLUDE 'INTERF/MOD_gbox.f'
      !INCLUDE 'INTERF/MOD_adnoise.f'

    SUBROUTINE CALGRA(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
                      noised,signoisd,val_ad_s,columnAD,rowAD,idx_sp_A)

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

      real(kind=8), DIMENSION(:), intent(inout)  :: val_ad_s
      integer, DIMENSION(:), intent(inout)   :: columnAD
      integer, DIMENSION(:), intent(inout)   :: rowAD
      integer, intent(inout)                 :: idx_sp_A
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

     character*3                         :: scount ! MP new     
     character*20                        :: file_name  ! MP new

 
           
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

               CALL GBOX(gFX,gFY,gFZ,xb1,yb1,zb1,xb2,yb2,zb2,rhoinit,&
                         g,sum)

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
                  idx_sp_A = (idp-1)*nbod + ipar
                  if(rhoinit.ne.0.d0) then
                     aderi(idp,ipar) = g/rhoinit
                     ! sparse Aderi density block 
                     val_ad_s(idx_sp_A) = g/rhoinit
                     columnAD(idx_sp_A) = ipar
                     rowAD(idx_sp_A)   = idp
                  else
                     aderi(idp,ipar) = sum*gamma*si2mg*km2m*1.d3
                     !sparse Aderi density block 
                     val_ad_s(idx_sp_A) = sum*gamma*si2mg*km2m*1.d3
                     columnAD(idx_sp_A) = ipar
                     rowAD(idx_sp_A)   = idp
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

      WRITE(scount,'(i3)') iiter
      file_name  = "grav.pred_iter_"//scount
      open(unit = 444,file = file_name,status ='replace') 

      do i=jbegin(1),jend(1)
         write(grav1,'(4f15.3)') FX(i),FY(i),FZ(i),caldata(i)
         write(444,'(4f15.3)') FX(i),FY(i),FZ(i),caldata(i)
      end do
      close(444)

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

!=====================================================================
!=====================================================================
!   SUBROUTINE GBOX
!=====================================================================
!=====================================================================
!> Computes the vertical attraction of a rectangular prism for one 
!! observation point. 
!!
!! Sides of prism are parallel to x,y,z axes. (BLAKELY)
!!
!! Z AXIS IS VERTICAL DOWN
!!
!! gFX,gFy,gFZ is the data point location (observation point).
!!
!! The prism extends from xb1 to xb2, from yb1 to yb2, and from zb1 to zb2.
!!
!! Density of prism is rhoinit.
!!
!! All distance parameters in units of km;
!!
!! RHO in units of kg/(m**3), G in mGal.
!
! 
! Called by CALGRA
! calls none
!=====================================================================

    SUBROUTINE GBOX(gFX,gFY,gFZ,xb1,yb1,zb1,xb2,yb2,zb2,rhoinit,g,sum)

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of GBOX
!=====================================================================
      real(kind=8),intent(in)              :: rhoinit
      real(kind=8),intent(in)              :: gFX,gFY,gFZ
      real(kind=8),intent(in)              :: xb1,yb1,zb1,xb2,yb2,zb2
      real(kind=8),intent(out)             :: g,sum
!=====================================================================
! Declaration of the dummy arguments of GBOX
!=====================================================================
      integer                              :: i,j,k,ijk
      integer,DIMENSION(2)                 :: isign

      real(kind=8),PARAMETER               :: gamma = 6.670d-11
      real(kind=8),PARAMETER               :: twopi = 6.2831853
      real(kind=8),PARAMETER               :: si2mg = 1.d5
      real(kind=8),PARAMETER               :: km2m  = 1.d3
      real(kind=8)                         :: rho,rijk,arg1,arg2,arg3
      real(kind=8),DIMENSION(2)            :: x,y,z
!=====================================================================
! Initialize some arrays and put the density in KG/M**3
!=====================================================================
      rho=rhoinit*1000.d0
      isign(1)=-1
      isign(2)=1
!=====================================================================
! Compute the distances between fieldpoint and body coordinates
!=====================================================================
      x(1)=gFX-xb1
      y(1)=gFY-yb1
      z(1)=gFZ-zb1
      x(2)=gFX-xb2
      y(2)=gFY-yb2
      z(2)=gFZ-zb2

      sum=0.

      do i=1,2
         do j=1,2
            do k=1,2
               rijk=dsqrt(x(i)**2+y(j)**2+z(k)**2)
               ijk=isign(i)*isign(j)*isign(k)
               arg1=datan2((x(i)*y(j)),(z(k)*rijk))
               if(arg1.lt.0.) then
                  arg1=arg1+twopi
               end if
               arg2=rijk+y(j)
               arg3=rijk+x(i)
               if(arg2.le.0.) then
                  arg2=1.
               end if
               if(arg3.le.0.) then
                  arg3=1.
               end if
               arg2=dlog(arg2)
               arg3=dlog(arg3)
               sum=sum+ijk*(z(k)*arg1-x(i)*arg2-y(j)*arg3)
            end do
         end do
      end do

      g=rho*gamma*sum*si2mg*km2m
!      print*,"dans gbox :: ",rho,gamma,sum,si2mg,km2m,g

    END SUBROUTINE GBOX
