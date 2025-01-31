!=====================================================================
!=====================================================================
!
!			SUBROUTINE OUTPUT
!
!> This subroutine writes the outputs :
!!         - velocity perturbation in percent
!!         - density perturbation
!!         - statistics
!
! Calls : none
! called by MAIN
!
!=====================================================================
!=====================================================================

    SUBROUTINE OUTPUT(iter,ddtot,grdtot,vdtot,bdtot,velvar,rovar,par,ibove,xb,yb,&
                     zb,vxnodes,vynodes,vznodes,vels,varv,varg,vargr,epsi,epsil,&
                     rmsg,rmst,rmsgr,errclt,errclgr,errclg,rtvar,ncomp)

      USE MOD_unit
      USE MOD_layer
      USE MOD_inv
      USE MOD_delim

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of OUTPUT
!=====================================================================
      integer,intent(in)                         :: iter
      integer,DIMENSION(:),intent(in)            :: ibove

      real(kind=8),intent(in)                    :: varv,varg,epsi,epsil
      real(kind=8),DIMENSION(6),intent(in)       :: vargr      
      real(kind=8),DIMENSION(:),intent(in)       :: par
      real(kind=8),DIMENSION(:),intent(in)       :: rmsg,rmst,errclt,errclg
      real(kind=8),DIMENSION(:),intent(in)       :: vxnodes,vynodes,vznodes
      real(kind=8),DIMENSION(:,:),intent(in)     :: ddtot,vdtot,bdtot
      real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb,rovar,velvar
      real(kind=8),DIMENSION(:,:,:,:),intent(in) :: vels
      real(kind=8),DIMENSION(:),intent(inout)    :: grdtot
      real(kind=8),DIMENSION(:,:),intent(inout)  :: rmsgr,errclgr

      real(kind=8),DIMENSION(:),intent(in)    :: rtvar         ! MP36
      integer,intent(inout)                   :: ncomp         ! MP36
      real(kind=8)      :: var1,var2,var3,var4,var5,var6
!=====================================================================
! Declaration of the dummy arguments of OUTPUT
!=====================================================================
      character(len=80)                 :: form1

      integer                           :: i,j,k,l,ier,ii,ipar,iterr,nbn

      real(kind=8)                      :: sum1,sum2,sum3,sum4,dvel
      real(kind=8)                      :: x,y,z,deltag,deltat,deltagr
   
!=====================================================================
! First output the mean differences.
!=====================================================================
      nbn=0     ! MP37 pour prendre en compte le nombre de méthode que l'on inverse en plus de la gradio
      if (INVD) then
        nbn=nbn+1
      endif
      if (INVV) then
        nbn=nbn+1
      endif

      write(*,*)''
      write(*,*)'********************************************************'
      write(*,*)'           FINAL RESULT OUTPUTS (iter=',iter,')'
      write(*,*)'********************************************************'
      write(*,*)''
      write(*,*)''
      write(*,*)'MEAN DIFFERENCES'
      write(*,'(''    Total should have been less than '',f7.5,'' to stop&
           &, or'')') epsi
      write(*,'(''    Difference between two iterations should have been&
           & less than '',f7.5,'' to stop'')') epsil
      write(*,*)'    ITER     Total        Data       &
           & Parameter       Smooth'
      
      write(inout,*)''
      write(inout,*)'********************************************************'
      write(inout,*)'           FINAL RESULT OUTPUTS (iter=',iter,')'
      write(inout,*)'********************************************************'
      write(inout,*)''
      write(inout,*)''
      write(inout,*)'MEAN DIFFERENCES'
      write(inout,'(''    Total should have been less than '',f5.3,'' to stop&
           &, or'')') epsi
      write(inout,'(''    Difference between two iterations should have been&
           & less than '',f5.3,'' to stop'')') epsil
      write(inout,*)'    ITER     Total        Data       &
           & Parameter       Smooth'

!  rajout du terme gradio. DDTOT est composé de 1) mean diff between 
!     obs and calc gravity data 2) mean diff between param 3) roughness 
!     from smoothing constraint 4)somme des 3

      do i=1,iter+1
         sum1 = ddtot(i,1) + vdtot(i,1) + grdtot(i)
         sum2 = ddtot(i,2) + vdtot(i,2) + bdtot(i,2)
         sum3 = ddtot(i,3) + vdtot(i,3)
         sum4 = ddtot(i,4) + vdtot(i,4) + bdtot(i,3) + grdtot(i)
         write(inout,'(5x,i3,1x,4(f12.3,1x))') i,sum4,sum1,sum2,sum3
         write(*,'(5x,i3,1x,4(f12.3,1x))') i,sum4,sum1,sum2,sum3
      end do
!=====================================================================
! Second output the Variances through the iterations
!=====================================================================
      form1='(5x,i3,1x,f12.3,''/'',f7.3,1x,f12.3,''/'',f7.3,1x,f12.3,''/'',f5.3)'
      write(*,*)''
      write(*,*)'VARIANCE OF CALCULATED DATA THROUGH ITERATIONS'
      write(*,*)'    ITER          GRAVITY           DELAYTIME           FTG (mean all comp.)'
      write(*,*)'                Var/Std.Dev.       Var/Std.Dev.              Var/Std.Dev.'

      write(inout,*)''
      write(inout,*)'VARIANCE OF CALCULATED DATA THROUGH ITERATIONS'
      write(inout,*)'    ITER          GRAVITY           DELAYTIME           FTG (mean all comp.)'
      write(inout,*)'                Var/Std.Dev.       Var/Std.Dev.              Var/Std.Dev.'

      do i=1,iter+1
         write(inout,form1) i,rovar(i,1),DSQRT(rovar(i,1)),velvar(i,1),&
              DSQRT(velvar(i,1)),rovar(i,3),DSQRT(rovar(i,3))
         write(*,form1) i,rovar(i,1),DSQRT(rovar(i,1)),velvar(i,1),&
              DSQRT(velvar(i,1)),rovar(i,3),DSQRT(rovar(i,3))
      end do

      write(*,*)''
      write(*,*)'VARIANCE AND STANDARD DEVIATION OF INITIAL DATA'
      write(*,*)'                  GRAVITY           DELAYTIME           FTG  (mean all comp.)'
      write(*,*)'                Var/Std.Dev.       Var/Std.Dev.              Var/Std.Dev.'
!      write(*,'(9x,f12.3,''/'',f7.3,1x,f12.3,''/'',f7.3,1x,f12.3,''/'',f5.3)') varg,DSQRT(varg),&
!           varv,DSQRT(varv),sum(vargr(:))/6,DSQRT(sum(vargr(:))/6)
      write(*,'(9x,f12.3,''/'',f7.3,1x,f12.3,''/'',f7.3,1x,f12.3,''/'',f5.3)') varg,DSQRT(varg),&           !MP36
           varv,DSQRT(varv),sum(vargr(:))/ncomp,DSQRT(sum(vargr(:))/ncomp)          

      write(inout,*)''
      write(inout,*)'VARIANCE AND STANDARD DEVIATION OF INITIAL DATA'
      write(inout,*)'                  GRAVITY           DELAYTIME           FTG  (mean all comp.)'
      write(inout,*)'                Var/Std.Dev.       Var/Std.Dev.              Var/Std.Dev.'
!      write(inout,'(9x,f12.3,''/'',f7.3,1x,f12.3,''/'',f5.3)') varg,&
!           DSQRT(varg),varv,DSQRT(varv),sum(vargr(:))/6,DSQRT(sum(vargr(:))/6)
      write(inout,'(9x,f12.3,''/'',f7.3,1x,f12.3,''/'',f5.3)') varg,&                  !MP36
           DSQRT(varg),varv,DSQRT(varv),sum(vargr(:))/ncomp,DSQRT(sum(vargr(:))/ncomp)

            iterr=0
            if (rtvar(nbn+1) .ne. 0) then
               iterr=iterr+1
               var1=vargr(iterr)
            else
               var1=0   
            endif
            if (rtvar(nbn+2) .ne. 0) then
               iterr=iterr+1
               var2=vargr(iterr)
            else
               var2=0   
            endif
            if (rtvar(nbn+3) .ne. 0) then
               iterr=iterr+1
               var3=vargr(iterr)
            else
               var3=0   
            endif
            if (rtvar(nbn+4) .ne. 0) then
               iterr=iterr+1
               var4=vargr(iterr)
            else
               var4=0   
            endif
            if (rtvar(nbn+5) .ne. 0) then
               iterr=iterr+1
               var5=vargr(iterr)
            else
               var5=0   
            endif
            if (rtvar(nbn+6) .ne. 0) then
               iterr=iterr+1
               var6=vargr(iterr)
            else
               var6=0   
            endif

      write(*,*)''
      write(*,*)'DETAILS FOR ACH COMPONENTS OF FTG'
      write(*,*)'          XX          XY          XZ           YY            YZ           ZZ'
      write(*,'(1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3)') &
         var1,var2,var3,var4,var5,var6             ! MP36

      write(inout,*)''
      write(inout,*)'DETAILS FOR ACH COMPONENTS OF FTG'
      write(inout,*)'          XX          XY          XZ           YY            YZ           ZZ'
      write(inout,'(1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3)') &
         var1,var2,var3,var4,var5,var6            ! MP36

      write(*,*)''
      write(*,*)'VARIANCES OF PARAMETERS THROUGH ITERATIONS'
      write(*,*)'    ITER     Var.dens.     Var.vel.'

      write(inout,*)''
      write(inout,*)'VARIANCES OF PARAMETERS THROUGH ITERATIONS'
      write(inout,*)'    ITER     Var.dens.     Var.vel.'
      do i=1,iter+1
         write(inout,'(5x,i3,1x,f12.3,1x,f12.3)') i,rovar(i,2),velvar(i,2)
         write(*,'(5x,i3,1x,f12.3,1x,f12.3)') i,rovar(i,2),velvar(i,2)
      end do

!=====================================================================
! Output the RMS of the obs-calc data
!=====================================================================
      write(*,*)''
      write(*,*)'RMS OF THE (OBS-CALC) DATA (HAS TO DECREASE)'
      write(*,*)'    ITER     RMS.GRAV      RMS.DT      RMS.GRADIO'

      write(inout,*)''
      write(inout,*)'RMS OF THE (OBS-CALC) DATA (HAS TO DECREASE)'
      write(inout,*)'    ITER     RMS.GRAV      RMS.DT      RMS.GRADIO'

      do i=1,iter+1
         write(*,'(5x,i3,1x,f12.3,1x,f12.3,1x,f12.3)') i,rmsg(i),rmst(i),rmsgr(i,1)
         write(inout,'(5x,i3,1x,f12.3,1x,f12.3,1x,f12.3)') i,rmsg(i),rmst(i),rmsgr(i,1)
      end do

      write(*,*)''
      write(*,*)'DETAILS FOR ACH COMPONENTS OF FTG (HAS TO DECREASE TOO)'
      write(*,*)'    ITER          XX          XY          XZ           YY            YZ           ZZ'

      write(inout,*)''
      write(inout,*)'DETAILS FOR ACH COMPONENTS OF FTG (HAS TO DECREASE TOO)'
      write(inout,*)'    ITER          XX          XY          XZ           YY            YZ           ZZ'

      do i=1,iter+1
         write(*,'(5x,i3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3)') i,&
         rmsgr(i,2),rmsgr(i,3),rmsgr(i,4),rmsgr(i,5),rmsgr(i,6),rmsgr(i,7)
         write(inout,'(5x,i3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3)') i,&
         rmsgr(i,2),rmsgr(i,3),rmsgr(i,4),rmsgr(i,5),rmsgr(i,6),rmsgr(i,7)         
      end do

      open(unit=333,file='RMS_FTG',status='REPLACE')      
      do i=1,iter+1
      write(333,'(5x,i3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3)') i,&
         rmsgr(i,1),rmsgr(i,2),rmsgr(i,3),rmsgr(i,4),rmsgr(i,5),rmsgr(i,6),rmsgr(i,7)
      end do
      close(333)

      open(unit=333,file='RMS_TOMO',status='REPLACE')      
      do i=1,iter+1
      write(333,'(5x,i3,1x,f12.3)') i,rmst(i)
      end do
      close(333)

      open(unit=333,file='RMS_GRAVI',status='REPLACE')      
      do i=1,iter+1
      write(333,'(5x,i3,1x,f12.3)') i,rmsg(i)
      end do
      close(333)

      if(rmsg(1).ne.0) then
         deltag = (rmsg(1)-rmsg(iter+1))*100/rmsg(1)
      end if
      if(rmst(1).ne.0) then
         deltat = (rmst(1)-rmst(iter+1))*100/rmst(1)
      end if

      if(rmsgr(1,1).ne.0) then
         deltagr = (rmsgr(1,1)-rmsgr(iter+1,1))*100/rmsgr(1,1)
      end if

      write(*,*)''
      if(INVD) then
         write(*,'(''TOTAL DECREASE OF THE GRAVITY DATA RMS:'',f7.2,''%'')') &
              deltag
      end if
      if(INVV) then
         write(*,'(''TOTAL DECREASE OF THE  DELAY TIME  RMS:'',f7.2,''%'')') &
              deltat
      end if
      
      if(INVGR) then
         write(*,'(''TOTAL DECREASE OF THE GRADIO DATA RMS:'',f7.2,''%'')') &
              deltagr
      end if

      write(inout,*)''
      if(INVD) then
         write(inout,'(''TOTAL DECREASE OF THE GRAVITY DATA RMS'',f7.2,''%'')') &
              deltag
      end if
      if(INVV) then
         write(inout,'(''TOTAL DECREASE OF THE  DELAY TIME  RMS'',f7.2,''%'')') &
              deltat
      end if
      
      if(INVGR) then
         write(inout,'(''TOTAL DECREASE OF THE GRADIO DATA RMS:'',f7.2,''%'')') &
              deltagr
      end if

!=====================================================================
! Open the file VELOCITY.RES and store the velocity contrasts in it
!=====================================================================
      if(INVV) then

         open(velout,file='velocity.res',status='unknown',iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in opening the file VELOCITY.RES, logical unit ',&
                 velout
            STOP 'in OUTPUT!'
         end if
         write(*,*)''
         write(*,*)'STORING THE PERTURBATIONS IN VELOCITY.RES FILE'
         write(*,*)'    format: x(km) y(km) z(km) par vel(%) vels'
         write(inout,*)''
         write(inout,*)'STORING THE PERTURBATIONS IN VELOCITY.RES FILE'
         write(inout,*)'    format: x(km) y(km) z(km) par vel(%) vels'

         ipar = 0
         do k=1,nznode-1
            do j=1,nynode
               do i=1,nxnode
                  ipar = ipar + 1
                  if(ibove(ipar).ne.0) then
                     dvel=par(ibove(ipar))/vels(i,j,k,1)*100
                  else
                     dvel=0.d0
                  end if
                  write(velout,'(4(f10.4,2x),10f10.3)') vxnodes(i),vynodes(j),&
                       vznodes(k),dvel,(vels(i,j,k,l),l=1,iter+1)
               end do
            end do
         end do

         close(velout,iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in closing the file VELOCITY.RES, logical unit ',&
                 velout
            STOP 'in OUTPUT!'
         end if

      end if

!=====================================================================
! Open the file DENSITY.RES and store the density variations in it
!=====================================================================
      if(INVD.or.INVGR) then

         open(denout,file='density.res',status='unknown',iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in opening the file DENSITY.RES, logical unit ',&
                 denout
            STOP 'in OUTPUT!'
         end if
      
         write(*,*)''
         write(*,*)'STORING THE VARIATIONS IN DENSITY.RES FILE'
         write(*,*)'    format: x(mid.km) y(mid.km) z(mid.km) Drho(g/cm3)'
         write(inout,*)''
         write(inout,*)'STORING THE VARIATIONS IN DENSITY.RES FILE'
         write(inout,*)'    format: x(mid.km) y(mid.km) z(mid.km) Drho(g/cm3)'
         
         ii = 0
         do i=ibegin(1),iend(1)
            ii = ii + 1
            x=(xb(ii,1)+xb(ii,2))*0.5
            y=(yb(ii,1)+yb(ii,2))*0.5
            z=(zb(ii,1)+zb(ii,2))*0.5
            write(denout,'(4(f10.3,2x))') x,y,z,par(i)
         end do
         
         close(denout,iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in closing the file DENSITY.RES, logical unit ',&
                 denout
            STOP 'in OUTPUT!'
         end if

      end if
    
!=====================================================================
! Finally put the B-value at screen and in PARAMETER.OUT file
!=====================================================================

      if((INV.and.INVV.and.INVD).or.(INV.and.INVV.and.INVGR)) then

         write(*,*)''
         write(*,*)'STORING THE B-RESULTS IN PARAMETER.OUT FILE'
         write(*,*)'    layer B(km.cm3/s.g)'
         write(inout,*)''
         write(inout,*)'STORING THE B-RESULTS IN PARAMETER.OUT FILE'
         write(inout,*)'    layer B(km.cm3/s.g)'

         ipar=ibegin(3)-1
         do k=1,nlayer
            ipar=ipar+1
            write(*,'(i4,3x,f10.3)') k,par(ipar)
            write(inout,'(i4,3x,f10.3)') k,par(ipar)
         end do

      end if

      write(*,*)''
      write(*,*)'ERROR FROM BASUYAU'
      write(*,*)'    ITER     ERRCLT     ERRCLGR. MEAN'

      write(inout,*)''
      write(inout,*)'ERROR FROM BASUYAU'
      write(inout,*)'    ITER     ERRCLT     ERRCLGR. MEAN'
           
      do i=1,iter+1
         write(inout,'(5x,i3,1x,f12.3,1x,f12.3)') i,errclt(i),errclgr(i,1)
         write(*,'(5x,i3,1x,f12.3,1x,f12.3)') i,errclt(i),errclgr(i,1)
      end do

    END SUBROUTINE OUTPUT

