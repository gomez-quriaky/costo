!=====================================================================
!   SUBROUTINE CALGRADIO
!=====================================================================

subroutine CALGRADIO(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
	  							noisegr,signoisgr,rtvar,ncomp)


USE MOD_unit
USE MOD_delim
USE MOD_inv

USE MOD_adnoise

IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of CALGRA
!=====================================================================
      integer,intent(in)                         :: iiter,nbod,noisegr
      real(kind=8),intent(in)                    :: signoisgr      
      real(kind=8),DIMENSION(:),intent(in)       :: par      
      real(kind=8),DIMENSION(:),intent(inout)    :: caldata
      real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
      real(kind=8),DIMENSION(:),pointer          :: FX,FY,FZ
      real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb,zb


!=====================================================================
! Declaration of the dummy arguments of CALGRA
!=====================================================================
      integer                            :: i,j,idp,ndp,ibod,ipar,iter
      integer                            :: idum,ier,nptsgr,nbn
      real(kind=8)                       :: rhoinit,sum
      real(kind=8)                       :: moyxx,moyxy,moyxz,moyyy,moyyz,moyzz
      real(kind=8)                       :: xb1,xb2,yb1,yb2,zb1,zb2
      real(kind=8)                       :: gFX,gFY,gFZ     
      real(kind=8)                       :: resxx,resxy,resxz,resyy,resyz,reszz
      real(kind=8)                       :: res2xx,res2xy,res2xz,res2yy,res2yz,res2zz

      character*3                         :: scount 
      character*20                        :: file_name  


      real(kind=8),DIMENSION(:),intent(in)    :: rtvar         
      integer,intent(inout)                   :: ncomp         
      real(kind=8)   :: XX,XY,XZ,YY,YZ,ZZ   

!=====================================================================
! Initialization of caldata to 0. and idum for random noise
!=====================================================================
      caldata(jbegin(3):jend(3)) = 0.d0
      idum=-1

      write(*,*)''
      write(*,*)'FORWARD CALCULATION OF SYNTHETIC FTG'
      write(inout,*)''
      write(inout,*)'FORWARD CALCULATION OF SYNTHETIC FTG'

!=====================================================================
! loop over the data points (idp)
!=====================================================================
         moyxx = 0.d0
         moyxy = 0.d0
         moyxz = 0.d0
         moyyy = 0.d0
         moyyz = 0.d0
         moyzz = 0.d0
         ndp = 0

! la en gravi on fait une boucle sur toutes les mesures. MAIS en gradio, les données
!  sont rangées de telle façon : tout d'abord tous les XX, ensuite tous les
!  XY, ensuite les XZ etc... jusq'à ZZ. Donc pas besoin de faire de boucle sur toutes
!  les données, uniquement sur le vrai nombre de point (donc la 1ere composante XX). Avec
!  ces coordonnées on pourra calculer l effet à chaque point et ensuite il faudra juste bien
!  reranger les réponses

      nbn=0    
      if (INVD) then
        nbn=nbn+1
      endif
      if (INVV) then
        nbn=nbn+1
      endif

         nptsgr=(jend(3)-jbegin(3)+1)/ncomp

         do idp = jbegin(3),jbegin(3)+nptsgr-1

            ndp = ndp+1
            gFX = FX(idp)
            gFY = FY(idp)
            gFZ = FZ(idp)

!=====================================================================
! loop over the bodies (ipar,ibod)
!=====================================================================
!            do ibod = 1,nbod
            do ibod = 1,nbod               
               ipar = ibegin(1)+ibod-1
               rhoinit = par(ipar)

               xb1 = xb(ibod,1)
               xb2 = xb(ibod,2)
               yb1 = yb(ibod,1)
               yb2 = yb(ibod,2)
               zb1 = zb(ibod,1)
               zb2 = zb(ibod,2)

! ICI ON APPELLE LE CODE DE CALCUL DIRECT DE LA GRADIO
! on peut laisser ça comme tel, il calculera l'ensemble des réponses du tenseur mais ensuite on ne garde que celle
!           qui nous intéresse. Comme ça pas besoin de modifier la subroutine tesseroid écrite en C
               call tesseroid(xb1,xb2,yb1,yb2,zb1,zb2,rhoinit,gFX,gFy,-gFZ,&
                     resxx,resxy,resxz,resyy,resyz,reszz,&
                     res2xx,res2xy,res2xz,res2yy,res2yz,res2zz)              

! Faudra faire gaffe ici avec le problème des 6 composantes
! en gravi, on récupère g à la sortie de gbox et pour chaque point, on va avoir un
!  g // en gradio c'est différent, car les données sont rangées de telle façon : 
!  tout d'abord tous les XX, ensuite tous les XY, ensuite les XZ etc... jusq'à ZZ

         iter=0
         if (rtvar(nbn+1) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + resxx
         endif
         if (rtvar(nbn+2) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + resxy
         endif
         if (rtvar(nbn+3) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + resxz
         endif
         if (rtvar(nbn+4) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + resyy
         endif
         if (rtvar(nbn+5) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + resyz
         endif
         if (rtvar(nbn+6) .ne. 0) then
            iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) + reszz
         endif

!=====================================================================
! In case of direct problem, can add gaussian random noise to the data
!=====================================================================
               if((.not.INV) .and. (noisegr.eq.1)) then
                  CALL ADNOISE(caldata(idp+nptsgr*0),signoisgr,idum)
                  CALL ADNOISE(caldata(idp+nptsgr*1),signoisgr,idum)
                  CALL ADNOISE(caldata(idp+nptsgr*2),signoisgr,idum)
                  CALL ADNOISE(caldata(idp+nptsgr*3),signoisgr,idum)
                  CALL ADNOISE(caldata(idp+nptsgr*4),signoisgr,idum)
                  CALL ADNOISE(caldata(idp+nptsgr*5),signoisgr,idum)
               end if

!=====================================================================
! Calculation of the derivatives (ADERI array) in mGal/(g/cm**3)
!=====================================================================
               if(INV) then
                  if(rhoinit.ne.0.d0) then
                     iter=0
                     if (rtvar(nbn+1) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=resxx/rhoinit
                     endif
                     if (rtvar(nbn+2) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=resxy/rhoinit
                     endif
                     if (rtvar(nbn+3) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=resxz/rhoinit
                     endif
                     if (rtvar(nbn+4) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=resyy/rhoinit
                     endif
                     if (rtvar(nbn+5) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=resyz/rhoinit
                     endif
                     if (rtvar(nbn+6) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=reszz/rhoinit
                     endif

                  else

                     iter=0
                     if (rtvar(nbn+1) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2xx*1.d3
                     endif
                     if (rtvar(nbn+2) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2xy*1.d3
                     endif
                     if (rtvar(nbn+3) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2xz*1.d3
                     endif
                     if (rtvar(nbn+4) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2yy*1.d3
                     endif
                     if (rtvar(nbn+5) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2yz*1.d3
                     endif
                     if (rtvar(nbn+6) .ne. 0) then
                        iter=iter+1
                     aderi(idp+nptsgr*(iter-1),ipar)=res2zz*1.d3
                     endif

                  end if
               end if

!=====================================================================
! End of loop over bodies (ipar)
!=====================================================================
            end do
!=====================================================================
! End of loop over data point (idp) and calcul of the mean
!=====================================================================
            iter=0
            if (rtvar(nbn+1) .ne. 0) then
               iter=iter+1
            moyxx = moyxx + caldata(idp+nptsgr*(iter-1))
            else
            moyxx=0.0
            endif

            if (rtvar(nbn+2) .ne. 0) then
               iter=iter+1
            moyxy = moyxy + caldata(idp+nptsgr*(iter-1))
            else
            moyxy=0.0
            endif

            if (rtvar(nbn+3) .ne. 0) then
               iter=iter+1
            moyxz = moyxz + caldata(idp+nptsgr*(iter-1))
            else
            moyxz=0.0
            endif

            if (rtvar(nbn+4) .ne. 0) then
               iter=iter+1
            moyyy = moyyy + caldata(idp+nptsgr*(iter-1))
            else
            moyyy=0.0
            endif

            if (rtvar(nbn+5) .ne. 0) then
               iter=iter+1
            moyyz = moyyz + caldata(idp+nptsgr*(iter-1))
            else
            moyyz=0.0
            endif

            if (rtvar(nbn+6) .ne. 0) then
               iter=iter+1
            moyzz = moyzz + caldata(idp+nptsgr*(iter-1))
            else
            moyzz=0.0
            endif

!=====================================================================
! End of loop over data point (idp) and calcul of the mean
!=====================================================================
         end do

moyxx = moyxx/nptsgr
moyxy = moyxy/nptsgr
moyxz = moyxz/nptsgr
moyyy = moyyy/nptsgr
moyyz = moyyz/nptsgr
moyzz = moyzz/nptsgr

!=====================================================================
! Remove the mean of the calculated data.
! As the average data input should be zero, the calculated ones should
! be normalized the same way.
!=====================================================================
do idp=jbegin(3),jbegin(3)+nptsgr-1

            iter=0
            if (rtvar(nbn+1) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyxx
            endif

            if (rtvar(nbn+2) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyxy
            endif

            if (rtvar(nbn+3) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyxz
            endif

            if (rtvar(nbn+4) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyyy
            endif

            if (rtvar(nbn+5) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyyz
            endif

            if (rtvar(nbn+6) .ne. 0) then
               iter=iter+1
               caldata(idp+nptsgr*(iter-1)) = caldata(idp+nptsgr*(iter-1)) - moyzz
            endif

end do

!=====================================================================
! Store the synthetic FTG data in file FTG.PRED
!=====================================================================
      open(ftg1,file='FTG.pred',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file FTG.PRED, logical unit ',&
                    ftg1,'.Stooooop in CALGRADIO!'
         STOP
      end if

      WRITE(scount,'(i3)') iiter
      file_name  = scount//"_iter_FTG.pred"
      open(unit = 444,file = file_name,status ='replace') 

      do idp=jbegin(3),jbegin(3)+nptsgr-1

            iter=0
            if (rtvar(nbn+1) .ne. 0) then
               iter=iter+1
               XX=caldata(idp+nptsgr*(iter-1))
            else
               XX=0   
            endif

            if (rtvar(nbn+2) .ne. 0) then
               iter=iter+1
               XY=caldata(idp+nptsgr*(iter-1))
            else
               XY=0   
            endif

            if (rtvar(nbn+3) .ne. 0) then
               iter=iter+1
               XZ=caldata(idp+nptsgr*(iter-1))
            else
               XZ=0   
            endif

            if (rtvar(nbn+4) .ne. 0) then
               iter=iter+1
               YY=caldata(idp+nptsgr*(iter-1))
            else
               YY=0   
            endif

            if (rtvar(nbn+5) .ne. 0) then
               iter=iter+1
               YZ=caldata(idp+nptsgr*(iter-1))
            else
               YZ=0   
            endif

            if (rtvar(nbn+6) .ne. 0) then
               iter=iter+1
               ZZ=caldata(idp+nptsgr*(iter-1))
            else
               ZZ=0   
            endif

         write(ftg1,'(1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9)') FX(idp),FY(idp),&
         FZ(idp),XX,XY,XZ,YY,YZ,ZZ

         write(444,'(1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9,1X,F15.9)') FX(idp),FY(idp),&
         FZ(idp),XX,XY,XZ,YY,YZ,ZZ

      end do

      close(444)
      close(ftg1,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file FTG.PRED, logical unit ',&
                    ftg1,'.Stooooop in CALGRADIO!'
         STOP
      end if

      write(*,*)'    Synthetic FTG data stored in FTG.PRED file'
      write(inout,*)'    Synthetic FTG data stored in FTG.PRED file'

! enregistrement de la trace du tenseur
      open(444,file='trace.pred',status='REPLACE')
      do idp=jbegin(3),jbegin(3)+nptsgr-1

            iter=0
            if (rtvar(nbn+1) .ne. 0) then
               iter=iter+1
               XX=caldata(idp+nptsgr*(iter-1))
            else
               XX=0   
            endif

            if (rtvar(nbn+2) .ne. 0) then
               iter=iter+1
               XY=caldata(idp+nptsgr*(iter-1))
            else
               XY=0   
            endif

            if (rtvar(nbn+3) .ne. 0) then
               iter=iter+1
               XZ=caldata(idp+nptsgr*(iter-1))
            else
               XZ=0   
            endif

            if (rtvar(nbn+4) .ne. 0) then
               iter=iter+1
               YY=caldata(idp+nptsgr*(iter-1))
            else
               YY=0   
            endif

            if (rtvar(nbn+5) .ne. 0) then
               iter=iter+1
               YZ=caldata(idp+nptsgr*(iter-1))
            else
               YZ=0   
            endif

            if (rtvar(nbn+6) .ne. 0) then
               iter=iter+1
               ZZ=caldata(idp+nptsgr*(iter-1))
            else
               ZZ=0   
            endif

         write(444,'(4f15.3)') FX(idp),FY(idp),FZ(idp),&
            XX+YY+ZZ
      enddo
      close(444)

!=====================================================================
! In case of direct problem, stop the program
!=====================================================================
      if(.not.INV) then
         write(*,*)''
         write(*,*)'End of forward calculation of gravity (FTG) anomaly'
         write(*,*)'Synthetic data stored in file FTG.PRED'
         write(*,*)'I''ve finished my job... BYE!'
         write(*,*)''
         write(*,*)'Number of data calculated (all component): ',jend(3)-jbegin(3)+1
         write(*,*)'Need to divide by 6 to have the real number of data'
         write(*,*)''
         write(inout,*)''
         write(inout,*)'End of forward calculation of gravity (FTG) anomaly'
         write(inout,*)'Synthetic data stored in file FTG.PRED'
         write(inout,*)'I''ve finished my job... BYE!'
         write(inout,*)''
         write(inout,*)'Number of data calculated (all component): ',jend(3)-jbegin(3)+1
         write(*,*)'Need to divide by 6 to have the real number of data'         
         write(inout,*)''
         STOP 'THE FORWARD CALCULATION'
      end if

end subroutine