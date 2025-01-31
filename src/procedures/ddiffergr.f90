!=====================================================================
!=====================================================================
!
!			SUBROUTINE DDIFFERGR
!
!> subroutine calculates the difference between observed and calculated
!! gradiometry data (USEDATA, CALDATA) and uses it to compute the Root Mean
!! Square (RMSGR)
!!
!! Also calculates the difference between previous (PAR0) and actual (PAR)
!! parameters if there is not gravity inversion
!
! This is useful for testing the threshold to stop the iteration. If the
! sum of differences is small enough, we stop.

!=====================================================================
!=====================================================================

SUBROUTINE DDIFFERGR(usedata,caldata,rmsgr,iiter,par,par0,punvar,rovar,&
                     varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
                     diff,dtot,grdtot,h1,bderi,npts,errclgr,rtvar,ncomp)

USE MOD_delim
USE MOD_unit
USE MOD_inv

USE MOD_dsmooth

USE MOD_ddiffer
USE MOD_grread

IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DIFFER
!=====================================================================
integer,intent(in)                               :: iiter
real(kind=8),DIMENSION(:),pointer                :: usedata,caldata
real(kind=8),DIMENSION(:,:),intent(inout)        :: rmsgr,errclgr
real(kind=8),DIMENSION(:),intent(in)             :: par,par0
real(kind=8),DIMENSION(:,:),intent(inout)        :: rovar
real(kind=8),DIMENSION(:),intent(inout)          :: diff,grdtot
real(kind=8),DIMENSION(:),pointer                :: punvar

! pour la partie modèle de densité
integer :: dcont
integer,intent(in)                               :: ismooth
integer,DIMENSION(:),intent(in)                  :: iside
integer,DIMENSION(:,:),intent(in)                :: jside
real(kind=8)                                     :: difs,dsum
real(kind=8),DIMENSION(:),intent(in)             :: smooth
real(kind=8),DIMENSION(:,:),intent(in)           :: xb,yb,zb
real(kind=8),DIMENSION(:),intent(inout)          :: h1
real(kind=8),DIMENSION(:,:),intent(inout)        :: bderi,dtot
real(kind=8),DIMENSION(:),intent(in)             :: varpar

integer,DIMENSION(:),intent(in)      :: npts

real(kind=8),DIMENSION(:),intent(in)    :: rtvar         ! MP36
integer,intent(inout)                   :: ncomp         ! MP36


!=====================================================================
! Declaration of the dummy arguments of DIFFER
!=====================================================================
integer                                :: i,id,j,ig,nbn
real(kind=8)                           :: maxdif,absdif
real(kind=8)                           :: dev,difgr,difd
real(kind=8)                           :: mean,s

real(kind=8),DIMENSION(:),allocatable :: rmsgr2
integer :: iter

real(kind=8),DIMENSION(6)              :: test          ! MP36 


!=====================================================================
! Writting some outputs
!=====================================================================
   write(*,*)''
   write(*,*)'COMPUTING THE DIFFERENCES BETWEEN OBS. AND CALC FTG DATA'
   if(iiter.eq.1) then
      write(inout,*)''
      write(inout,*)'DDIFFERGR CALLED FOR COMPUTING INITIAL SDT DEVIATION'
      write(inout,*)'ITERATION = 0'
   else
      write(inout,*)''
      write(inout,*)'DDIFFERGR CALLED FOR DIFFERENCE BETWEEN OBS. and CALC'
      write(inout,*)'ITERATION = ',iiter
   end if
   
   if(jbegin(3).eq.0) then
      write(*,*)''
      write(*,*)'Inconsistency in DDIFFERGR: data file is &
                 &read, but no data are stored.'
      STOP 'in DDIFFERGR.'
   end if

!=====================================================================
! Differences between observed and calculated FTG DATA
! Attention, à la différence de la GRAVI, il y a cette fois çi 6 composantes.
! Etant donné que le RMS est juste utilisé comme information, et non dans un calcul,
! on peut calculer un rms moyen et un rms par composante, pour voir si telle ou telle composante
! est mieux fittée
!          _________________________________
! RMSGR(2:7) = V ( Sum[(obs(2:7)-calc(2:7))**2] / (Ndata-1) 
! RMSGR(1)= sum RMSGR(2:7) / 6
! rangé dans l'ordre suivant : moyenne / XX / XY / XZ / YY / YZ / ZZ
! VARIANCE = ROVAR = Sum[x - moy]**2/(Ndata - 1)
!=====================================================================
   maxdif = 0.d0
   mean = 0.d0
   s = 0.d0

   if(.not.INVD) then
      rovar(iiter,:) = 0.0d0
   endif

! dimension de rovar ? à la base c'est ALLOCATE (rovar(iter+1,2)) 
! colonne 1 : les données GRAVI
! colonne 2 : les paramètres
! ajout colonne 3 => les données gradio comme ça pas besoin de changer les autres codes
!	qui sont déjà callé avec les colonnes 1 et 2

! calcul du RMS pour l'ensemble des composantes
   do i=jbegin(3),jend(3)
      mean = mean + caldata(i)
      diff(i) = usedata(i) - caldata(i)
      absdif = dabs(diff(i))
      if(absdif.gt.maxdif) maxdif = absdif
      rmsgr(iiter,1) = rmsgr(iiter,1) + diff(i)*diff(i)
   end do

   if(jend(3).gt.jbegin(3)) then
      mean = mean/(jend(3)-jbegin(3)+1)
      rmsgr(iiter,1) = dsqrt(rmsgr(iiter,1)/(jend(3)-jbegin(3)))
   else
      write(*,*)'    Can''t calculate the standard deviation'
      STOP ' in DDIFFERGR: jend(3)-jbegin(3)=0 !'
   end if

! calcul d'un RMS pour chaque composante'
allocate(rmsgr2(ncomp))
rmsgr2(:)=0

      nbn=0     ! pour prendre en compte le nombre de méthode que l'on inverse en plus de la gradio
      if (INVD) then
        nbn=nbn+1
      endif
      if (INVV) then
        nbn=nbn+1
      endif

    do j=1,ncomp
       do i=jbegin(3)+(j-1)*((jend(3)-jbegin(3)+1)/ncomp),jbegin(3)+(j)*((jend(3)-jbegin(3)+1)/ncomp)-1 
          rmsgr(iiter,j+1) = rmsgr(iiter,j+1) + diff(i)*diff(i)
       end do
          rmsgr(iiter,j+1) = dsqrt(rmsgr(iiter,j+1)/(((jend(3)-jbegin(3)+1)/ncomp)-1)) 

    end do

    do j=1,ncomp
      rmsgr2(j)=rmsgr(iiter,j+1)
    enddo
    
    iter=0
    do j=1,6
      if (rtvar(nbn+j) .ne. 0) then
        rmsgr(iiter,j+1)=rmsgr2(j-iter)
      else
        iter=iter+1
        rmsgr(iiter,j+1)=0
      endif
    enddo

   do i=jbegin(3),jend(3)
      s = caldata(i) - mean
      rovar(iiter,3) = rovar(iiter,3) + s*s
   end do

   if(jend(3).gt.jbegin(3)) then
      rovar(iiter,3) = rovar(iiter,3)/(jend(3)-jbegin(3))
      write(inout,'(5x,''Variance of the calculated FTG:'',22(''.''),&
           &f12.3)') rovar(iiter,3)
      write(inout,'(5x,''Mean of the calculated FTG:'',26(''.''),&
           &f12.3)') mean
   else
      write(*,*)'    Can''t calculate the variance'
      STOP ' in DDIFFERGR : jend(3)-jbegin(3)=0 !'
   end if

   write(inout,'(5x,''Max. difference between obs. and calc. FTG &
        &data:'',5(''.''),f12.3,'' E'')') maxdif

   write(*,'(5x,''Max. difference between obs. and calc. FTG &
        &data:'',5(''.''),f12.3,'' E'')') maxdif

!=====================================================================
! Differences between observed and calculated data points
!=====================================================================
   difgr = 0.d0
   ig = 0
   do i=jbegin(3),jend(3)
      ig=ig+1
      difgr = difgr + diff(i)*diff(i)*punvar(i)
   end do
   difgr = difgr/ig

! les 3 prochaines parties ne sont calculées que si il n'y a pas d'inversion gravi
!  sinon pas besoin c est deja fait dans la partie gravi (DDIFFER)
!=====================================================================
! Differences between previous and actual parameters
!=====================================================================
   if(.not.INVD) then
      id = 0
      dev = 0.d0
      difd = 0.d0
      mean=0.d0
      s=0.d0
      do i=ibegin(1),iend(1)
         mean=mean+par(i)
         dev = par(i) - par0(i)
         id = id+1
         difd = difd + dev*dev*varpar(i)
      end do
      difd = difd/id
      mean = mean/id
   
!=====================================================================
! Variances of the actual parameters
!=====================================================================
      do i=ibegin(1),iend(1)
         s=par(i)-mean
         rovar(iiter,2)=rovar(iiter,2) + s*s
      end do

      rovar(iiter,2)=rovar(iiter,2)/(iend(1)-ibegin(1))

!=====================================================================
! If Smoothing constraint is present and iteration .gt. 1,
! it returns the difference between parameters of contiguous blocks
! weighted by their distance and smoothing factor.
!=====================================================================
      if(ismooth.ne.0 .and. iiter.gt.1) then
         CALL DSMOOTH(smooth,iside,jside,xb,yb,par,dsum,dcont,&
                      bderi,h1)
         if(dcont.ne.0) difs=dsum/dcont
      else
         difs = 0.d0
      end if

      dtot(iiter,1) = 0
      dtot(iiter,2) = difd
      dtot(iiter,3) = difs
      dtot(iiter,4) = difd+difs
   endif

      grdtot(iiter)=difgr

      write(*,'(5x,''RMS of obs-calc FTG data (mean of all components):'',28(''.''),&
           &f12.3,'' E'')') rmsgr(iiter,1)
      write(*,'(5x,''    - for XX components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,2)
      write(*,'(5x,''    - for XY components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,3)
      write(*,'(5x,''    - for XZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,4)
      write(*,'(5x,''    - for YY components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,5)
      write(*,'(5x,''    - for YZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,6)
      write(*,'(5x,''    - for ZZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,7)

      write(inout,'(5x,''RMS of obs-calc FTG data (mean of all components) :'',28(''.''),&
           &f12.3,'' E'')') rmsgr(iiter,1)
      write(inout,'(5x,''    - for XX components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,2)
      write(inout,'(5x,''    - for XY components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,3)
      write(inout,'(5x,''    - for XZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,4)
      write(inout,'(5x,''    - for YY components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,5)
      write(inout,'(5x,''    - for YZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,6)
      write(inout,'(5x,''    - for ZZ components:'',28(''.''),f12.3,'' E'')') rmsgr(iiter,7)

      write(inout,'(5x,''Sum of mean diff. for gravity/density:'',19(''.''),&
           &f12.3)') dtot(iiter,4)

      write(inout,*)'    whereof:'

      write(inout,'(13x,''Mean diff. between obs. and calc. FTG data:..&
           &'',f12.3,'' E'')') difgr

!      write(inout,'(13x,''toto'',f12.3)') dtot(iiter,4)      

      write(inout,'(13x,''Mean diff. between act. and prev. &
           &parameters:'',4(''.''),f12.3)') difd

      write(inout,'(5x,''Roughness from smoothing constraint:'',21(''.''),&
           &f12.3)') difs

! calcul de l'erreur comme Clémence sans la partie pondération de chaque méthode (alpha)
! une erreur pour l'ensemble des données
do i=jbegin(3),jend(3)
  diff(i)=usedata(i)-caldata(i)
  errclgr(iiter,1) = errclgr(iiter,1) + ((diff(i)*diff(i))/((1/punvar(i))*((jend(3)-jbegin(3)+1))))

enddo

END SUBROUTINE