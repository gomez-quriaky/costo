!=====================================================================
!=====================================================================
!
!			SUBROUTINE DDIFFER
!
!> subroutine calculates the difference between observed and calculated
!! density data (USEDATA, CALDATA) and uses it to compute the Root Mean
!! Square (RMSG)
!!
!! Also calculates the difference between previous (PAR0) and actual (PAR)
!! parameters
!
! This is useful for testing the threshold to stop the iteration. If the
! sum of differences is small enough, we stop.
!
! 04/07/2001: Fixing a bug in calling the DSMOOTH subroutine. It needs to
!             have the h1 and bderi dimension because they are optional
!             parameters, and thus should be present in the declaration of 
!             DDIFFER.
!
! Calls: DSMOOTH
! called by MAIN
!
!=====================================================================
!=====================================================================

    SUBROUTINE DDIFFER(usedata,caldata,rmsg,iiter,par,par0,punvar,rovar,&
                       varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
                       diff,dtot,h1,bderi,errclg)
        
      USE MOD_delim
      USE MOD_unit

      USE MOD_dsmooth

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DIFFER
!=====================================================================
      integer,intent(in)                               :: iiter,ismooth
      integer,DIMENSION(:),intent(in)                  :: iside
      integer,DIMENSION(:,:),intent(in)                :: jside

      real(kind=8),DIMENSION(:),intent(inout)          :: rmsg,errclg
      real(kind=8),DIMENSION(:),intent(in)             :: par,par0,varpar
      real(kind=8),DIMENSION(:),intent(in)             :: smooth
      real(kind=8),DIMENSION(:),intent(inout)          :: h1
      real(kind=8),DIMENSION(:,:),intent(inout)        :: bderi
      real(kind=8),DIMENSION(:,:),intent(in)           :: xb,yb,zb
      real(kind=8),DIMENSION(:,:),intent(inout)        :: dtot,rovar
      real(kind=8),DIMENSION(:),intent(inout)          :: diff
      real(kind=8),DIMENSION(:),pointer                :: punvar
      real(kind=8),DIMENSION(:),pointer                :: usedata,caldata
!=====================================================================
! Declaration of the dummy arguments of DIFFER
!=====================================================================
      integer                                :: i,id,ig,dcont

      real(kind=8)                           :: maxdif,absdif,dsum
      real(kind=8)                           :: dev,difg,difd,difs
      real(kind=8)                           :: mean,s
!=====================================================================
! Writting some outputs
!=====================================================================
      write(*,*)''
      write(*,*)'COMPUTING THE DIFFERENCES BETWEEN OBS. AND CALC GRAVITY'
      if(iiter.eq.1) then
         write(inout,*)''
         write(inout,*)'DDIFFER CALLED FOR COMPUTING INITIAL SDT DEVIATION'
         write(inout,*)'ITERATION = 0'
      else
         write(inout,*)''
         write(inout,*)'DDIFFER CALLED FOR DIFFERENCE BETWEEN OBS. and CALC'
         write(inout,*)'ITERATION = ',iiter
      end if

      if(jbegin(1).eq.0) then
         write(*,*)''
         write(*,*)'Inconsistency in DDIFFER: data file is &
                    &read, but no data are stored.'
         STOP 'in DDIFFER.'
      end if
!=====================================================================
! Differences between observed and calculated GRAVITY DATA
! 
!          _________________________________
! RMSG = V ( Sum[(obs-calc)**2] / (Ndata-1) 
! 
! VARIANCE = ROVAR = Sum[x - moy]**2/(Ndata - 1)
!=====================================================================
      maxdif = 0.d0
      mean = 0.d0
      s = 0.d0
      rovar(iiter,:) = 0.0d0

      do i=jbegin(1),jend(1)
         mean = mean + caldata(i)
         diff(i) = usedata(i) - caldata(i)
         absdif = dabs(diff(i))
         if(absdif.gt.maxdif) maxdif = absdif
         rmsg(iiter) = rmsg(iiter) + diff(i)*diff(i)
      end do

      if(jend(1).gt.jbegin(1)) then
         mean = mean/(jend(1)-jbegin(1)+1)
         rmsg(iiter) = dsqrt(rmsg(iiter)/(jend(1)-jbegin(1)))
      else
         write(*,*)'    Can''t calculate the standard deviation'
         STOP ' in DDIFFER: jend(1)-jbegin(1)=0 !'
      end if

      do i=jbegin(1),jend(1)
         s = caldata(i) - mean
         rovar(iiter,1) = rovar(iiter,1) + s*s
      end do
      if(jend(1).gt.jbegin(1)) then
         rovar(iiter,1) = rovar(iiter,1)/(jend(1)-jbegin(1))
         write(inout,'(5x,''Variance of the calculated gravity:'',22(''.''),&
              &f12.3)') rovar(iiter,1)
         write(inout,'(5x,''Mean of the calculated gravity:'',26(''.''),&
              &f12.3)') mean
      else
         write(*,*)'    Can''t calculate the variance'
         STOP ' in DDIFFER: jend(1)-jbegin(1)=0 !'
      end if
      
      write(inout,'(5x,''Max. difference between obs. and calc. gravity &
           &data:'',5(''.''),f12.3,'' mGal'')') maxdif

      write(*,'(5x,''Max. difference between obs. and calc. gravity &
           &data:'',5(''.''),f12.3,'' mGal'')') maxdif
!=====================================================================
! Differences between observed and calculated data points
!=====================================================================
      difg = 0.d0
      ig = 0
      do i=jbegin(1),jend(1)
         ig=ig+1
         difg = difg + diff(i)*diff(i)*punvar(i)
      end do
      difg = difg/ig
!=====================================================================
! Differences between previous and actual parameters
!=====================================================================
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
!=====================================================================
! It's time to write some outputs in inout logical unit (parameter.out)
!=====================================================================
      dtot(iiter,1) = difg
      dtot(iiter,2) = difd
      dtot(iiter,3) = difs
      dtot(iiter,4) = difg+difd+difs

      write(*,'(5x,''RMS of obs-calc gravity data:'',28(''.''),&
           &f12.3,'' mGal'')') rmsg(iiter)

      write(inout,'(5x,''RMS of obs-calc gravity data:'',28(''.''),&
           &f12.3,'' mGal'')') rmsg(iiter)
      write(inout,'(5x,''Sum of mean diff. for gravity/density:'',19(''.''),&
           &f12.3)') dtot(iiter,4)
      write(inout,*)'    whereof:'
      write(inout,'(13x,''Mean diff. between obs. and calc. gravity data:..&
           &'',f12.3,'' mGal'')') difg
      write(inout,'(13x,''Mean diff. between act. and prev. &
           &parameters:'',4(''.''),f12.3)') difd
      write(inout,'(5x,''Roughness from smoothing constraint:'',21(''.''),&
           &f12.3)') difs

! MP34 : calcul de l'erreur comme C. Basuyau et al., 2011
!print*,jbegin(1),jend(1)
      do i=jbegin(1),jend(1)
        diff(i) = usedata(i) - caldata(i)
        errclg(iiter)=errclg(iiter)+((diff(i)*diff(i))/((1/punvar(i))*((jend(1)-jbegin(1)+1))))
!        print*,i,jbegin(1),jend(1),diff(i),errclg(iiter)
      enddo

!pause

    END SUBROUTINE DDIFFER
