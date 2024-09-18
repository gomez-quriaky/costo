!=====================================================================
!=====================================================================
!
!			SUBROUTINE VDIFFER
!
!> subroutine calculates the difference between observed and calculated
!! delaytime data (USEDATA, CALDATA) and uses it to compute the Root Mean
!! Square (RMST)
!!
!! Also calculates the difference between previous (PAR0) and actual (PAR)
!! velocity parameters
!
! 04/07/2001: Fixing a bug in calling the VSMOOTH subroutine. It needs to
!             have the h1 and bderi dimensions because they are optional
!             parameters, and thus should be present in the declaration of 
!             VDIFFER.
!
! Calls: VSMOOTH
! called by MAIN
!
!=====================================================================
!=====================================================================

    SUBROUTINE VDIFFER(usedata,caldata,rmst,iiter,par,par0,punvar,&
                       varpar,vxnodes,vynodes,ismooth,smooth,ivside,&
                       jvside,ibove,diff,dtot,velvar,h1,bderi,errclt)
        
      USE MOD_delim
      USE MOD_inv
      USE MOD_unit
      USE MOD_size
      
      USE MOD_vsmooth

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of VDIFFER
!=====================================================================
      integer,intent(in)                       :: iiter,ismooth
      integer,DIMENSION(:),intent(in)          :: ivside,ibove
      integer,DIMENSION(:,:),intent(in)        :: jvside

      real(kind=8),DIMENSION(:),intent(inout)  :: rmst,errclt
      real(kind=8),DIMENSION(:),intent(in)     :: par,par0,varpar
      real(kind=8),DIMENSION(:),intent(in)     :: smooth
      real(kind=8),DIMENSION(:),intent(inout)  :: h1
      real(kind=8),DIMENSION(:),intent(in)     :: vxnodes,vynodes
      real(kind=8),DIMENSION(:),intent(inout)  :: diff
      real(kind=8),DIMENSION(:,:),intent(inout):: bderi
      real(kind=8),DIMENSION(:,:),intent(inout):: dtot,velvar
      real(kind=8),DIMENSION(:),pointer        :: punvar
      real(kind=8),DIMENSION(:),pointer        :: usedata,caldata
!=====================================================================
! Declaration of the dummy arguments of VDIFFER
!=====================================================================
      integer                                :: i,iv,idt,vcont

      real(kind=8)                           :: maxdif,absdif,vsum
      real(kind=8)                           :: dev,dift,difv,difs
      real(kind=8)                           :: s,mean
!=====================================================================
! Writting some outputs
!=====================================================================
      write(*,*)''
      write(*,*)'COMPUTING THE DIFFERENCES BETWEEN OBS. AND CALC VELOCITY'
      
      if(iiter.eq.1) then
         write(inout,*)''
         write(inout,*)'VDIFFER CALLED FOR COMPUTING INITIAL SDT DEVIATION'
         write(inout,*)'ITERATION = 0'
      else
         write(inout,*)''
         write(inout,*)'VDIFFER CALLED FOR DIFFERENCE BETWEEN OBS. and CALC'
         write(inout,*)'ITERATION = ',iiter
      end if

      maxdif = 0.d0
      mean = 0.d0
      s = 0.d0
      velvar(iiter,:) = 0.d0

      if(jbegin(2).eq.0) then
         write(*,*)''
         write(*,*)'Inconsistency in VDIFFER: data file is &
                    &read, but no data are stored. Stooooop!'
         STOP
      end if
!=====================================================================
! Sum of the differences for each data point and
! Calcul of the standard deviation for delaytime data:
!          ___________________________________
! RMST = V ( Sum[(obs-calc)**2] / (Ndata-1) )
!
! VARIANCE = VELVAR = Sum[x - moy]**2 / (Ndata-1)
!=====================================================================
      do i=jbegin(2),jend(2)
         mean = mean + caldata(i)
         diff(i) = usedata(i) - caldata(i)
         absdif = dabs(diff(i))
         if(absdif.gt.maxdif) maxdif = absdif
         rmst(iiter) = rmst(iiter) + diff(i)*diff(i)

!         print*,jbegin(2),jend(2)
!         print*,usedata(i),caldata(i),diff(i)

      end do

      if(jend(2).gt.jbegin(2)) then
         mean = mean / (jend(2)-jbegin(2) +1)
         rmst(iiter) = dsqrt(rmst(iiter)/(jend(2)-jbegin(2)))
      else
         write(*,*)'    Can''t calculate the standard deviation'
         STOP ' in VDIFFER: jend(2)-jbegin(2)=0 !'
      end if

      do i=jbegin(2),jend(2)
         s = caldata(i) - mean
         velvar(iiter,1) = velvar(iiter,1) + s*s
      end do
      if(jend(2).gt.jbegin(2)) then
         velvar(iiter,1) = velvar(iiter,1)/(jend(2)-jbegin(2))
         write(inout,'(5x,''Variance of the calculated delay times:'',&
              &18(''.''),f12.3)') velvar(iiter,1)
         write(inout,'(5x,''Mean of the calculated delay times:'',22(''.''),&
              &f12.3)') mean
      else
         write(*,*)'    Can''t calculate the variance'
         STOP ' in VDIFFER: jend(2)-jbegin(2)=0 !'
      end if

      write(*,'(5x,''Max. difference between obs. and calc. delay &
           &time data:..'',f12.3,'' s'')') maxdif

      write(inout,'(5x,''Max. difference between obs. and calc. delay &
           &time data:..'',f12.3,'' s'')') maxdif
!=====================================================================
! Differences between observed and calculated data points
!=====================================================================
      idt = 0
      dift = 0.d0
      do i=jbegin(2),jend(2)
         idt=idt+1
         dift = dift + diff(i)*diff(i)*punvar(i)
      end do
      dift = dift/idt
!=====================================================================
! Differences between previous and actual parameters
!=====================================================================
      iv = 0
      dev = 0.d0
      difv = 0.d0
      mean = 0.d0
      s = 0.d0
      do i=ibegin(2),iend(2)
         mean = mean + par(i)
         dev = par(i) - par0(i)
         iv = iv+1
         difv = difv + dev*dev*varpar(i)
      end do
      difv = difv/iv
      mean = mean/iv
!=====================================================================
! Variance of the actual parameters
!=====================================================================
      do i=ibegin(2),iend(2)
         s=par(i)-mean
         velvar(iiter,2)=velvar(iiter,2) + s*s
      end do

      velvar(iiter,2)=velvar(iiter,2)/(iend(2)-ibegin(2))

!=====================================================================
! If Smoothing constraint is present and iteration .gt. 1
!=====================================================================
      if(ismooth.ne.0 .and. iiter.gt.1) then
         CALL VSMOOTH(smooth,ivside,jvside,par,ibove,vxnodes,&
                      vynodes,vsum,vcont,bderi,h1)
         if(vcont.ne.0) difs=vsum/vcont
      else
         difs = 0.d0
      end if
!=====================================================================
! It's time to write some outputs in inout logical unit
!=====================================================================
      dtot(iiter,1) = dift
      dtot(iiter,2) = difv
      dtot(iiter,3) = difs
      dtot(iiter,4) = dift+difv+difs

      write(*,'(5x,''RMS of obs-calc delay time data:'',25(''.''),f12.3,&
           &'' s'')') rmst(iiter)

      write(inout,'(5x,''RMS of obs-calc delay time data:'',25(''.''),f12.3,&
           &'' s'')') rmst(iiter)
      write(inout,'(5x,''Sum of mean diff. for delaytime/velocity:'',&
           &16(''.''),f12.3)') dtot(iiter,4)
      write(inout,*)'    whereof:'
      write(inout,'(13x,''Mean diff. between obs. and calc. time data:'',&
           &5(''.''),f12.3,'' s'')') dift
      write(inout,'(13x,''Mean diff. between act. and prev. parameters:'',&
           &4(''.''),f12.3)') difv
      write(inout,'(5x,''Roughness from smoothing constraint:'',21(''.''),&
           &f12.3)') difs


      do i=jbegin(2),jend(2)
        diff(i) = usedata(i) - caldata(i)
        errclt(iiter)=errclt(iiter)+((diff(i)*diff(i))/((1/punvar(i))*((jend(2)-jbegin(2)+1))))
      enddo


    END SUBROUTINE VDIFFER
