!=====================================================================
!=====================================================================
!   SUBROUTINE BCALC
!
!> This subroutine computes the B-coefficient of the linear relation
!! between density and velocity.
!!
!! The mean coefficient is calculated from the densities and velocities
!! only constrained by 3D raytracing.
!!
!! A coefficient of correlation is calculated between velocity and
!! density to indicate if the linear relation is far or not.
!! (i.e. a good correlation mean a good linear relationship between
!! both parameters, a bad one means that there is a big dispersion
!! between the parameters).
!
!=====================================================================
!=====================================================================
! 
! Called by MAIN
! calls none
!=====================================================================

    SUBROUTINE BCALC(nbod,ilay,par,ibove,iiter,thrhold)

      USE MOD_delim
      USE MOD_layer
      USE MOD_unit
      USE MOD_inv
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of BCALC
!=====================================================================
      integer,intent(in)                      :: nbod,iiter
      integer,DIMENSION(:),intent(in)         :: ilay,ibove

      real(kind=8),intent(in)                 :: thrhold
      real(kind=8),DIMENSION(:),intent(inout) :: par
!=====================================================================
! Declaration of the dummy arguments of BCALC
!=====================================================================
      logical                                :: constrain

      integer                                :: ier,m,nunode,hitbloc
      integer                                :: i,j,vpn,ii,nl,ib
      integer,DIMENSION(:),ALLOCATABLE       :: nrho

      real(kind=8)                           :: bcoef,corre,parho
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: sumrho,sumv,parvel
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: sumrho2,sumvr,sumb
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: covv,covrho,covrv
!=====================================================================
! Initializatio of some variables and arrays
!=====================================================================

      ALLOCATE (parvel(nbod))
      ALLOCATE (sumrho(nlayer))
      ALLOCATE (sumrho2(nlayer))
      ALLOCATE (sumv(nlayer))
      ALLOCATE (sumvr(nlayer))
      ALLOCATE (sumb(nlayer))
      ALLOCATE (covv(nlayer))
      ALLOCATE (covrv(nlayer))
      ALLOCATE (covrho(nlayer))
      ALLOCATE (nrho(nlayer))

      parvel(:) = 0.d0
      sumrho(:) = 0.d0
      sumrho2(:) = 0.d0
      sumv(:) = 0.d0
      sumvr(:) = 0.d0
      sumb(:) = 0.d0
      covrv(:) = 0.d0
      covrho(:) = 0.d0
      covv(:) = 0.d0
      nrho(:) = 0.d0

      write(*,*)''
      write(*,*)'FORWARD CALCULATION OF B-COEFFICIENT FROM BCALC'
      write(inout,*)''
      write(inout,*)'FORWARD CALCULATION OF B-COEFFICIENT FROM BCALC'
!=====================================================================
! Opening the blonod file (BLOC.NOD) for constrained blocks
!=====================================================================
      open(blocnod,file='bloc.nod',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file BLOC.NOD, logical unit ',&
                   blocnod,'.Stooooop in BCALC!'
         STOP 'in BCALC!'
      end if
!=====================================================================
! Search only for the constrained density blocks and the mean velocity
! associated
!=====================================================================
      do i=ibegin(1),iend(1)

         read(blocnod,*) m,nunode
         ii=i-ibegin(1)+1
         nl=ilay(ii)
         if(ii.ne.m) then
            write(*,*)'Inconsistency in BCALC. Stooooop!'
            close(blocnod)
            STOP
         end if

         constrain = .FALSE.

         if(nunode.ne.0) then
            hitbloc=0
            do j=1,nunode
               read(blocnod,*) vpn
               if(ibove(vpn).ne.0) then
                  constrain = .TRUE.
                  parvel(i)=parvel(i) + par(ibove(vpn))
                  hitbloc=hitbloc+1
               end if
            end do
         end if
         
         if(constrain) then
            nrho(nl) = nrho(nl) + 1
            parvel(i)=parvel(i)/hitbloc
            sumrho(nl) = sumrho(nl) + par(i)
            sumrho2(nl)=sumrho2(nl) + par(i)*par(i)
            sumv(nl)=sumv(nl)+parvel(i)
            sumvr(nl)=sumvr(nl)+parvel(i)*par(i)
            if(par(i).ne.0) then
               sumb(nl) = sumb(nl) + parvel(i)/par(i)
            end if
         end if

      end do
!=====================================================================
! Computing the mean density and velocity of each layer
!=====================================================================
      do nl=1,nlayer
         if(nrho(nl).ne.0) then
            sumrho(nl)=sumrho(nl)/nrho(nl)
            sumrho2(nl)=sumrho2(nl)/nrho(nl)
            sumv(nl)=sumv(nl)/nrho(nl)
            sumvr(nl)=sumvr(nl)/nrho(nl)
            sumb(nl)=sumb(nl)/nrho(nl)
         end if
      end do
!=====================================================================
! Compute the Covariance for each layer and store the velocity density
! values for each block in CORREL.OUT file
!=====================================================================
      rewind(blocnod)

      do i=ibegin(1),iend(1)
         read(blocnod,*) m,nunode
         ii=i-ibegin(1)+1
         nl=ilay(ii)
         if(ii.ne.m) then
            write(*,*)'Inconsistency in BCALC. Stooooop!'
            close(blocnod)
            STOP
         end if

         constrain = .FALSE.

         if(nunode.ne.0) then
            hitbloc=0
            do j=1,nunode
               read(blocnod,*) vpn
               if(ibove(vpn).ne.0) then
                  constrain = .TRUE.
               end if
            end do
         end if

         if(constrain) then
            covrv(nl)=covrv(nl)+(par(i) - sumrho(nl))*(parvel(i) - sumv(nl))
            covrho(nl) = covrho(nl) + (par(i) - sumrho(nl))**2
            covv(nl) = covv(nl) + (parvel(i) - sumv(nl))**2
            write(correl,'(i2,1x,i2,f10.4,1x,f10.4)') iiter,nl,&
                 par(i),parvel(i)
         end if
      end do
      

!=====================================================================
! Closing the file blocnod
!=====================================================================
      close(blocnod,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file BLOC.NOD, logical unit ',&
                   blocnod,'.Stooooop in BCALC!'
         STOP
      end if

!=====================================================================
! Testing calculated B-coeff.
!
! Mean B-coeff = (mean V*rho)/(mean rho**2)
! Correlation coefficient is the Pearson's one:
!             S(xi-xm)(yi-ym)
! Corr = -------------------------
!         _______________________
!        V (S(xi-xm)2*S(yi-ym)2)
!
! If calculated B-coeff do not coincide, the corresponding parameter
! is changed to the real one
!=====================================================================

      write(bcoeff,'(''ITERATION no '',i2)') iiter

      do nl=1,nlayer
         ib=ibegin(3)-1+nl

         if(sumrho(nl).ne.0 .and. covv(nl).ne.0 .and. covrho(nl).ne.0) then
            bcoef=sumvr(nl)/sumrho2(nl)
            corre=covrv(nl)/DSQRT(covv(nl)*covrho(nl))
         else
            bcoef=0.d0
            corre=0.d0
         end if

!=====================================================================
! May be I have to change that... Instead of corre.ne.0, I can enter
! a given threshold, so that the B will only be change if this is passed.
! For ex. if corre.gt.0.80, then par(ib)=bcoef. else, it's the same
! value.
!=====================================================================
!         if(corre.ne.0.) then
         if(DABS(corre).gt.thrhold) then
            par(ib)=bcoef
            write(inout,'(5x,''LAYER '',i2,'': B= '',f6.2,'' or '',f6.2)') &
                 nl,bcoef,sumb(nl)
            write(inout,'(15x,''Correlation Coef= '',f6.3)') corre
            write(inout,'(15x,''Number of points= '',i6)') nrho(nl)
         else
            write(inout,'(5x,''LAYER '',i2,'': B-coeff not directly defined; &
                 &taking the inversion result'')') nl
            write(inout,'(15x,''Correlation Coef ('',f6.3,'') < thrhold&
	    & ('',f5.3,'')'')') corre,thrhold
         end if

         write(bcoeff,'(i3,1x,f10.4,f10.4)') nl,par(ib),corre

      end do
      
      DEALLOCATE (parvel)
      DEALLOCATE (sumrho)
      DEALLOCATE (sumrho2)
      DEALLOCATE (sumv)
      DEALLOCATE (sumvr)
      DEALLOCATE (sumb)
      DEALLOCATE (covrho)
      DEALLOCATE (covv)
      DEALLOCATE (covrv)
      DEALLOCATE (nrho)

    END SUBROUTINE BCALC
