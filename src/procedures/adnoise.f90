!=====================================================================
!=====================================================================
!   SUBROUTINE ADNOISE
!=====================================================================
!=====================================================================
!> This SUBROUTINE finds gaussian random numbers and adds it to
!! synthetic datum.
!!
!! Taken from Steck and Prothero programm
!
! Called by CALGRA
! Calls: RAN1(internal function)
!=====================================================================
    !INCLUDE 'INTERF/MOD_ran1.f'

    SUBROUTINE ADNOISE(syndata,signois,idum)

      USE MOD_ran1
      
      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of ADNOISE
!=====================================================================
      integer,intent(inout)       :: idum

      real(kind=8),intent(in)     :: signois
      real(kind=8),intent(inout)  :: syndata
!=====================================================================
! Declaration of the dummy arguments of ADNOISE
!=====================================================================
      integer                     :: iset
      real(kind=8)                :: noise,fac,gset,rsq,v1,v2,rana,ranb

      DATA iset/0/
!=====================================================================
! Initialization of some variables
!=====================================================================
      rsq = 0.d0

      if(iset.eq.0) then
         do while(rsq.ge.1..or.rsq.eq.0.)
            CALL RAN1(idum,rana)
            v1=2.d0*rana - 1.d0
            CALL RAN1(idum,ranb)
            v2=2.d0*ranb - 1.d0
            rsq=v1**2+v2**2
         end do
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         noise=v2*fac
         iset=1
      else
         noise=gset
         iset=0
      end if
!=====================================================================
! Add the random gaussian noise to the synthetic datum
!=====================================================================
      syndata = syndata + noise*signois


    END SUBROUTINE ADNOISE
!=====================================================================
!=====================================================================
!   SUBROUTINE RAN1
!=====================================================================
!=====================================================================
!> This subroutine returns a random number between 0 and 1. Taken from the
!! Numerical Recipes (ran1).
!
! "Minimal" random number generator of Park and Miller with Bays-Durham 
! shuffle and added safeguards. Returns a uniform random deviate 
! between 0.0 and 1.0 (exclusive of the endpoint values). Call with 
! idum a negative integer to initialize; thereafter, do not alter idum 
! between successive deviates in a sequence. RNMX should approximate 
! the largest oating value that is less than 1.
!
! Called by ADNOISE
! calls none
!=====================================================================
      SUBROUTINE RAN1(idum,ran)
      
        IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of RAN1
!=====================================================================
        integer,intent(inout)    :: idum

        real(kind=8),intent(out) :: ran
!=====================================================================
! Declaration of the dummy arguments of RAN1
!=====================================================================
        integer                 :: j,k,iy
        integer,PARAMETER       :: IA=16807,IM=2147483647,IQ=127773,IR=2836
        integer,PARAMETER       :: NTAB=32,NDIV=1+(IM-1)/NTAB
        integer,DIMENSION(NTAB) :: iv

        real(kind=8),PARAMETER  :: AM=1./IM,EPS=1.2e-7,RNMX=1.-EPS

        DATA iv/NTAB*0/,iy/0/
!=====================================================================
! Initialize to be sure to prevent idum = 0
!=====================================================================
        if(idum.le.0.or.iy.eq.0) then 
           idum=max0(-idum,1)
!=====================================================================
! Load the shuffle table (after 8 warm-ups)
!=====================================================================
           do j=NTAB+8,1,-1
              k=idum/IQ 
              idum=IA*(idum-k*IQ)-IR*k 
              if (idum.lt.0) idum=idum+IM 
              if (j.le.NTAB) iv(j)=idum 
           end do
           iy=iv(1) 
        end if
!=====================================================================
! Start here when not initializing.
!=====================================================================
        k=idum/IQ
!=====================================================================
! Computes idum=mod(IA*idum,IM) without overflows by Schrage's method
!=====================================================================
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
!=====================================================================
! Will be in the range 1:NTAB.
!=====================================================================
        j=1+iy/NDIV
!=====================================================================
! Output previously stored value and refill the shuffle table
!=====================================================================
        iy=iv(j)
        iv(j)=idum
!=====================================================================
! Because users don't expect endpoint values.
!=====================================================================
        ran=dmin1(AM*iy,RNMX)

      END SUBROUTINE RAN1


