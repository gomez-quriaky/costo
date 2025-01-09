!=====================================================================
!=====================================================================
!       SUBROUTINE DSMOOTH
!
!> This subroutine computes the smoothing constrained for BDERI matrix
!! for the density bodies, term Cs**(-1).Ds in BLDMAT.
!!
!! It returns the mean square difference between two contiguous blocks
!! parameters, weighted by their distance and the smoothing factor.
!
! Calls none
! Called by MAIN, DDIFFER
!
!=====================================================================
!=====================================================================
    SUBROUTINE DSMOOTH(smooth,iside,jside,xb,yb,par,dsum,&
                         dcont,bderi,h1,ismooth)
        
      USE MOD_delim

      IMPLICIT NONE

!=====================================================================
! Declaration of in/out arguments of DSMOOTH
!=====================================================================
      integer,intent(out)                                :: dcont
      integer,OPTIONAL,intent(in)                        :: ismooth
      integer,DIMENSION(:),intent(in)                    :: iside
      integer,DIMENSION(:,:),intent(in)                  :: jside

      real(kind=8),intent(out)                           :: dsum
      real(kind=8),DIMENSION(:,:),intent(in)             :: xb,yb
      real(kind=8),DIMENSION(:),intent(in)               :: smooth,par
      real(kind=8),DIMENSION(:),intent(inout)            :: h1
      real(kind=8),DIMENSION(:,:),intent(inout)          :: bderi
!=====================================================================
! Declaration of dummy arguments of DSMOOTH
!=====================================================================
      integer                 :: i,id,j,js,jj

      real(kind=8)            :: fs,xmi1,ymi1,xmi2,ymi2,dx,dy,dev
      real(kind=8)            :: deltarho

!=====================================================================
! Initialization of some variables
! fs is the diagonal element of the Cs**(-1) matrix
!=====================================================================
      dsum = 0.d0
      dcont = 0

      if(smooth(1).ne.0) then
         fs = 1.d0/(smooth(1)*smooth(1))
!=====================================================================
! Loop over density bodies
!=====================================================================
         i=0
         do id=ibegin(1),iend(1)
            i = i+1
            if(i.gt.size(iside)) then
               write(*,*)'     Out of range subscript for iside'
               STOP 'in DSMOOTH'
            end if
!=====================================================================
! compute the midle of the body
!=====================================================================
            xmi1=(xb(i,1)+xb(i,2))/2.d0
            ymi1=(yb(i,1)+yb(i,2))/2.d0
!=====================================================================
! iside(i) = total number of contiguous constrained blocks
! jside(i,4) = corresponding number of the contiguous blocks
! jj = number of the parameter which corresponds to jside
!=====================================================================
            if(iside(i).gt.0) then

               do j=1,iside(i)
                  
                  js=jside(i,j)
                  jj=ibegin(1)+js-1
                  
                  if(js.gt.size(xb,1)) then
                     write(*,*)'    Out of range subscript for xb'
                     STOP 'in DSMOOTH'
                  end if
                  xmi2=(xb(js,1)+xb(js,2))/2.d0
                  ymi2=(yb(js,1)+yb(js,2))/2.d0

                  dx=xmi2-xmi1
                  dy=ymi2-ymi1

                  dev=fs/(dx*dx+dy*dy)

                  deltarho=par(jj)-par(id)
!=====================================================================
! if not h1, smoothing after inversion
! if h1, filling matrix BDERI before inversion
!=====================================================================
                  if(PRESENT(ismooth)) then
                     bderi(id,id)=bderi(id,id)+dev
                     bderi(id,jj)=bderi(id,jj)-dev
                     h1(id) = h1(id) + deltarho*dev
                  else
                     dsum=dsum+deltarho*deltarho*dev
                     dcont=dcont+1
                  end if
               end do
            end if
         end do
      end if

    END SUBROUTINE DSMOOTH
