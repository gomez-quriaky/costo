!=====================================================================
!=====================================================================
!       SUBROUTINE VSMOOTH
!
!> This subroutine computes the smoothing constrained for BDERI matrix
!! for the velocity bodies, term Cs**(-1).Ds in BLDMAT
!
! Calls none
! Called by MAIN, VDIFFER
!
!=====================================================================
!=====================================================================
    SUBROUTINE VSMOOTH(smooth,ivside,jvside,par,ibove,vxnodes,vynodes,&
                         vsum,vcont,bderi,h1,ismooth)
        
      USE MOD_delim
      USE MOD_layer

      IMPLICIT NONE

!=====================================================================
! Declaration of in/out arguments of VSMOOTH
!=====================================================================
      integer,intent(out)                               :: vcont
      integer,OPTIONAL,intent(in)                       :: ismooth
      integer,DIMENSION(:),intent(in)                   :: ivside,ibove
      integer,DIMENSION(:,:),intent(in)                 :: jvside

      real(kind=8),intent(out)                          :: vsum
      real(kind=8),DIMENSION(:),intent(in)              :: smooth,par
      real(kind=8),DIMENSION(:),intent(inout)           :: h1
      real(kind=8),DIMENSION(:),intent(in)              :: vxnodes,vynodes
      real(kind=8),DIMENSION(:,:),intent(inout)         :: bderi
!=====================================================================
! Declaration of dummy arguments of VSMOOTH
!=====================================================================
      integer                 :: i,j,k,inod,ipar,j1,js,jj
      integer                 :: ik,iy,ix,n

      real(kind=8)            :: xnod1,ynod1,xnod2,ynod2
      real(kind=8)            :: dx,dy,dev,fs,deltav
!=====================================================================
! Initialization of some variables
! fs is the diagonal element of the Cs**(-1) matrix
!=====================================================================
      vcont=0
      vsum=0
      fs=1./(smooth(2)*smooth(2))
!=====================================================================
! Loop over velocity parameters (nodes)
!=====================================================================
      inod=0
      do k=1,nznode-1
         do j=1,nynode
            do i=1,nxnode

               inod=inod+1
!=====================================================================
! If the node is constrained  AND has contiguous nodes,
! ipar=number of the corresponding parameter in PAR array
! jj=number of the parameter in PAR array that is contiguous to ipar
!=====================================================================
               if((ibove(inod).ne.0) .and. (ivside(inod).ne.0)) then

                  ipar=ibove(inod)
                  xnod1=vxnodes(i)
                  ynod1=vynodes(j)

                  do j1=1,ivside(inod)
                     js=jvside(inod,j1)
                     jj=ibove(js)
!=====================================================================
! If the adjacent node is constrained (min. rays pass through it)
! it compute the model roughness 
!=====================================================================
                     if(jj.ne.0) then

                        n=0
                        do ik=1,nznode-1
                           do iy=1,nynode
                              do ix=1,nxnode

                                 n=n+1
                                 if(n.eq.js) then
                                    xnod2=vxnodes(ix)
                                    ynod2=vynodes(iy)
                                    dx=xnod2-xnod1
                                    dy=ynod2-ynod1
                                    dev=fs/(dx*dx+dy*dy)
                                    deltav=par(jj)-par(ipar)

                                    if(PRESENT(ismooth)) then
                                       bderi(ipar,ipar)=bderi(ipar,ipar)+dev
                                       bderi(ipar,jj)=bderi(ipar,jj)-dev
                                       h1(ipar)=h1(ipar)+deltav*dev
                                    else
                                       vsum=vsum+deltav*deltav*dev
                                       vcont=vcont+1
                                    end if

                                 end if

                              end do
                           end do
                        end do

                     end if
                  end do

               end if

            end do
         end do
      end do



    END SUBROUTINE VSMOOTH
