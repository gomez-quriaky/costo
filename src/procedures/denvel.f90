!=====================================================================
!=====================================================================
!
!			SUBROUTINE DENVEL
!
!> Subroutine prepares the linear set of equations for
!! connection of density and velocity by a linear correlation factor
!
! Calls none
! called by BLDMAT
!
!=====================================================================
!=====================================================================

    SUBROUTINE DENVEL(bderi,h1,ilay,ibove,par,ddvr,nbod)
        
      USE MOD_delim
      USE MOD_unit


      use ISO_Fortran_env, only: i4=>int32, i8=>int64
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of DENVEL
!=====================================================================
      integer(i8),intent(in)                         :: nbod
      integer(i8),DIMENSION(:),intent(in)            :: ilay,ibove
      
      real(kind=8),DIMENSION(:),intent(in)       :: par,ddvr
      real(kind=8),DIMENSION(:),intent(inout)    :: h1
      real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
!=====================================================================
! Declaration of the dummy arguments of DENVEL
!=====================================================================
      integer                                    :: ier,id,iv,ib,ibod
      integer                                    :: m,nunode,inod,vpn
      integer                                    :: nl,vpar

      real(kind=8)                               :: cov,diff1,diff2,coef

!=====================================================================
! Opening of the blocnod file for constrained blocks
!=====================================================================
      open(blocnod,file='bloc.nod',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file BLOC.NOD, logical unit ',&
                   blocnod,'. Stooooop in DENVEL!'
         STOP
      end if
!=====================================================================
! Define some variables for starting values
! As the subroutine is done only if INVD.and.INVV, ibegin(1,2) .ne. 0
!=====================================================================
      id=ibegin(1)-1
!=====================================================================
! Loop over density bodies
!=====================================================================
      do ibod=1,nbod

         id=id+1
         nl=ilay(ibod)
         ib=ibegin(3)+nl-1
         read(blocnod,*) m,nunode
         if(m.ne.ibod) then
            write(*,*)'Inconsistency in DENVEL subroutine. Stooooop!'
            STOP
         end if

         if(nunode.ne.0) then
!=====================================================================
! Loop over the velocity nodes constraining the density block
!=====================================================================
            do inod=1,nunode
               read(blocnod,*) vpn
!=====================================================================
! If the node is constrained (i.e. rays pass through it)
! cov = covariance?
! diff1 = difference between inverted velocity and calculated one with
!         actual B-coef (par(ib))
!=====================================================================
               if(ibove(vpn).ne.0) then
                  iv=ibove(vpn)

                  cov=1./(ddvr(nl)*ddvr(nl))
                  diff1=(par(iv) - par(ib)*par(id))
                  diff2=(diff1-par(ib)*par(id))*cov

                  diff1=diff1*cov
                  coef=par(ib)*cov

                  bderi(id,id) = bderi(id,id)+par(ib)*coef
                  bderi(id,iv) = bderi(id,iv) - coef
                  bderi(iv,id) = bderi(iv,id) - coef
                  bderi(iv,iv) = bderi(iv,iv) + cov
                  h1(id) = h1(id)+par(ib)*diff1

                  bderi(ib,ib) = bderi(ib,ib) + cov*par(id)*par(id)
                  bderi(ib,id) = bderi(ib,id) - diff2
                  bderi(id,ib) = bderi(id,ib) - diff2
                  bderi(ib,iv) = bderi(ib,iv) - cov*par(id)
                  bderi(iv,ib) = bderi(iv,ib) - cov*par(id)
                  h1(ib) = h1(ib) + diff1*par(id)
                  h1(iv) = h1(iv) - diff1

               end if
            end do

         end if

      end do

!=====================================================================
! Closing of the blocnod file
!=====================================================================
      close(blocnod,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file BLOC.NOD, logical unit ',&
                   blocnod,'. Stooooop in DENVEL!'
         STOP
      end if


    END SUBROUTINE DENVEL
