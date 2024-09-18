!=====================================================================
!=====================================================================
!
!			SUBROUTINE ZEROAV
!
!> subroutine reduces the mean of velocity and density variations for
!! each layer to 0.
!!
!! PAR array will then have a zero average
!!
!! Originally, only density blocks or velocity nodes hitten by rays are 
!! taken into account. The average is calculated only about them.
!! But I think it's rubbish, since all the density blocks are inverted!
!! We must take only the constrained velocity nodes, as information is
!! contained only in them, but we must take ALL the density blocks, as
!! they all account for gravity signal.
!!
!! First, compute the average density and then remove the average
!! from the density parameters
!! 
!! Same for velocity
! Calls: none
! called by MAIN
!
!=====================================================================
!=====================================================================
    SUBROUTINE ZEROAV(x_block,y_block,par,drmax,ibove)

      USE MOD_unit
      USE MOD_layer
      USE MOD_delim
      USE MOD_inv

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of ZEROAV
!=====================================================================
      integer,DIMENSION(:),intent(in)         :: x_block,y_block
      integer,DIMENSION(:),intent(in)         :: ibove

      real(kind=8),DIMENSION(:),intent(in)    :: drmax
      real(kind=8),DIMENSION(:),intent(inout) :: par
        
!=====================================================================
! Declaration of the dummy arguments of ZEROAV
!=====================================================================
      integer                               :: id,id2,ncount,ncount2
      integer                               :: ier,k,j,i,m,nunode,l
      integer                               :: ndens,nvelo,vpn

      real(kind=8)                          :: fac
      real(kind=8),DIMENSION(:),ALLOCATABLE :: averd,averv
!=====================================================================
! BC = Density Block counter
! NC = Velocity nodes counter
! ID = delimiter for the density parameters
! IV = delimiter for the velocity parameters
!=====================================================================
      id = ibegin(1)-1
      id2= ibegin(1)-1
      ncount = 0
      ncount2= 0
      
      ALLOCATE (averd(nlayer))
      ALLOCATE (averv(nznode-1))
      averd(:) = 0.d0
      averv(:) = 0.d0

      write(*,*)''
      write(*,*)'PUT THE PARAMETER AVERAGE TO 0'
      write(inout,*)''
      write(inout,*)'PUT THE PARAMETER AVERAGE TO 0'

!=====================================================================
! ndens = number of density parameters accounted for the average
! nvelo = number of velocity parameters accounted for the average
! averd = average of density parameter for each layer
! averv = average of velocity parameter for each layer
!=====================================================================

!=====================================================================
! First, compute the average density
! All blocks have to be taken into account
!=====================================================================
      if(INVD.or.INVGR) then
         do k=1,nlayer

	    ndens = 0

            if(ibegin(1).eq.0) then
               write(*,*)'Error in consistency in ZEROAV: INVD=TRUE but &
                    &ibegin(1)=0'
               STOP 'in ZEROAV!'
            end if

            do j=1,y_block(k)
               do i=1,x_block(k)
                  id = id+1
                  averd(k) = averd(k) + par(id)
                  ndens = ndens+1
               end do
            end do

            if(ndens.gt.0) then
               averd(k) = averd(k)/ndens
               write(*,'(5x,''Average density for layer '',i3,'': &
                    &'',f10.3)') k,averd(k)
               write(inout,'(5x,''Average density for layer '',i3,'': &
                    &'',f10.3)') k,averd(k)
            else
               write(*,*)''
               write(*,*)'Inconsistency in ZEROAV: INVD=.true. but ndens=0'
               STOP 'in ZEROAV!'
            end if
!=====================================================================
! Remove the average from the density parameters
!=====================================================================

            do j=1,y_block(k)
               do i=1,x_block(k)
                  id2 = id2+1
                  if(dabs(par(id2)-averd(k)).le.drmax(k)) then
                     par(id2) = par(id2) - averd(k)
                  else
                     if(averd(k).ne.0d0) then
                        fac = (drmax(k)-dabs(par(id2)))/averd(k)
                        par(id2) = par(id2) - averd(k)*fac
                     end if
                  end if
               end do
            end do
         end do
      end if
!=====================================================================
! Second, compute the average velocity
! Only constrained blocks are taken into account
!=====================================================================
      if(INVV) then

         do k=1,nznode-1

            nvelo=0

            if(ibegin(2).eq.0) then
               write(*,*)''
               write(*,*)'Inconsistency in ZEROAV. INVV=.true. but &
                    &ibegin(2)=0'
               STOP 'in ZEROAV!'
            end if
           
            do j=1,nynode
               do i=1,nxnode
                  ncount = ncount+1
                  if(ibove(ncount).ne.0) then
                     averv(k) = averv(k) + par(ibove(ncount))
                     nvelo = nvelo+1
                  end if
               end do
            end do

            if(nvelo.gt.0) then
               averv(k) = averv(k)/nvelo
               write(*,'(5x,''Average perturbation of velocity for &
                    &layer '',i3,'': '',f10.3)') k,averv(k)
               write(inout,'(5x,''Average perturbation of velocity for &
                    &layer '',i3,'': '',f10.3)') k,averv(k)
            else
               write(*,*)''
               write(*,*)'Inconsistency in ZEROAV: INVV=.true. but nvelo=0'
               STOP 'in ZEROAV!'
            end if
!=====================================================================
! Remove average from the velocity parameters
!=====================================================================
            do j=1,nynode
               do i=1,nxnode
                  ncount2=ncount2+1
                  if(ibove(ncount2).ne.0) then
                     par(ibove(ncount2)) = par(ibove(ncount2)) - averv(k)
                  end if
               end do
            end do

         end do
      end if

      DEALLOCATE(averd)
      DEALLOCATE(averv)

    END SUBROUTINE ZEROAV
