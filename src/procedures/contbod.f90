!=====================================================================
!=====================================================================
!   SUBROUTINE CONTBOD
!=====================================================================
!=====================================================================
!> This subroutine looks for contiguous bodies, in order to store
!! the parameters of the smoothing matrix Cs in ISIDE and JSIDE
!! It is called if ismooth = 1 and for the first iteration (iiter = 1)
!
! 05/02/2001 I fix something weird: it needed blocnod file, but did
! not use the information in it... The initial subroutine wrote by M.J.
! was even more weird...
!
! Called by MAIN
! calls none
!=====================================================================

    SUBROUTINE CONTBOD(xb,yb,nbod,ilay,iside,jside)

      USE MOD_unit
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of SMODEL
!=====================================================================
      integer,intent(in)                     :: nbod
      integer,DIMENSION(:),intent(in)        :: ilay
      integer,DIMENSION(:),intent(inout)     :: iside
      integer,DIMENSION(:,:),intent(inout)   :: jside

      real(kind=8),DIMENSION(:,:),intent(in) :: xb,yb
!=====================================================================
! Declaration of the dummy arguments of SMODEL
!=====================================================================
      integer                     :: i,j,k,m,nunode,vpn,ier
!=====================================================================
! Initialization of the arrays ISIDE,JSIDE
!=====================================================================
      iside(:) = 0
      jside(:,:) = 0

      write(inout,*)'    Call of CONTBOD subroutine'
      write(inout,*)'    Looking for contiguous density bodies...'
      write(inout,*)'    Store it in arrays ISIDE,JSIDE'
!=====================================================================
! Loop over the density bodies (i)
!=====================================================================
      do i = 1,nbod

!=====================================================================
! Look for the contiguous blocks in the following ones
!=====================================================================
         do j=i+1,nbod
            if(ilay(j).eq.ilay(i)) then
               if(xb(j,1).eq.xb(i,2) .and. yb(j,1).eq.yb(i,1)) then
                        
                  iside(i) = iside(i) + 1
                  jside(i,iside(i)) = j

                  iside(j) = iside(j) + 1
                  jside(j,iside(j)) = i

               else if(yb(j,1).eq.yb(i,2) .and. xb(j,1).eq.xb(i,1)) then

                  iside(i) = iside(i) + 1
                  jside(i,iside(i)) = j

                  iside(j) = iside(j) + 1
                  jside(j,iside(j)) = i

               end if
            end if
         end do

      end do


    END SUBROUTINE CONTBOD
