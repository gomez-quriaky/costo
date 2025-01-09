!=====================================================================
!=====================================================================
!
!			SUBROUTINE FINDNOD
!
!> subroutine creates file blocnod where for each density block
!! the corresponding number of the velocity nodes including
!! information are stored.
!!
!! It should be used only when velocity AND density inversion is 
!! required. Otherwise, it is nonesense to use it.
!
! Calls: none
! called by MAIN used in DENVEL,SMODEL,ZEROAV ..
!
!=====================================================================
!=====================================================================
    SUBROUTINE FINDNOD(xb,yb,zb,vxnodes,vynodes,vznodes,nbod,ibove)

      USE MOD_unit
      USE MOD_layer

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of FINDNOD
!=====================================================================
      integer,intent(in)                     :: nbod
      integer,DIMENSION(:),intent(in)        :: ibove

      real(kind=8),DIMENSION(:,:),intent(in) :: xb,yb,zb
      real(kind=8),DIMENSION(:),intent(in)   :: vxnodes,vynodes,vznodes
!=====================================================================
! Declaration of the dummy arguments of FINDNOD
!=====================================================================
      integer                                :: ier,m,node,nunode
      integer                                :: i,j,k
      integer,DIMENSION(:),ALLOCATABLE       :: HIT
!=====================================================================
! Open the file with the corresponding blocks(density)/nodes(velocity)
!
! This file is organized by blocks with in the header 
!      - the number of the density block constrained according to dens.mod
!      - the number of the density block (incremented by 1)
! and for each block
!      - the numbers of nodes contained in it
!=====================================================================
      open(blocnod,file='bloc.nod',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file BLOC.NOD, logical unit ',&
              blocnod,'. Stooooop in FINDNOD!'
         STOP
      end if
!=====================================================================
! Allocate maximum memory for HIT array
!=====================================================================
      ALLOCATE(HIT((nznode-1)*nxnode*nynode))
      HIT(:) = 0
!=====================================================================
! Loop over the density blocks
!=====================================================================
      do m=1,nbod
         node=0
         nunode=0
!=====================================================================
! Loop over the velocity nodes
!=====================================================================
         do k=1,nznode-1
            do j=1,nynode
               do i=1,nxnode
!=====================================================================
! check if the node is constrained first (with ibove array)
!=====================================================================
                  node=node+1
                  if(ibove(node).ne.0) then
!=====================================================================
! check if the block contains the node
! if yes, store it in HIT array
! nunode = number of nodes that constrain the density block
!=====================================================================
                     if(xb(m,1).le.vxnodes(i).and.xb(m,2).ge.vxnodes(i)&
                        .and.yb(m,1).le.vynodes(j).and.yb(m,2).ge.vynodes(j)&
                        .and.(zb(m,1)).le.vznodes(k).and.(zb(m,2)).ge.vznodes(k))&
                        then
                        nunode=nunode+1
                        HIT(nunode)=node
                     end if
                  end if
               end do
            end do
         end do

         write(blocnod,'(i5,i5)') m,nunode
         do i=1,nunode
            write(blocnod,'(i5)') HIT(i)
         end do

      end do

      DEALLOCATE(HIT)
!=====================================================================
! Close the file with the corresponding blocks(density)/nodes(velocity)
!=====================================================================
      close(blocnod,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file BLOC.NOD, logical unit ',&
              blocnod,'. Stooooop in FINDNOD!'
         STOP
      end if
    END SUBROUTINE FINDNOD
