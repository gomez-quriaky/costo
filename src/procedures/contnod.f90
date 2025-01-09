!=====================================================================
!=====================================================================
!   SUBROUTINE CONTNOD
!=====================================================================
!=====================================================================
!
!> This subroutine looks for contiguous velocity nodes, in order to store
!! the parameters of the smoothing matrix Cs in IVSIDE and JVSIDE
!
! Called by MAIN
! calls none
!=====================================================================

    SUBROUTINE CONTNOD(ibove,ivside,jvside)

      USE MOD_unit
      USE MOD_layer
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of SMODEL
!=====================================================================
      integer,DIMENSION(:),intent(in)        :: ibove
      integer,DIMENSION(:),intent(inout)     :: ivside
      integer,DIMENSION(:,:),intent(inout)   :: jvside
!=====================================================================
! Declaration of the dummy arguments of SMODEL
!=====================================================================
      integer                     :: i,j,k,ii,jj,kk,nunode,folnod
!=====================================================================
! Initialization of the arrays ISIDE,JSIDE
!=====================================================================
      ivside(:) = 0
      jvside(:,:) = 0

      write(inout,*)''
      write(inout,*)'    Call of CONTNOD subroutine'
      write(inout,*)'    Looking for contiguous velocity nodes...'
      write(inout,*)'    Store it in arrays IVSIDE,JVSIDE'

!=====================================================================
! Loop over the velocity nodes 
!=====================================================================
      nunode = 0
      do k=1,nznode-1
         do j=1,nynode
            do i=1,nxnode
               nunode=nunode+1
!=====================================================================
! Look for a constrained velocity node (ibove(nunode).ne.0)
!=====================================================================
               if(ibove(nunode).ne.0) then
!=====================================================================
! Loop over the following nodes to find the neighbouring ones
!=====================================================================
                  folnod=0
		  do kk=1,nznode-1
                     do jj=1,nynode
                        do ii=1,nxnode
                           folnod=folnod+1
!=====================================================================
! Choose a following constrained node
!=====================================================================
                           if(ibove(folnod).ne.0) then
!=====================================================================
! Choose the right layer and the contiguous node in y and x directions
!=====================================================================
                              if((kk.eq.k).and.&
                                (jj.eq.(j+1).or.jj.eq.(j-1).or.jj.eq.j).and.&
                                (ii.eq.(i+1).or.ii.eq.(i-1).or.ii.eq.i)) then
!=====================================================================
! Eliminate the originale node from the first loop
!=====================================================================
                                 if((i.eq.ii .and. j.ne.jj).or.&
                                   (i.ne.ii.and.j.eq.jj)) then

                                    ivside(nunode)=ivside(nunode)+1
                                    jvside(nunode,ivside(nunode)) = folnod

                                 end if
                              end if
                           end if
!=====================================================================
! End of loop over the following nodes
!=====================================================================
                        end do
                     end do
                  end do
                  
               end if
!=====================================================================
! End of loop over the nodes
!=====================================================================
            end do
         end do
      end do


    END SUBROUTINE CONTNOD
