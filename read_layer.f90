!=====================================================================
!=====================================================================
!
!			SUBROUTINE READ_LAYER
!
!> Read the maximum number of layer, x, y, blocks or nodes
!! in files dens.mod, vel.mod
!
! Calls: none
! Called by MAIN
!=====================================================================
!=====================================================================
    SUBROUTINE READ_LAYER

      USE MOD_unit
      USE MOD_layer
      USE MOD_inv

      IMPLICIT NONE
!=====================================================================
! Declaration of the dummy arguments of READ_LAYER
!=====================================================================
        
      integer                :: ier
      character(len=80)      :: dummy

!=====================================================================
! Open,read and close the density file
!=====================================================================

      if(INVD.or.INVGR) then
         open(indens,file=densmod,status='OLD',iostat=ier)
         if(ier.ne.0) then
            write(*,*) 'Error in opening file ',densmod
            write(*,*) 'Stooooop in READ_LAYER!'
            stop
         end if

         read(indens,'(a80)') dummy
         read(indens,*) nxbloc,nybloc

         close(indens)
      else
         nxbloc = 0
         nybloc = 0
      end if

!=====================================================================
! Open,read and close the velocity file
!=====================================================================
      if(INVV) then
         open(invelm,file=velmod,status='OLD',iostat=ier)
         if(ier.ne.0) then
            write(*,*) 'Error in opening file ',velmod
            write(*,*) 'Stooooop in READ_LAYER!'
            stop
         end if

         read(invelm,'(a80)') dummy
         read(invelm,*) nxnode,nynode,nznode
         
         if(nznode.ne.(nlayer+1)) then
            write(*,*)''
            write(*,*)'Inconsistency between PARAMETER.INP and ',velmod
            write(*,*)'nznode in ',velmod,' should be nlayer+1 (i.e. ',&
                 (nlayer+1),'). Stooooop in READ_LAYER!'
            STOP
         end if

         close(invelm)
      else
         nxnode = 0
         nynode = 0
         nznode = 0
      end if


    END SUBROUTINE READ_LAYER

