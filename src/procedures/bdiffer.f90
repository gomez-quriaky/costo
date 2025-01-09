!=====================================================================
!=====================================================================
!
!			SUBROUTINE BDIFFER
!
!> subroutine calculates the difference between observed and calculated
!! data (USEDATA, CALDATA) and uses it to compute the standard
!! deviation (SDEV)
!!
!! Also calculates the difference between previous (PAR0) and actual (PAR)
!! parameters
!
! Calls: BLOCKMEAN(VEL3(INTMAP))
! called by MAIN
!
!=====================================================================
!=====================================================================
   !INCLUDE 'INTERF/MOD_blockmean.f'

    SUBROUTINE BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,&
                       velco,v0,nbod,dtot,ibove)
        
      USE MOD_delim
      USE MOD_unit

      USE MOD_blockmean

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of BDIFFER
!=====================================================================
      integer,intent(in)                        :: iiter,nbod
      integer,DIMENSION(:),intent(in)           :: ilay,ibove

      real(kind=8),DIMENSION(:),intent(in)      :: par,par0,varpar
      real(kind=8),DIMENSION(:),intent(in)      :: velco,v0,varb
      real(kind=8),DIMENSION(:,:),intent(in)    :: xb,yb,zb
      real(kind=8),DIMENSION(:,:),intent(inout) :: dtot
!=====================================================================
! Declaration of the dummy arguments of BDIFFER
!=====================================================================
      logical                                :: constrain

      integer                                :: i,j,ier,ib,i1,conbod
      integer                                :: m,nunode,vpn

      real(kind=8)                           :: Vm,dev,difr,difb
!=====================================================================
! Writting some outputs
!=====================================================================
      write(*,*)''
      write(*,*)'COMPUTING THE DIFFERENCES BETWEEN OBS. AND CALC B-COEFF'


      if(iiter.eq.1) then
         write(inout,*)''
         write(inout,*)'BDIFFER CALLED FOR COMPUTING INITIAL SDT DEVIATION'
         write(inout,*)'ITERATION = 0'
      else
         write(inout,*)''
         write(inout,*)'BDIFFER CALLED FOR DIFFERENCE BETWEEN OBS. and CALC'
         write(inout,*)'ITERATION = ',iiter
      end if

      if(ibegin(3).eq.0 .or. ibegin(1).eq.0 .or. ibegin(2).eq.0) then
         write(*,*)''
         write(*,*)'Inconsistency in BDIFFER: inversion of linear &
                   &relation is asked, but ibegin(i)=0. Stooooop!'
         STOP
      end if
!=====================================================================
! Differences between previous and actual parameters
!=====================================================================
      ib = 0
      dev = 0.d0
      difb = 0.d0
      do i=ibegin(3),iend(3)
         ib = ib+1
         dev = par(i) - par0(i)
         difb = difb + dev*dev*varpar(i)
      end do
      difb = difb/ib
!=====================================================================
! Differences between parameters with linear relation
!=====================================================================
      difr = 0.d0
      dev = 0.d0
      
!      i1=ibegin(1)-1
!      i3=ibegin(3)-1

      open(blocnod,file='bloc.nod',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)''
         write(*,*)'Error in opening the file BLOC.NOD, &
              &logical unit ', blocnod
         STOP 'in BDIFFER subroutine.'
      end if

!      do i=1,nbod
!         i1=i1+1
!         nl = i3+ilay(i)
      conbod = 0
      do i=ibegin(1),iend(1)
         read(blocnod,*) m,nunode
         i1=i-ibegin(1)+1
         ib=ibegin(3)-1+ilay(i1)
         if(i1.ne.m) then
            write(*,*)'Inconsistency in BDIFFER. Stooooop!'
            close(blocnod)
            STOP 'in BDIFFER'
         end if

         constrain = .FALSE.

         if(nunode.ne.0) then
            do j=1,nunode
               read(blocnod,*) vpn
               if(ibove(vpn).ne.0) then
                  constrain = .TRUE.
               end if
            end do
            if(constrain) then
               conbod = conbod + 1
!=====================================================================
! Compute the average velocity Vm in the density block from 25 velocity
! nodes from vel.mod file and VEL3 subroutine
!
! del = Vmean - (Bcoef * density)
!=====================================================================
               CALL BLOCKMEAN(xb(i1,1),xb(i1,2),yb(i1,1),yb(i1,2),&
                                 zb(i1,1),zb(i1,2),Vm,velco)
               Vm = Vm - v0(ilay(i1))
               dev = Vm - par(ib)*par(i)
               difr = difr + dev*dev/(varb(ilay(i1))*varb(ilay(i1)))
            end if
         end if
      end do
!      difr = difr/nbod
      if(conbod.ne.0) then
         difr = difr/conbod
      else
         write(*,*)'Number of constrained density blocks = 0'
      end if

      close(blocnod,iostat=ier)
      if(ier.ne.0) then
         write(*,*)''
         write(*,*)'Error in closing the file BLOC.NOD, &
              &logical unit ', blocnod
         STOP 'in BDIFFER subroutine.'
      end if
!=====================================================================
! It's time to write some outputs in inout logical unit
!=====================================================================
      dtot(iiter,1) = difr
      dtot(iiter,2) = difb
      dtot(iiter,3) = difb + difr

      write(inout,'(5x,''Sum of mean diff. for B-parameters:'',22(''.''),&
           &f12.3)') dtot(iiter,3)
      write(inout,*)'    whereof:'
      write(inout,'(13x,''Mean diff. between act. and prev. B-parameters:..&
           &'',f12.3)') difb
      write(inout,'(13x,''Mean deviation from linearity law:'',15(''.''),&
           &f12.3)') difr

      write(*,'(5x,''Sum of mean diff. for B-parameters:'',22(''.''),&
           &f12.3)') dtot(iiter,3)
      write(*,*)'    whereof:'
      write(*,'(13x,''Mean diff. between act. and prev. B-parameters:..&
           &'',f12.3)') difb
      write(*,'(13x,''Mean deviation from linearity law:'',15(''.''),&
           &f12.3)') difr


    END SUBROUTINE BDIFFER
           
!=====================================================================
!=====================================================================
!       SUBROUTINE BLOCKMEAN
!
!> This subroutine computes the average velocity for a density block,
!! from VEL3 subroutine and velocity stored in vel.mod file
!
! Calls VEL3(INTMAP)
! Called by BDIFFER
!
!=====================================================================
!=====================================================================
    SUBROUTINE BLOCKMEAN(xb1,xb2,yb1,yb2,zb1,zb2,Vm,velco)

      USE MOD_vel3

      IMPLICIT NONE
        
!=====================================================================
! Declaration of in/out arguments of BLOCKMEAN
!=====================================================================
      real(kind=8),intent(out)              :: Vm
      real(kind=8),intent(in)               :: xb1,xb2,yb1,yb2,zb1,zb2
      real(kind=8),DIMENSION(:),intent(in)  :: velco
!=====================================================================
! Declaration of dummy arguments of BLOCKMEAN
!=====================================================================
      integer                           :: i,ip,jp,kp

      real(kind=8)                      :: deltaz,zl1,xhalf,yhalf
      real(kind=8),DIMENSION(5,5)       :: v
!=====================================================================
! Divide the block into 5 layers
! Search 5 velocity points for these 5 sub-layers
!=====================================================================
      vm = 0.d0
      deltaz = DABS((zb2-zb1)/4.d0)
      zl1 = zb1
      xhalf = (xb1+xb2)/2.d0
      yhalf = (yb1+yb2)/2.d0

      do i=1,5
         CALL VEL3(xb1,yb1,zl1,v(1,i),velco,ip,jp,kp)
         CALL VEL3(xb1,yb2,zl1,v(2,i),velco,ip,jp,kp)
         CALL VEL3(xb2,yb1,zl1,v(3,i),velco,ip,jp,kp)
         CALL VEL3(xb2,yb2,zl1,v(4,i),velco,ip,jp,kp)
         CALL VEL3(xhalf,yhalf,zl1,v(5,i),velco,ip,jp,kp)

         zl1 = zb1+(i*deltaz)
      end do

      do i=1,5
         Vm = Vm+v(1,i)+v(2,i)+v(3,i)+v(4,i)+v(5,i)
      end do

      Vm = Vm/25.d0
      
    END SUBROUTINE BLOCKMEAN
