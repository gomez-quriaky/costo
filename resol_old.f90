!=====================================================================
!=====================================================================
!
!			SUBROUTINE RESOL
!
! This subroutine calculates the resolution terms for BDERI matrix
! The formula given in Zeyen & Achauer (1997) seems to be wrong for me
! as the dimensions of matrices are not respected in this way
! (equation (6), p.158 of their paper).
! Instead, I take inspiration from ACH method resolution, and it seems
! to better fit.
!
! Calls : DGEMM
! called by MAIN
!
!=====================================================================
!=====================================================================

    SUBROUTINE RESOL(bderi,cderi,npar,xb,yb,ilay,vxnodes,vynodes,ibove)

      USE MOD_unit
      USE MOD_layer
      USE MOD_inv
      USE MOD_delim

      USE MOD_dgemm

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of RESOL
!=====================================================================
      integer,intent(in)                       :: npar
      integer,DIMENSION(:),intent(in)          :: ilay,ibove

      real(kind=8),DIMENSION(:),intent(in)     :: vxnodes,vynodes
      real(kind=8),DIMENSION(:,:),intent(in)   :: bderi,cderi,xb,yb
!=====================================================================
! Declaration of the dummy arguments of RESOL
!=====================================================================
      integer                             :: ier,i,j,k,ibloc,inode

      real(kind=8)                        :: x,y
      real(kind=8),DIMENSION(npar,npar)   :: matres
      real(kind=8),PARAMETER              :: one = 1.0d0
      real(kind=8),PARAMETER              :: zero = 0.0d0
!=====================================================================
! Initialization of some variables
!=====================================================================
      matres(:,:) = 0.0
      ibloc = 0
      inode = 0
!=====================================================================
! Calculation of the resolution matrix, but is it correct ?
! In Zeyen & Achauer (97) it is: 
! (At.Cv**(-1).A + Cp**(-1) - Cb**(-1).Db - Cs**(-1).Ds)**(-1).A
! which is wrong, even for the dimensions of matrices!!
! Instead it seems safer to take: 
! (At.Cv**(-1).A + Cp**(-1) - Cb**(-1).Db - Cs**(-1).Ds)**(-1).(At.Cv**(-1).A)
!=====================================================================
!      matres(:,:) = MATMUL(bderi,cderi)
      CALL DGEMM('N','N',npar,npar,npar,one,bderi,npar,cderi,npar,zero,matres,npar)
!=====================================================================
! Write the results in the output files RESOL.VEL
!=====================================================================
      if(INVV) then

         open(resolvel,file='resol.vel',status='UNKNOWN',iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in opening file: RESOL.RES, logical unit ', &
                 resolvel
            STOP 'in RESOL!'
         end if

         write(*,*)'    Diagonal elts for velocity resolution stored &
              &in RESOL.VEL file'
         write(inout,*)'    Diagonal elts for velocity resolution stored &
              &in RESOL.VEL file'
         write(resolvel,*)'Diagonal elements for the resolution matrix'
         write(resolvel,*)'VELOCITY PARAMETERS (x,y,layer,res) NI=Non Inverted'

         do k=1,nznode-1
            do j=1,nynode
               do i=1,nxnode
                  inode=inode+1
                  if (ibove(inode).ne.0) then
                     write(resolvel,'(f10.3,2x,f10.3,2x,i3,2x,f10.6)') vxnodes(i),&
                          vynodes(j),k,matres(ibove(inode),ibove(inode))
                  else
                     write(resolvel,'(f10.3,2x,f10.3,2x,i5,''  NI'')') vxnodes(i),&
                          vynodes(j),k
                  end if
               end do
            end do
         end do

      end if
!=====================================================================
! Write the results in the output files RESOL.DEN
!=====================================================================
      if(INVD) then

         open(resolden,file='resol.den',status='UNKNOWN',iostat=ier)
         if(ier.ne.0) then
            write(*,*)'Error in opening file: RESOL.DEN, logical unit ', &
                 resolden
            STOP 'in RESOL!'
         end if

         write(*,*)'    Diagonal elts for density resolution stored &
              &in RESOL.DEN file'
         write(inout,*)'    Diagonal elts for density resolution stored &
              &in RESOL.DEN file'
         write(resolden,*)'Diagonal elements for the resolution matrix'
         write(resolden,*)'DENSITY PARAMETERS (x,y,layer,res)'

         do i=ibegin(1),iend(1)
            ibloc=ibloc+1
            x=(xb(ibloc,1)+xb(ibloc,2))*0.5
            y=(yb(ibloc,1)+yb(ibloc,2))*0.5
            write(resolden,'(f10.3,2x,f10.3,2x,i5,2x,f10.6)') x,y,ilay(ibloc),&
                 matres(i,i)
         end do

      end if
!=====================================================================
! Closing the files RESOL.DEN and RESOL.VEL
!=====================================================================
      close(resolvel,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file: RESOL.VEL, logical unit ', &
              resolvel
         STOP 'in RESOL!'
      end if
      close(resolden,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file: RESOL.DEN, logical unit ', &
              resolden
         STOP 'in RESOL!'
      end if


    END SUBROUTINE RESOL
