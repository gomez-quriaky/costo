!=====================================================================
!=====================================================================
!   SUBROUTINE BCALC
!
! This subroutine computes the B-coefficient of the linear relation
! between density and velocity.
! The mean coefficient is calculated from the densities and velocities
! only constrained by 3D raytracing.
! A coefficient of correlation is calculated between velocity and
! density to indicate if the linear relation is far or not.
! (i.e. a good correlation mean a good linear relationship between
! both parameters, a bad one means that there is a big dispersion
! between the parameters).
!
!=====================================================================
!=====================================================================
! 
! Called by MAIN
! calls none
!=====================================================================

    SUBROUTINE BCALC(nbod,nnod,ilay,par,dvrho,ibove)

      USE MOD_delim
      USE MOD_layer
      USE MOD_unit
      USE MOD_inv
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of BCALC
!=====================================================================
      integer,intent(in)                      :: nbod,nnod
      integer,DIMENSION(:),intent(in)         :: ilay,ibove

      real(kind=8),DIMENSION(:),intent(inout) :: par,dvrho
!=====================================================================
! Declaration of the dummy arguments of BCALC
!=====================================================================
      logical                                :: constrain

      integer                                :: ier,m,nunode,node,const_nod
      integer                                :: i,j,k,vpn,ii,nl,ib,ibloc
      integer                                :: parpos
      integer,DIMENSION(:),ALLOCATABLE       :: nrho,nvel,hitbloc

      real(kind=8)                           :: bcoef,corre,parho,parvel
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: sumrho,sumrho2,sumv,sumv2
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: sumvr
!=====================================================================
! Initializatio of some variables and arrays
!=====================================================================

      ALLOCATE (sumrho(nlayer))
      ALLOCATE (sumrho2(nlayer))
      ALLOCATE (sumv(nlayer))
      ALLOCATE (sumv2(nlayer))
      ALLOCATE (sumvr(nlayer))
      ALLOCATE (nrho(nlayer))
      ALLOCATE (nvel(nlayer))
      ALLOCATE (hitbloc(nnod))

      sumrho(:)=0.d0
      sumrho2(:)=0.d0
      sumv(:)=0.d0
      sumv2(:)=0.d0
      sumvr(:)=0.d0
      nrho(:)=0
      nvel(:)=0

      write(*,*)''
      write(*,*)'FORWARD CALCULATION OF B-COEFFICIENT'
      write(inout,*)''
      write(inout,*)'FORWARD CALCULATION OF B-COEFFICIENT'
!=====================================================================
! Opening the blonod file (BLOC.NOD) for constrained blocks
!=====================================================================
      open(blocnod,file='bloc.nod',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file BLOC.NOD, logical unit ',&
                   blocnod,'.Stooooop in BCALC!'
         STOP 'in BCALC!'
      end if
!=====================================================================
! First getting through the density bodies
!=====================================================================
      do i=ibegin(1),iend(1)

         read(blocnod,*) m,nunode
         ii=i-ibegin(1)+1
         nl=ilay(ii)
         if(ii.ne.m) then
            write(*,*)'Inconsistency in BCALC. Stooooop!'
            close(blocnod)
            STOP
         end if

         constrain = .FALSE.

         if(nunode.ne.0) then
            do j=1,nunode
               read(blocnod,*) vpn
               if(ibove(vpn).ne.0) then
                  constrain = .TRUE.
               end if
            end do
         end if
         
         if(constrain) then
            sumrho(nl) = sumrho(nl) + par(i)
            sumrho2(nl) = sumrho2(nl) + par(i)*par(i)
            nrho(nl) = nrho(nl) + 1
         end if

      end do
!=====================================================================
! Second getting through the velocity nodes
!=====================================================================
      node = 0

      do k=1,nlayer
         do j=1,nynode
            do i=1,nxnode

               rewind(blocnod)
               node=node+1
               const_nod = 0
!=====================================================================
! Read the file blocnod
!=====================================================================
               do ibloc=1,nbod
                  read(blocnod,*) m,nunode
                  if(m.ne.ibloc) then
                     write(*,*)'Inconsistency in BCALC. Stooooop!'
                     close(blocnod)
                     STOP
                  end if

                  if(nunode.ne.0) then
                     do ii=1,nunode
                        read(blocnod,*) vpn
                        if(vpn.eq.node) then
                           const_nod=const_nod+1
                           if(const_nod.gt.nnod) then
                              write(*,*)'Out of range subscript for &
                                       &array HITBLOC. Stooooop in BCALC!'
                              STOP
                           end if
                           hitbloc(const_nod) = m
                        end if
                     end do
                  end if
               end do
!=====================================================================
! End of reading the file blocnod
!
! Then, calculate the mean density (parho) associated to each node
! parpos=indice of the constrained density parameter
!=====================================================================
               if(const_nod.ge.1) then
                  parho=0.d0
                  do ii=1,const_nod
                     parpos=ibegin(1)-1+hitbloc(const_nod)
                     parho=parho+par(parpos)
                  end do
                  parho=parho/const_nod
               end if

               if(ibove(node).ne.0) then
                  parvel=par(ibove(node))
                  sumv(k)=sumv(k)+parvel
                  sumv2(k)=sumv2(k)+parvel*parvel
                  sumvr(k)=sumvr(k)+parvel*parho
                  nvel(k)=nvel(k)+1
               end if

            end do
         end do
      end do
!=====================================================================
! Closing the file blocnod
!=====================================================================
      close(blocnod,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file BLOC.NOD, logical unit ',&
                   blocnod,'.Stooooop in BCALC!'
         STOP
      end if

!=====================================================================
! Testing calculated B-coeff.
!
! Mean B-coeff = (mean V.rho)/(mean rho**2)
! Correlation coefficient is the Pearson's one:
!                  Sxy - Sx.Sy/N
! Corr = --------------------------------
!         _______________________________
!        V (Sx2 - (Sx)2/N)(Sy2 - (Sy)2/N)
!
! If calculated B-coeff do not coincide, the corresponding parameter
! is changed to the real one
!=====================================================================
      write(*,*)''
      write(*,*)'VELOCITY-DENSITY RELATION, B-COEF FROM BCALC'
      write(*,*)''
      write(inout,*)''
      write(inout,*)'VELOCITY-DENSITY RELATION, B-COEF FROM BCALC'
      write(inout,*)''
      do nl=1,nlayer
         ib=ibegin(3)-1+nl

         if(sumrho2(nl).gt.0.00001 .and. nvel(nl).ne.0 .and. sumv(nl).ne.0) then
            bcoef=(sumvr(nl)/nvel(nl))/(sumrho2(nl)/nrho(nl))
            corre=(sumvr(nl)-sumv(nl)*sumrho(nl)/nvel(nl))&
                 /dsqrt(DABS((sumv2(nl)-sumv(nl)*sumv(nl)/nvel(nl)))*&
                 DABS((sumrho2(nl)-sumrho(nl)*sumrho(nl))/nrho(nl)))
         else
            bcoef=0.d0
            corre=0.d0
         end if

         if(corre.ne.0.) then
            dvrho(nl)=bcoef
            par(ib)=bcoef
            write(inout,*)'LAYER ',nl
            write(inout,*)'    Really: ',bcoef,'*RHO. Corr. Coef.: ',corre
            write(inout,*)'    Number of points: ',nrho(nl)
         else
            write(inout,*)'LAYER ',nl
            write(inout,*)'    Velocity-density relation not defined'
            write(inout,*)'    Corr. Coef. = 0.d0'
         end if

      end do
      
      DEALLOCATE (sumrho)
      DEALLOCATE (sumrho2)
      DEALLOCATE (sumv)
      DEALLOCATE (sumv2)
      DEALLOCATE (sumvr)
      DEALLOCATE (nrho)
      DEALLOCATE (nvel)
      DEALLOCATE (hitbloc)

    END SUBROUTINE BCALC
