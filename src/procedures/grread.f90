!! ne pas oublier de changer le module MOD_grread si on change la subroutine !!!
!=====================================================================
!=====================================================================
! SUBROUTINE GRREAD
!=====================================================================
    SUBROUTINE GRREAD(fdata,ifil,modvar,rtvar,npts,grFX,grFY,grFZ,&
                        XX,XY,XZ,YY,YZ,ZZ,ptvar3,ncomp)

        USE MOD_unit
        USE MOD_layer
        USE MOD_size
        USE MOD_inv

      IMPLICIT NONE

      character(len=80),intent(in)          :: fdata

      integer,intent(in)                      :: ifil
      integer,DIMENSION(:),intent(inout)      :: npts
      integer,intent(inout)                   :: ncomp
      integer,DIMENSION(:),intent(in)         :: modvar
      real(kind=8),DIMENSION(:),intent(in)    :: rtvar
      real(kind=8),DIMENSION(:),pointer       :: ptvar3

      real(kind=8),DIMENSION(:),pointer       :: grFX,grFY,grFZ
      real(kind=8),DIMENSION(:),pointer       :: XX,XY,XZ,YY,YZ,ZZ


!=====================================================================
! Declaration of dummy arguments for GRREAD
!=====================================================================
      character(len=80)                        :: dummy

      integer                                  :: num,nptot,i,j,itake,iter,nbn
      real(kind=8),DIMENSION(:,:),ALLOCATABLE  :: datatot
      real(kind=8)                             :: xmin,xmax,ymin,ymax
      real(kind=8)                             :: perc
      real(kind=8),DIMENSION(:),ALLOCATABLE    :: ptvar

!      print*,"on est dans la subroutine GRREAD"   

!=====================================================================
! Read the general boundaries and the used ones in data file
!=====================================================================
      inquire(file=fdata,number=num)
      write(*,*)'    logical unit currently read:',num,',',fdata
      read(num,*) dummy
      read(num,*) nptot

!=====================================================================
! Allocation of the data and position array
!=====================================================================
      ALLOCATE (datatot(nptot,9))

!=====================================================================
! Read coordinates of the inversion area 
!=====================================================================   
      read(num,*) dummy
      read(num,*) xmin,xmax,ymin,ymax
      read(num,*) dummy

! calcul du nombre de composante à inverser (ou modéliser)

      ncomp=0
      nbn=0     ! pour prendre en compte le nombre de méthode que l'on inverse en plus de la gradio

      if (INVD) then
        nbn=nbn+1
      endif

      if (INVV) then
        nbn=nbn+1
      endif

      do i=1,6
        if (rtvar(nbn+i) .ne. 0) then
          ncomp=ncomp+1
        endif
      enddo

!=====================================================================
! Read position of all x,y,z and FTG data
!=====================================================================   
      do i=1,nptot
         read(num,*) datatot(i,1),datatot(i,2),datatot(i,3),datatot(i,4),datatot(i,5),&
            datatot(i,6),datatot(i,7),datatot(i,8),datatot(i,9)

         if(datatot(i,1).ge.xmin.and.datatot(i,1).le.xmax) then
            if(datatot(i,2).ge.ymin.and.datatot(i,2).le.ymax) then
                  npts(ifil)=npts(ifil)+1
            end if
         end if
      end do

!=====================================================================
! Allocation of the data and position array for used data
!=====================================================================
      write(*,*)'    Number of used data for FTG file: ',npts(ifil)
      write(inout,*)''
      write(inout,*)'    Number of used data for FTG file: ',&
           npts(ifil)
      ALLOCATE (grFX(npts(ifil)))
      ALLOCATE (grFY(npts(ifil)))
      ALLOCATE (grFZ(npts(ifil)))
      ALLOCATE (XX(npts(ifil)),XY(npts(ifil)),XZ(npts(ifil)))
      ALLOCATE (YY(npts(ifil)),YZ(npts(ifil)),ZZ(npts(ifil)))
      ALLOCATE (ptvar(npts(ifil)*6))
      ALLOCATE (ptvar3(npts(ifil)*ncomp))   ! MP36

!=============================================================================
! Storage of the used data in grFX, grFY, grFZ, XX, XY, XZ, YY, YZ, ZZ arrays
!=============================================================================

    itake=0
      do i=1,nptot
         if(datatot(i,1).ge.xmin.and.datatot(i,1).le.xmax) then
            if(datatot(i,2).ge.ymin.and.datatot(i,2).le.ymax) then

              itake=itake+1 

              grFX(itake)= datatot(i,1)
              grFY(itake)= datatot(i,2)
              grFZ(itake)= datatot(i,3)

 
! ajout d'une hypothèse, si le RTVAR est nulle, et donc qu'on ne souhaite pas inverser cette composante, on dit 
!   que les données sont nulles.

              if (rtvar(nbn+1) == 0) then     
                XX(itake)=-9999
                datatot(i,4)=-9999
              else
                XX(itake)= datatot(i,4)
              endif

              if (rtvar(nbn+2) == 0) then     
                XY(itake)=-9999
                datatot(i,5)=-9999
              else
                XY(itake)= datatot(i,5)
              endif

              if (rtvar(nbn+3) == 0) then     
                XZ(itake)=-9999
                datatot(i,6)=-9999
              else
                XZ(itake)= datatot(i,6)
              endif

              if (rtvar(nbn+4) == 0) then     
                YY(itake)=-9999
                datatot(i,7)=-9999
              else
                YY(itake)= datatot(i,7)
              endif

              if (rtvar(nbn+5) == 0) then     
                YZ(itake)=-9999
                datatot(i,8)=-9999
              else
                YZ(itake)= datatot(i,8)
              endif

              if (rtvar(nbn+6) == 0) then     
                ZZ(itake)=-9999
                datatot(i,9)=-9999
              else
                ZZ(itake)= datatot(i,9)
              endif

            end if
         end if
      end do

         write(*,*)'    There is no need to substracte the mean value to FTG data'
         write(inout,*)'    There is no need to substracte the mean value to FTG data'

!=====================================================================
! Calculation of the covariance of the data
! At that moment, we only take constant covariance over the whole data
! MP : est ce que l'on prend une covariance contante pour les 6 composantes du tenseur ?
!=====================================================================
      if(INV) then
         select case (modvar(ifil))
         case (0)
            write(inout,*)''
            write(inout,*)'    One covariance for each components'
!            ptvar3(:) = 1./(rtvar(ifil)*rtvar(ifil))
! avant on calculait ptvar3 pour toutes les données avec le même RTVAR, maintenant on
!          on veut calculer pour chaque composante, le ptvar3 avec le bon rtvar
! le RTVAR correspond maintenant à un pourcentage de la data => NON ça marchait pas pour les valeurs proches de 0
!           DU COUP on prend le pourcentage par rapport à l'amplitude min et max des données

            do j=1,6
              perc=maxval(datatot(:,j+3))-minval(datatot(:,j+3))
              perc=perc*(nbn+rtvar(nbn+ifil+j-1))/100     ! MP37
                do i=1,nptot 
                  if (datatot(i,j+3) .ne. -9999) then
                    ptvar(i+(j-1)*nptot)=1./(perc*perc)
                  else
                    ptvar(i+(j-1)*nptot)=-9999            
                  endif          
                enddo   
            enddo

            iter=0
            do j=1,npts(ifil)*6
              if (ptvar(j) .ne. -9999) then
                iter=iter+1
                ptvar3(iter)=ptvar(j)
              endif
            enddo

         case default
            write(inout,*)''
            write(*,*)''
            write(*,*)'No constant covariance not yet define'
            write(*,*)'Maybe it will be interresting to define different values &
                &according to the FTG component'
            STOP
         end select
      end if

    END SUBROUTINE GRREAD