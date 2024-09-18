!=====================================================================
!=====================================================================
!             SUBROUTINE RDATA
!=====================================================================
!=====================================================================
!> Read the field point location and the data from the data files 
!! of density and/or velocity
!! 
!! According to the kind of data set (codata value), this routine calls `DATA_REG` (codata=1), `DATA_IRR` 
!! (codata=2) and / or `VREAD` (codata=3)
!
! Called by main
! Calls DATA_REG,DATA_IRR,VREAD
!=====================================================================
    INCLUDE 'INTERF/MOD_datareg.f'
    INCLUDE 'INTERF/MOD_datairr.f'
    INCLUDE 'INTERF/MOD_vread.f'
    INCLUDE 'INTERF/MOD_grread.f' ! MP18

    SUBROUTINE RDATA(nfil,codata,dfile,modvar,rtvar,&
                     FX,FY,FZ,DVAL,AMASK,usedata,npts,rayparm,bazin,&
                     weight,ndat,ist,stc,ieq,ttt,punvar,varv,varg,vargr,ncomp)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_delim
      USE MOD_vdata
      USE MOD_inv

      USE MOD_datareg
      USE MOD_datairr
      USE MOD_vread
      USE MOD_grread ! MP18

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments for RDATA
!=====================================================================
      character(len=80),DIMENSION(:),intent(in) :: dfile

      integer,intent(in)                 :: nfil
      integer,DIMENSION(:),intent(in)    :: codata,modvar
      integer,DIMENSION(:),pointer       :: ieq,ist
      integer,DIMENSION(nfil),intent(out):: npts
      integer,intent(out)                :: ndat

      real(kind=8),intent(out)                :: varv,varg
      real(kind=8),DIMENSION(6),intent(out)   :: vargr      
      real(kind=8),DIMENSION(:),intent(in)    :: rtvar
      real(kind=8),DIMENSION(:),pointer       :: FX,FY,FZ,punvar
      real(kind=8),DIMENSION(:),intent(inout) :: DVAL,AMASK
      real(kind=8),DIMENSION(:),pointer       :: rayparm,bazin,weight
      real(kind=8),DIMENSION(:),pointer       :: usedata
      real(kind=8),DIMENSION(:,:),pointer     :: stc,ttt

      integer,intent(inout)                   :: ncomp

!=====================================================================
! Declaration of dummy arguments for RDATA
!=====================================================================
      integer                             :: i,ier,j,ifil,nptsgr,iter,nbn
        
      character(len=80)                   :: dummy,fdata

      real(kind=8)                        :: s,moy
      real(kind=8),dimension(6)           :: sgr,moygr
      real(kind=8),DIMENSION(:),pointer   :: ptvar1,ptvar2,ptvar3 ! MP22
      real(kind=8),DIMENSION(:),pointer   :: gdata,gFX,gFY,gFZ
      real(kind=8),DIMENSION(:),pointer       :: grFX,grFY,grFZ
      real(kind=8),DIMENSION(:),pointer       :: XX,XY,XZ,YY,YZ,ZZ 

!=====================================================================
! Initialization of the array delimiters of data
! jbegin(1), jend(1)= density
! jbegin(2), jend(2)= velocity
! MP18 jbegin(3), jend(3)= FTG
! MP18 jbegin(4), jend(4)= B parameter
!=====================================================================
      jbegin(:) = 0
      jend(:) = 0
      npts(:) = 0
!=====================================================================
! Checking and writing some outputs in parameter.out
!=====================================================================

      write(*,*)''
      write(*,*)'********************************************************'
      write(*,*)'                   READING DATA FILES'
      write(*,*)'********************************************************'
      write(*,*)''
      write(*,*)'Number of data files =',nfil
      write(inout,*)''
      write(inout,*)'***************************************************&
           &*****'
      write(inout,*)'                   READING DATA FILES'
      write(inout,*)'***************************************************&
           &*****'
      write(inout,*)' '
      write(inout,*)'Number of data files =',nfil

      if(nfil.eq.0) then
         write(*,*)'No Data file read in (nfil=0). Stooooop in RDATA!'
         STOP
      end if
        
      do ifil=1,nfil
!=====================================================================
! Open and reading the data file
!=====================================================================
         select case (codata(ifil))
!=====================================================================
! codata = 1 : Irregular gravity data file
!=====================================================================
         case (1)
            write(*,*)''
            write(*,*)'IRREGULAR GRAVITY DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,*)' '
            write(inout,*)'IRREGULAR GRAVITY DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,'(5x,''Mode of Variance: '',i1)') modvar(ifil)
            write(inout,'(5x,''Standard Deviation: '',f10.3)') rtvar(ifil)
               
            open(inddat,FILE=dfile(ifil),STATUS='OLD',IOSTAT=ier)
            fdata=dfile(ifil)
            if(ier.ne.0) then
               write(*,*)'Error in opening file ',dfile(ifil)
               write(*,*)'Stooooop in RDATA'
               STOP
            end if
            CALL DATA_IRR(fdata,ifil,modvar,rtvar,&
                 npts,ptvar1,gdata,gFX,gFY,gFZ,DVAL,AMASK)

            jbegin(1) = 1
            jend(1) = npts(ifil)

            write(inout,*)'    Number of data read: ',npts(ifil)
            write(inout,*)'    Storage of gravity data from i=',&
                 jbegin(1),' to i=',jend(1)

            close(inddat,IOSTAT=ier)
            if(ier.ne.0) then
               write(*,*)'Error in closing file ',dfile(ifil)
               write(*,*)'Stooooop in RDATA'
               STOP
            end if
!=====================================================================
! codata = 2 : Regular gravity data file
!=====================================================================
         case (2)
            write(*,*)''
            write(*,*)'REGULAR GRAVITY DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,*)' '
            write(inout,*)'REGULAR GRAVITY DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,'(5x,''Mode of Variance: '',i1)') modvar(ifil)
            write(inout,'(5x,''Standard Deviation: '',f10.3)') rtvar(ifil)
            
            open(inddat,FILE=dfile(ifil),STATUS='OLD',IOSTAT=ier)
            fdata=dfile(ifil)
            if(ier.ne.0) then
               write(*,*)'Error in opening file ',dfile(ifil)
               STOP
            end if

            call DATA_REG(fdata,ifil,modvar,rtvar,&
                 npts,ptvar1,gdata,gFX,gFY,gFZ,DVAL,AMASK)
                
            jbegin(1) = 1
            jend(1) = npts(ifil)

            write(inout,*)'    Number of data read: ',npts(ifil)
            write(inout,*)'    Storage of gravity data from i=',&
                 jbegin(1),' to i=',jend(1)

            close(inddat,IOSTAT=ier)
            if(ier.ne.0) then
               write(*,*)'Error in closing file ',dfile(ifil)
               write(*,*)'Stooooop in RDATA'
               STOP
            end if

!=====================================================================
! codata = 3 : Irregular delay times data file
!=====================================================================
         case (3)
            write(*,*)''
            write(*,*)'IRREGULAR DELAY TIMES INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,*)' '
            write(inout,*)'IRREGULAR DELAY TIMES INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,'(5x,''Mode of Variance: '',i1)') modvar(ifil)
            write(inout,'(5x,''Standard Deviation: '',f10.3)') rtvar(ifil)
            
            open(invdat,FILE=dfile(ifil),STATUS='OLD',IOSTAT=ier)
            fdata=dfile(ifil)
            if(ier.ne.0) then
               write(*,*)'Error in opening file ',dfile(ifil)
               STOP
            end if

            CALL VREAD(fdata,ifil,rayparm,bazin,weight,ttt,npts,&
                 ist,ieq,stc,ptvar2,modvar,rtvar)
            jbegin(2) = jend(1)+1
            jend(2) = jbegin(2) + npts(ifil) -1
            
            write(inout,*)''
            write(inout,*)'    Number of rays read: ',npts(ifil)
            write(inout,*)'    Storage of velocity data from i=',&
                 jbegin(2),' to i=',jend(2)
                
            close(invdat,IOSTAT=ier)
            if(ier.ne.0) then
               write(*,*)'Error in closing file ',dfile(ifil)
               write(*,*)'Stooooop in RDATA!'
               STOP
            end if

!=====================================================================
! codata = 4 : Gravity Gradients data
!=====================================================================
         case (4)
            write(*,*)''
            write(*,*)'GRAVITY GRADIENTS DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,*)' '
            write(inout,*)'GRAVITY GRADIENTS DATA INVERTED, &
                 &DATA FILE: ',dfile(ifil)
            write(inout,'(5x,''Mode of Variance: '',i1)') modvar(ifil)
            write(inout,'(5x,''Standard Deviation: '',f10.3 f10.3 f10.3 f10.3 f10.3 f10.3)') &
               rtvar(ifil),rtvar(ifil+1),&
               rtvar(ifil+2),rtvar(ifil+3),rtvar(ifil+4),rtvar(ifil+5)

            open(invdat,FILE=dfile(ifil),STATUS='OLD',IOSTAT=ier)
            fdata=dfile(ifil)
            if(ier.ne.0) then
               write(*,*)'Error in opening file ',dfile(ifil)
               STOP
            end if         

            CALL GRREAD(fdata,ifil,modvar,rtvar,npts,grFX,grFY,grFZ,&
                        XX,XY,XZ,YY,YZ,ZZ,ptvar3,ncomp)
            
            jbegin(3) = jend(2)+1
            jend(3) = jbegin(3) + npts(ifil)*ncomp -1      ! MP36

            nptsgr=npts(ifil)

            write(inout,*)'    Number of data read: ',npts(ifil)
            write(inout,*)'    Storage of FTG data from i=',&
                 jbegin(3),' to i=',jend(3)

            close(inddat,IOSTAT=ier)
            if(ier.ne.0) then
               write(*,*)'Error in closing file ',dfile(ifil)
               write(*,*)'Stooooop in RDATA'
               STOP
            end if

          nbn=0
          if (INVD) then
            nbn=nbn+1
          endif
          if (INVV) then
            nbn=nbn+1
          endif

!=====================================================================
! Else : Do not know
!=====================================================================
         case default
            write(*,*)''
            write(*,*)'Invalid data code: Should be between 1 and &
                 &4, I don''t know how to do better for the moment.'
            STOP 'in RDATA'
         end select

      end do

!=====================================================================
! Calculate total number of data and allocate memory for usedata array
!=====================================================================
!  Il y a une petite subtitlité ici car en gradio, pour une mesure, il y a 6 données !
! donc lorsque l'on va compter le nombre de points total, si on a de la gradio, on multiplie directement par 6

      ndat=0
      do i=1,nfil 
         if ( codata(i) == 4 ) then         
            ndat = ndat + npts(i)*ncomp      

         else
            ndat = ndat + npts(i)
         endif
      end do

      write(*,*)''
      write(*,*)'TOTAL NUMBER OF DATA USED (of all origins):',ndat
      write(inout,*)''
      write(inout,*)'TOTAL NUMBER OF DATA USED (of all origins):',ndat

      ALLOCATE (usedata(ndat))
      ALLOCATE (FX(ndat))
      ALLOCATE (FY(ndat))
      ALLOCATE (FZ(ndat))
      ALLOCATE (punvar(ndat))

!=====================================================================
! Stack of the data in the same file USEDATA
!=====================================================================

      if (jbegin(1).ne.0) then
         usedata(jbegin(1):jend(1)) = gdata(:)
         FX(jbegin(1):jend(1)) = gFX(:)
         FY(jbegin(1):jend(1)) = gFY(:)
         FZ(jbegin(1):jend(1)) = gFZ(:)
         punvar(jbegin(1):jend(1)) = ptvar1(:)

         moy=0.d0
         varg=0.d0
         do i=jbegin(1),jend(1)
            moy=moy+usedata(i)
         end do
         moy = moy/(jend(1)-jbegin(1)+1)
         do i=jbegin(1),jend(1)
            s=usedata(i)-moy
            varg=varg+s*s
         end do
         varg=varg/(jend(1)-jbegin(1))
         write(*,'(5x,''Variance of the gravity data: '',f8.3)') varg
         write(inout,'(5x,''Variance of the gravity data: '',f8.3)') varg
            
      end if
      if (jbegin(2).ne.0) then
         usedata(jbegin(2):jend(2)) = ttt(3,:)
         FX(jbegin(2):jend(2)) = stc(ist(:),1)
         FY(jbegin(2):jend(2)) = stc(ist(:),2)
         FZ(jbegin(2):jend(2)) = stc(ist(:),3)
         punvar(jbegin(2):jend(2)) = ptvar2(:)

         moy=0.d0
         varv=0.d0
         do i=jbegin(2),jend(2)
            moy=moy+usedata(i)
         end do
         moy = moy/(jend(2)-jbegin(2)+1)
         do i=jbegin(2),jend(2)
            s=usedata(i)-moy
            varv=varv+s*s
         end do
         varv=varv/(jend(2)-jbegin(2))
         write(*,'(5x,''Variance of the delaytime data: '',f8.3)') varv
         write(inout,'(5x,''Variance of the delaytime data: '',f8.3)') varv
      end if
        
      if (jbegin(3).ne.0) then

         iter=0     

        if (rtvar(nbn+1) .ne. 0) then
            iter=iter+1
            print*,"pour XX :: ",iter
           usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = XX(:)         
        endif

        if (rtvar(nbn+2) .ne. 0) then
           iter=iter+1
            print*,"pour XY :: ",iter
           usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = XY(:)
        endif

        if (rtvar(nbn+3) .ne. 0) then
           iter=iter+1
            print*,"pour XZ :: ",iter
           usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = XZ(:)
        endif

        if (rtvar(nbn+4) .ne. 0) then
           iter=iter+1
            print*,"pour YY :: ",iter
           usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = YY(:)
        endif

        if (rtvar(nbn+5) .ne. 0) then
           iter=iter+1
            print*,"pour YZ :: ",iter
           usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = YZ(:)
        endif

        if (rtvar(nbn+6) .ne. 0) then  
           iter=iter+1       
            print*,"pour ZZ :: ",iter
            usedata(jbegin(3)+nptsgr*(iter-1):jbegin(3)+iter*nptsgr-1) = ZZ(:)
         endif

      do i=1,nptsgr
         do j=1,ncomp
         FX(jbegin(3)+(j-1)*nptsgr+(i-1)) = grFX(i)
         FY(jbegin(3)+(j-1)*nptsgr+(i-1)) = grFY(i)
         FZ(jbegin(3)+(j-1)*nptsgr+(i-1)) = grFZ(i)
         enddo
      enddo

         moygr=0.d0
         vargr=0.d0
         punvar(jbegin(3):jend(3)) = ptvar3(:)

         write(*,'(5x,''Variance of the FTG data : '',f8.3)')
         write(inout,'(5x,''Variance of the FTG data : '',f8.3)')

      do i=1,nptsgr
         do j=1,ncomp
               moygr(j)=moygr(j)+usedata(jbegin(3)+(i-1)+nptsgr*(j-1))
         enddo
      enddo
      moygr(:)=moygr(:)/nptsgr

      do i=1,nptsgr
         do j=1,ncomp
            sgr(j)=usedata(jbegin(3)+(i-1)+nptsgr*(j-1))-moygr(j)
            vargr(j)=vargr(j)+sgr(j)*sgr(j)
         enddo
      enddo
      vargr(:)=vargr(:)/(nptsgr-1)


         iter=0     

         if (rtvar(nbn+1) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the XX component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the XX component : '',f8.3)') vargr(iter)   
         else
            write(*,'(5x,'' - for the XX component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the XX component : '',f8.3)') 0.0
         endif

         if (rtvar(nbn+2) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the XY component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the XY component : '',f8.3)') vargr(iter)
         else
            write(*,'(5x,'' - for the XY component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the XY component : '',f8.3)') 0.0
         endif

         if (rtvar(nbn+3) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the XZ component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the XZ component : '',f8.3)') vargr(iter)
         else
            write(*,'(5x,'' - for the XZ component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the XZ component : '',f8.3)') 0.0
         endif

         if (rtvar(nbn+4) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the YY component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the YY component : '',f8.3)') vargr(iter)
         else
            write(*,'(5x,'' - for the YY component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the YY component : '',f8.3)') 0.0
         endif

         if (rtvar(nbn+5) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the YZ component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the YZ component : '',f8.3)') vargr(iter)
         else
            write(*,'(5x,'' - for the YZ component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the YZ component : '',f8.3)') 0.0
         endif

         if (rtvar(nbn+6) .ne. 0) then
            iter=iter+1
            write(*,'(5x,'' - for the ZZ component : '',f8.3)') vargr(iter)
            write(inout,'(5x,'' - for the ZZ component : '',f8.3)') vargr(iter)
         else
            write(*,'(5x,'' - for the ZZ component : '',f8.3)') 0.0
            write(inout,'(5x,'' - for the ZZ component : '',f8.3)') 0.0
         endif

      end if

    END SUBROUTINE RDATA

!=====================================================================
!=====================================================================
!   SUBROUTINE DATA_REG
!=====================================================================
!=====================================================================
! Read a Regular data file (gravity) and store the fieldpoints and
! data in an array storage
!
! Called from RDATA subroutine
! Calls none
!=====================================================================
    SUBROUTINE DATA_REG(fdata,ifil,modvar,rtvar,&
                        npts,ptvar1,gdata,gFX,gFY,gFZ,DVAL,AMASK)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_inv

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments for DATA_REG
!=====================================================================
      character(len=80),intent(in)          :: fdata

      integer,intent(in)                    :: ifil
      integer,DIMENSION(:),intent(inout)    :: npts
      integer,DIMENSION(:),intent(in)       :: modvar

      real(kind=8),DIMENSION(:),pointer     :: ptvar1
      real(kind=8),DIMENSION(:),pointer     :: gdata,gFX,gFY,gFZ
      real(kind=8),DIMENSION(:),intent(inout) :: DVAL,AMASK
      real(kind=8),DIMENSION(:),intent(in)  :: rtvar

!=====================================================================
! Declaration of dummy arguments for DATA_REG
!=====================================================================
      integer                 :: NL2,NC2,num,resx,resy,&
                                 dx2e,dy2e,dx3,dy3,l1,l2,c1,c2,&
                                 i,j,NC,NL,count,take
        
      character(len=80)       :: dummy

      real(kind=8)            :: x2,y2,z2,dx2,dy2,xe,ye,x2e,y2e
      real(kind=8)            :: X0,Y0,Z0,DX,DY,mean,difx,dify
        
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: datatot,datacut
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: dataz0,dataz
      
!=====================================================================
! Read the general boundaries and the used ones in data file
!=====================================================================
      inquire(file=fdata,number=num)
      write(*,*)'    logical unit currently read:',num,',',fdata
      read(num,*) dummy
      read(num,*) x0,y0,z0,NC,NL,dx,dy,dval(ifil),amask(ifil)
      read(num,*) dummy
      read(num,*) x2,y2,z2,NC2,NL2,dx2,dy2
      write(inout,*)''
      write(inout,*)'    general boundaries of data file:'
      write(inout,*)'      x0      y0      z0      NC  NL    dx     dy'
      write(inout,'(4x,3f8.2,2i4,2f8.2)') x0,y0,z0,NC,NL,dx,dy
      write(inout,*)''
      write(inout,*)'    Used boundaries:'
      write(inout,'(4x,3f8.2,2i4,2f8.2)') x2,y2,z2,NC2,NL2,dx2,dy2
      write(inout,*)''
      if(dval(ifil).ne.0) then
         write(*,'(4x,''Value added to each data (correction): '',&
              &f8.2)') dval(ifil)
         write(inout,'(4x,''Value added to each data (correction): '',&
              &f8.2)') dval(ifil)
      else
         write(*,*)'    Mean value of measured data is substracted'
         write(inout,*)'    Mean value of measured data is substracted'
      end if
!=====================================================================
! Allocate memory for the total data storage
!=====================================================================
      ALLOCATE (datatot(NC,NL))
      ALLOCATE (dataz0(NC,NL))
!=====================================================================
! Calculation of the most SouthEastern point, and number of points
!=====================================================================
      xe=x0+dx*(NC-1)
      ye=y0-dy*(NL-1)
!=====================================================================
! Checking of the good consistency of area limits
!=====================================================================
      resx = IDINT(DABS(DMOD(dx2,dx)))
      resy = IDINT(DABS(DMOD(dy2,dy)))
      if((resx.ne.0) .or. (resy.ne.0)) then
         write(*,*)''
         write(*,*)'Extracted grid width not multiple of area width.'
         write(*,*)'Stooooop in RDATA (DATA_REG)!'
         STOP
      end if
        
      difx = x2-x0
      dify = y0-y2
      if((difx.lt.0).or.(dify.lt.0)) then
         write(*,*)''
         write(*,'(''Wrong origin coordinates, x0,y0 = '',f8.2,'','',f8.2,&
              &'' and x2,y2 = '',f8.2,'','',f8.2)') x0,y0,x2,y2
         STOP 'in RDATA (DATA_REG)!'
      end if
      if(DMOD(difx,dx).ne.0. .or. DMOD(dify,dy).ne.0.) then
         write(*,*)''
         write(*,*)'Wrong origin coordinates x2,y2 compared to dx,dy'
         STOP 'in RDATA (DATA_REG)!'
      end if

      x2e = x2+dx2*(NC2-1)
      y2e = y2-dy2*(NL2-1)
      dx2e = IDINT(xe-x2e)
      dy2e = IDINT(y2e-Ye)
      if((dx2e.lt.0).or.(dy2e.lt.0)) then
         write(*,*)''
         write(*,*)'Too big area relative to initial one'
         STOP
      end if
!=====================================================================
! Reading the data with the used boundaries x2,y2,dx2,dy2
! data reading spacing: dx3,dy3
! number of line to begin at : l1
! number of line to end at : l2
! number of column to begin at : c1
! number of column to end at : c2
!=====================================================================
      dx3= IDINT(dx2/dx)
      dy3= IDINT(dy2/dy)
      l1 = IDINT((y0-y2)/dy) + 1
      l2 = IDINT((y0-y2e)/dy) + 1
      c1 = IDINT((x2-x0)/dx) + 1
      c2 = IDINT((x2e-x0)/dx) + 1
!=====================================================================
! Allocate memory for the used (cut) data storage
!=====================================================================
      ALLOCATE (datacut(NC2,NL2))
      ALLOCATE (dataz(NC2,NL2))
!=====================================================================
! Read and store the data
!=====================================================================

      write(inout,*)' '
      write(inout,*)'    Only read from line ',l1,' to line ',l2,&
           ' with a step of ',dy3
      write(inout,*)'    Only read from column ',c1,' to column ',c2,&
           ' with a step of ',dx3
      write(*,*)' '
      read(num,*) dummy
      do j=1,NL
         read(num,*) (datatot(i,j),i=1,NC)
      end do
         
      datacut(:,:) = datatot(c1:c2:dx3,l1:l2:dy3)

!=====================================================================
! If non uniform elevation exists (ie z=99), have to read the elevation 
! value after the data
!=====================================================================
      if((z0.eq.99) .or. (z2.eq.99)) then
         read(num,*) dummy
         do j=1,NL
            read(num,*) (dataz0(i,j),i=1,NC)
         end do
         dataz(:,:) = dataz0(c1:c2:dx3,l1:l2:dy3)
      end if
!=====================================================================
! Must get rid of the mask value (amask where there is no data)
! The first loop is just to get the size of gdata and to allocate
!=====================================================================

      take = 0
      do j=1,NL2
         do i=1,NC2
            if(datacut(i,j).ne.amask(ifil)) then
               take = take+1
            end if
         end do
      end do

      npts(ifil) = take
      ALLOCATE(gdata(npts(ifil)))
      ALLOCATE(gFX(npts(ifil)))
      ALLOCATE(gFY(npts(ifil)))
      ALLOCATE(gFZ(npts(ifil)))
      ALLOCATE(ptvar1(npts(ifil)))

      take = 0
      if((z0.ne.99) .or. (z2.ne.99)) then
         do j=1,NL2
            do i=1,NC2
               if(datacut(i,j).ne.amask(ifil)) then
                  take = take+1
                  gdata(take) = datacut(i,j)
                  gFX(take) = x2 + (i-1)*dx2
                  gFY(take) = y2 - (j-1)*dy2
                  gFZ(take) = z2
               end if
            end do
         end do
      else
         do j=1,NL2
            do i=1,NC2
               if(datacut(i,j).ne.amask(ifil)) then
                  take = take+1
                  gdata(take) = datacut(i,j)
                  gFX(take) = x2 + (i-1)*dx2
                  gFY(take) = y2 - (j-1)*dy2
                  gFZ(take) = dataz(i,j)
               end if
            end do
         end do
      end if

!=====================================================================
! If DVAL=0, the mean value is calculated and substract to data
!        /=0, DAVLA is added to the data
!=====================================================================
      if(dval(ifil).ne.0.) then

         do i=1,npts(ifil)
            if(gdata(i).ne.amask(ifil)) then
               gdata(i) = gdata(i)+dval(ifil)
            end if
         end do

      else if(dval(ifil).eq.0.) then

         mean = 0.d0
         count = 0
         do i=1,npts(ifil)
            if(gdata(i).ne.amask(ifil)) then
               mean = mean + gdata(i)
               count = count+1
            end if
         end do
         mean = mean/count
    	 write(inout,'(5x,''Substract mean='',f12.3,'' mGal'')') mean
         do i=1,npts(ifil)
            if(gdata(i).ne.amask(ifil)) then
               gdata(i) = gdata(i) - mean
            end if
         end do

      end if

      DEALLOCATE (datatot)
      DEALLOCATE (datacut)
      DEALLOCATE (dataz0)
      DEALLOCATE (dataz)

!=====================================================================
! Calculation of the covariance of the data
! At that moment, we only take constant covariance over the whole data
!=====================================================================
      if(INV) then
         select case (modvar(ifil))
         case (0)
            write(inout,*)''
            write(inout,*)'    Constant covariances'
            ptvar1(:) = 1./(rtvar(ifil)*rtvar(ifil))
         case default
            write(*,*)''
            write(*,*)'MODVAR does not correspond to allowed &
                 &value in this case.'
            write(*,*)'When dealing with gravity data, only constant &
                 &covariance is allowed, Sorry!'
            STOP 'in RDATA (DATA_REG)!'
         end select
      end if


    END SUBROUTINE DATA_REG

!=====================================================================
!=====================================================================
!   SUBROUTINE DATA_IRR
!=====================================================================
!=====================================================================
! Read an irregular data file (gravity) and store the fieldpoints and
! data in an array storage
!
! Called from RDATA subroutine
! Calls none
!=====================================================================
    SUBROUTINE DATA_IRR(fdata,ifil,modvar,rtvar,&
                        npts,ptvar1,gdata,gFX,gFY,gFZ,DVAL,AMASK)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_inv

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments for DATA_IRR
!=====================================================================
      character(len=80),intent(in)            :: fdata
      
      integer,intent(in)                      :: ifil
      integer,DIMENSION(:),intent(inout)      :: npts
      integer,DIMENSION(:),intent(in)         :: modvar

      real(kind=8),DIMENSION(:),pointer       :: ptvar1
      real(kind=8),DIMENSION(:),pointer       :: gdata,gFX,gFY,gFZ
      real(kind=8),DIMENSION(:),intent(inout) :: DVAL,AMASK
      real(kind=8),DIMENSION(:),intent(in)    :: rtvar
!=====================================================================
! Declaration of dummy arguments for DATA_IRR
!=====================================================================
      character(len=80)                        :: dummy

      integer                                  :: num,i,nptot,itake

      real(kind=8)                             :: xmin,xmax,ymin,ymax,mean
      real(kind=8),DIMENSION(:,:),ALLOCATABLE  :: datatot
    
      real(kind=8)                             :: perc

      inquire(file=fdata,number=num)
      write(*,*)'    logical unit currently read: ',num,' ',fdata
      read(num,*) dummy
      read(num,*) nptot
!=====================================================================
! Allocation of the data and position array
!=====================================================================
      ALLOCATE (datatot(nptot,4))
!=====================================================================
! Read position of all x,y,z and gravity data
! and determination of the number of used data (npts)
! which is 0 at the beginning of this subroutine
!=====================================================================        
      read(num,*) dummy
      read(num,*) xmin,xmax,ymin,ymax
      read(num,*) dummy
      read(num,*) dval(ifil), amask(ifil)
      read(num,*) dummy
!=====================================================================
! If the data point belongs to the studied area and value /= amask,
! then count the point and use it
!=====================================================================
      do i=1,nptot
         read(num,*) datatot(i,1),datatot(i,2),datatot(i,3),datatot(i,4)
         if(datatot(i,1).ge.xmin.and.datatot(i,1).le.xmax) then
            if(datatot(i,2).ge.ymin.and.datatot(i,2).le.ymax) then
               if(datatot(i,4).ne.amask(ifil)) then
                  npts(ifil)=npts(ifil)+1
               end if
            end if
         end if
      end do
!=====================================================================
! Allocation of the data and position array for used data
!=====================================================================
      write(*,*)'    Number of used data for irregular file: ',npts(ifil)
      write(inout,*)''
      write(inout,*)'    Number of used data for irregular file: ',&
           npts(ifil)
      ALLOCATE (gdata(npts(ifil)))
      ALLOCATE (gFX(npts(ifil)))
      ALLOCATE (gFY(npts(ifil)))
      ALLOCATE (gFZ(npts(ifil)))
      ALLOCATE (ptvar1(npts(ifil)))
!=====================================================================
! Storage of the used data in gdata,fx,fy,fz arrays
!=====================================================================
      itake=0
      do i=1,nptot
         if(datatot(i,1).ge.xmin.and.datatot(i,1).le.xmax) then
            if(datatot(i,2).ge.ymin.and.datatot(i,2).le.ymax) then
               if(datatot(i,4).ne.amask(ifil)) then
                  itake=itake+1
                  gdata(itake) = datatot(i,4)
                  gFX(itake) = datatot(i,1)
                  gFY(itake) = datatot(i,2)
                  gfZ(itake) = datatot(i,3)
               end if
            end if
         end if
      end do

!=====================================================================
! If DVAL =0 then mean of data is substracted to the data
!        /=0 then dval is added to each data
!=====================================================================
      if(dval(ifil).ne.0.) then
         write(*,*)'    Value added to gravi data: ',dval(ifil)
         write(inout,*)'    Value added to gravi data: ',dval(ifil)
         do i=1,itake
            gdata(i) = gdata(i) + dval(ifil)
         end do
      else if(dval(ifil).eq.0.) then
         write(*,*)'    Mean is substracted to gravi data (dval=0)'
         write(inout,*)'    Mean is substracted to gravi data (dval=0)'
         mean = 0.d0
         do i=1,itake
            mean = mean + gdata(i)
         end do
         mean = mean/itake
	 write(inout,'(5x,''Substracted Mean= '',f12.3,'' mGal'')') mean
         do i=1,itake
            gdata(i) = gdata(i) - mean
         end do
      end if
!=====================================================================
! Calculation of the covariance of the data
! At that moment, we only take constant covariance over the whole data
!=====================================================================
      if(INV) then
         select case (modvar(ifil))
         case (0)
            write(inout,*)''
            write(inout,*)'    Constant covariances'
!            ptvar1(:) = 1./(rtvar(ifil)*rtvar(ifil))

            perc=maxval(datatot(:,4))-minval(datatot(:,4))
            perc=perc*(rtvar(ifil))/100
            ptvar1(:)=1./(perc*perc)
            
         case default
            write(inout,*)''
            write(*,*)''
            write(*,*)'MODVAR does not correspond to allowed &
                 &value in this case.'
            write(*,*)'When dealing with gravity data, only constant &
                 &covariance is allowed, Sorry!'
            STOP
         end select
      end if
!=====================================================================
! Deallocation of the arrays
!=====================================================================
      DEALLOCATE (datatot)

    END SUBROUTINE DATA_IRR

!=====================================================================
!=====================================================================
!   SUBROUTINE VREAD
!=====================================================================
!=====================================================================
! Read an irregular data file (velocity) and store the fieldpoints and
! data in an array storage
!
! Called from RDATA subroutine
!> Calls STATCOORD (only if input coordinates are in lat-lon)
!=====================================================================
    INCLUDE 'INTERF/MOD_statcoord.f'

    SUBROUTINE VREAD(fdata,ifil,rayparm,bazin,weight,ttt,npts,ist,&
                     ieq,stc,ptvar2,modvar,rtvar)

      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_vdata
      USE MOD_inv
      
      USE MOD_statcoord

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of VREAD
!=====================================================================
      character(len=80),intent(in)          :: fdata

      integer,intent(in)                    :: ifil
      integer,DIMENSION(:),intent(in)       :: modvar
      integer,DIMENSION(:),intent(inout)    :: npts
      integer,DIMENSION(:),pointer          :: ieq,ist

      real(kind=8),DIMENSION(:),intent(in)  :: rtvar
      real(kind=8),DIMENSION(:),pointer     :: rayparm,bazin,weight,ptvar2
      real(kind=8),DIMENSION(:,:),pointer   :: ttt,stc
!=====================================================================
! Declaration of dummy arguments of VREAD
!=====================================================================
      integer                                   :: num,i,j,maxi,latmod
      integer,DIMENSION(:),ALLOCATABLE          :: count

      character(len=80)                         :: dummy
      character(len=1)                          :: qual
      character(len=4),DIMENSION(:),ALLOCATABLE :: stn
      character(len=17)                         :: forma

      real(kind=8)                              :: qa,qb,qc,qd,tha
      real(kind=8),DIMENSION(:),ALLOCATABLE     :: off,mean

!=====================================================================
! Checking the logical unit of velocity data file.
! Should be 9
!=====================================================================
      inquire(file=fdata,number=num)
      write(*,*)'    logical unit currently read: ',num,' ',fdata
!=====================================================================
! Read the latitude and longitude of the origin
!=====================================================================
      read(num,*) dummy
      read(num,*) orlon,orlat
      write(inout,*)''
      write(inout,'(5x,''Coordinates of the origin: '',f6.2,''E'',&
           &f6.2,''N'')') orlon,orlat
!=====================================================================
! Read the station coordinates
!=====================================================================
      read(num,*) dummy
      read(num,*) latmod,tha
      read(num,*) dummy
      read(num,*) nstat
      ALLOCATE (stc(nstat,3))
      ALLOCATE (stn(nstat))
      ALLOCATE (off(nstat))
      write(inout,*)'    Number of stations (nstat) =',nstat
      do i=1,nstat
         read(num,*) stn(i),stc(i,1),stc(i,2),stc(i,3),off(i)
      end do
!=====================================================================
! If the coordinates are in latitude-longitude system, convert them
! into kilometric system centered on the origin previously read
!=====================================================================
      if(latmod.eq.1) then
         write(inout,'(5x,''Rotating angle from North (counter&
              &clockwise): '',f6.2)') tha
         CALL STATCOORD(stc,tha)
      end if

      write(inout,*)'    Stations information'
      write(inout,*)'    stn  |  stcx  |   stcy   |   stcz   |   off'
      write(inout,*)'    -------------------------------------------'
      forma='(a9,3f10.4,f10.4)'
      do i=1,nstat
         write(inout,forma) stn(i),stc(i,1),stc(i,2),&
              stc(i,3),off(i)
      end do
      read(num,*) dummy
      read(num,*) ioff
!=====================================================================
! Read number of data, earthquake, qualities
!=====================================================================
      read(num,*) dummy
      read(num,*) n_data,neq
      npts(ifil) = n_data
      read(num,*) dummy
      read(num,*) qa,qb,qc,qd
!=====================================================================
! Allocation of some arrays
!=====================================================================
      ALLOCATE (rayparm(n_data))
      ALLOCATE (bazin(n_data))
      ALLOCATE (weight(n_data))
      ALLOCATE (ttt(3,n_data))
      ALLOCATE (ieq(n_data))
      ALLOCATE (ist(n_data))
      ALLOCATE (ptvar2(n_data))
!=====================================================================
! Read earthquake informations
!=====================================================================
      read(num,*) dummy

      do i=1,n_data
         read(num,*) ieq(i),ist(i),rayparm(i),bazin(i),&
              (ttt(j,i),j=1,3),qual
         select case (qual)
         case ('a')
            weight(i)=qa
         case ('b')
            weight(i)=qb
         case ('c')
            weight(i)=qc
         case ('d')
            weight(i)=qd
         case default
            weight(i)=1.0D0
         end select
         if (ioff.eq.1) then
            ttt(1,i) = ttt(1,i) - dble(ioff)*off(ist(i))
            ttt(3,i) = ttt(1,i) - ttt(2,i)
         end if
      end do
        
!=====================================================================
! Determine the mean for each event and substract it to the data
!=====================================================================
      maxi = MAXVAL(ieq)
      if(maxi.ne.neq) then
	write(*,*)''
	write(*,*)'Error of consistency in RDATA. maxi.ne.neq!'
	STOP 'in RDATA'
      end if
      write(*,*)'    Number of event = ',maxi
      write(inout,*)'    Number of event = ',maxi
      
      write(*,*)'    Mean value is substracted to the delay times'
      write(inout,*)'    Mean value is substracted to the delay times'
      ALLOCATE(mean(maxi))
      ALLOCATE(count(maxi))

      mean(:) = 0.d0
      count(:) = 0
      do i=1,n_data
         mean(ieq(i)) = mean(ieq(i)) + ttt(3,i)
         count(ieq(i)) = count(ieq(i)) + 1
      end do
      do i=1,neq
         mean(i) = mean(i)/count(i)
      end do
      do i=1,n_data
         ttt(3,i) = ttt(3,i) - mean(ieq(i))
      end do
!=====================================================================
! Calculation of the covariance of the data
! At that moment, we only take constant covariance over the whole data
!=====================================================================
      if(INV) then
         select case (modvar(ifil))
         case (0)
            write(inout,*)''
            write(inout,*)'    Constant covariances'
            ptvar2(:) = 1./(rtvar(ifil)*rtvar(ifil))
         case (1)
            write(inout,*)''
            write(inout,*)'    Weighting applied for each data'
            do i=1,n_data
               ptvar2(i) = 1./(weight(i)*weight(i))
            end do
         case default
            write(*,*)'MODVAR does not correspond to allowed &
                 &values (0 or 1). Stooooop in RDATA!'
            STOP
         end select
      end if

      DEALLOCATE (mean)
      DEALLOCATE (count)
      DEALLOCATE (stn)
      DEALLOCATE (off)

    END SUBROUTINE VREAD   