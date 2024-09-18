!=====================================================================
!=====================================================================
!   SUBROUTINE INFO3D
!=====================================================================
!=====================================================================
!> Read the last information for the 3D ray tracing velocity model
!! form the file raytrac.param
!!
!!--------
!!
!! _raytrac.param_ : 
!! * scale1 : step length along raypath
!! * stepl : step length for velocity model partial derivatives
!! * h,amp,ar,cf : number of harmonics, amplitude, amplitude ratio, cutoff
!! * nsweep = maximum number of sweep to do for the simplex algorithm
!! * tmin = minimum time difference between two deformed rays to stop the simplex algorithm     
!! * i3d : 3D raytracing 1=yes, 0=no
!v ioutext : extented output 1=yes, 0=no
!! * min_hit : minimum rays to pass through a node to invert it
!
! Called from main
! Calls none
!=====================================================================

    SUBROUTINE INFO3D

      USE MOD_unit
      USE MOD_inv
      USE MOD_vdata

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of INFO3D
!=====================================================================
! NONE in/out arguments, all variables are declared in MOD_vdata

!=====================================================================
! Declaration of the dummy arguments of INFO3D
!=====================================================================
      integer            :: ier

      character(len=80)  :: dummy
!=====================================================================
! Opening the "raytrac.param" file with the informations
! 
! scale1 : step length along raypath
! stepl : step length for velocity model partial derivatives
! h,amp,ar,cf : number of harmonics, amplitude, amplitude ratio, cutoff
! nsweep = maximum number of sweep to do for the simplex algorithm
! tmin = minimum time difference between two deformed rays to stop
!        the simplex algorithm
! i3d : 3D raytracing 1=yes, 0=no
! ioutext : extented output 1=yes, 0=no
! min_hit : minimum rays to pass through a node to invert it
!=====================================================================

      open(vpara,file='raytrac.param',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*) 'Error in opening file: , raytrac.param &
              &logical unit ', vpara,'. Stooooop in INFO3D!'
         stop
      end if

      read(vpara,*) dummy
      read(vpara,*) scale1
      read(vpara,*) dummy
      read(vpara,*) stepl
      read(vpara,*) dummy
      read(vpara,*) h,amp,ar,cf,nsweep,tmin
      read(vpara,*) dummy
      read(vpara,*) i3d
      read(vpara,*) dummy
      read(vpara,*) ioutext
      read(vpara,*) dummy
      read(vpara,*) min_hit
!=====================================================================
! Write the output in file parameter.out
!=====================================================================
      if (INVV) then
         write(*,*)''
         write(*,*)'READING PARAMETERS FOR THE RAYTRACING'
         write(inout,*)' '
         write(inout,*)'READING PARAMETERS FOR THE RAYTRACING'
         select case (ioutext)
         case (1)
            write(inout,*)' '
            write(inout,*)'    Extented output for 3D raytracing'
         case (0)
            write(inout,*)' '
            write(inout,*)'    No extented output for 3D raytracing'
         end select
         write(inout,*)' '
         write(inout,*)'    CONCERNING DELAYTIMES DATA'
         write(inout,*)'       -There are ',nstat,' stations, ',neq,' &
              &earthquakes and ',n_data,' data'
         select case (ioff)
         case (1)
            write(inout,*)'       -Correcting from the geologic factor'
         case (0)
            write(inout,*)'       -NO Correction from the geologic &
                 &factor'
         end select
      end if

      write(inout,*)' '
      write(inout,*)'    PARAMETERS FOR 3D RAYTRACING (Steck and &
           &Prothero''s programm)'
      write(inout,*)'       -Minimum ray passes to take into account &
           &nodes: ',min_hit
      select case (i3d)
      case (1)
         write(inout,*) '       -Three dimensional raytracing is ON'
      case (0)
         write(inout,*) '       -Three dimensional raytracing is OFF'
      end select
        
      write(inout,*)'       -Step length of raypath (km): ',scale1
      write(inout,*)'       -Number of Harmonics: ',h
      write(inout,*)'       -Initial Amplitude of Harmonics: ',amp
      write(inout,*)'       -Vertical/Horizontal Amplitude Ratio: ',ar
      write(inout,'(''        -Optimization Cutoff: '',f9.6)') cf
      write(inout,*)'       -Number of sweeps: ',nsweep
      write(inout,'(''        -Minimum time difference to exit the &
           &sweep: '',f9.6)') tmin
!=====================================================================
! Closing the "raytrac.param" file
!=====================================================================
      close(vpara,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the raytrac.param file, logical &
              &unit ',vpara,'. Stoooooop!'
         STOP
      end if

    END SUBROUTINE INFO3D
