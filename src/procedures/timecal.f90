!=====================================================================
!=====================================================================
!   SUBROUTINE TIMECAL
!=====================================================================
!=====================================================================
! This subroutine calculate the time consumming of the calculations
! between date1 and date2
!
! Called from MAIN
! calls none
!=====================================================================

    SUBROUTINE TIMECAL(date1,date2)

      USE MOD_unit

      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of TIMECAL
!=====================================================================
      integer,DIMENSION(8),intent(inout)  :: date1,date2
!=====================================================================
! Declaration of the dummy arguments of TIMECAL
!=====================================================================
      integer                             :: dt
!=====================================================================
! Calculate the time difference in second between the two dates
! DATE(3) = day
! DATE(5) = hour (0-23)
! DATE(6) = minutes (0-59)
! DATE(7) = seconds (0-60)
! DATE(8) = milliseconds (0-999)
! dt = time in seconds between date1 and date2
!=====================================================================

      dt=(date2(3)-date1(3))*86400 + (date2(5)-date1(5))*3600 + &
           (date2(6)-date1(6))*60 + (date2(7)-date1(7))


      write(inout,*)''
      write(inout,'('' EXECUTING TIME: '',i3,''D'',i3,'':'',i3,'':'',f6.3)') &
           dt/86400,(dt-(dt/86400)*86400)/3600,&
!           (dt-(dt/86400)*86400)/60,&
           ((dt-(dt/86400)*86400)-((dt-(dt/86400)*86400)/3600)*3600)/60,&
           (1000.*(dt-(dt/86400)*86400-((dt-(dt/86400)*86400)/60)*60) &
           +(date2(8)-date1(8)))/1000.

      write(*,*)''
      write(*,'('' EXECUTING TIME: '',i3,''D'',i3,'':'',i3,'':'',f6.3)') &
           dt/86400,(dt-(dt/86400)*86400)/3600,&
!           (dt-(dt/86400)*86400)/60,&
           ((dt-(dt/86400)*86400)-((dt-(dt/86400)*86400)/3600)*3600)/60,&
           (1000.*(dt-(dt/86400)*86400-((dt-(dt/86400)*86400)/60)*60) &
           +(date2(8)-date1(8)))/1000.

    END SUBROUTINE TIMECAL
