MODULE mynrutil
	IMPLICIT NONE

	INTERFACE swap
		MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
			swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
			masked_swap_rs,masked_swap_rv,masked_swap_rm
	END INTERFACE
	
	INTERFACE imaxloc
		MODULE PROCEDURE imaxloc_r,imaxloc_i
	END INTERFACE
	
	INTERFACE assert_eq
		MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
	END INTERFACE
	
	INTERFACE outerprod
		MODULE PROCEDURE outerprod_r,outerprod_d
	END INTERFACE
	
CONTAINS
!BL
        SUBROUTINE nrerror(string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        write (*,*) 'nrerror: ',string
        STOP 'program terminated by nrerror'
        END SUBROUTINE nrerror
!BL
!BL
	SUBROUTINE swap_i(a,b)
	INTEGER, INTENT(INOUT) :: a,b
	INTEGER                :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_i
!BL
	SUBROUTINE swap_r(a,b)
	REAL(kind=8), INTENT(INOUT) :: a,b
	REAL(kind=8)                :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_r
!BL
	SUBROUTINE swap_rv(a,b)
	REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(kind=8), DIMENSION(SIZE(a))          :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_rv
!BL
	SUBROUTINE swap_c(a,b)
	COMPLEX(kind=4), INTENT(INOUT) :: a,b
	COMPLEX(kind=4)                :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_c
!BL
	SUBROUTINE swap_cv(a,b)
	COMPLEX(kind=4), DIMENSION(:), INTENT(INOUT) :: a,b
	COMPLEX(kind=4), DIMENSION(SIZE(a))          :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cv
!BL
	SUBROUTINE swap_cm(a,b)
	COMPLEX(kind=4), DIMENSION(:,:), INTENT(INOUT)  :: a,b
	COMPLEX(kind=4), DIMENSION(size(a,1),size(a,2)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_cm
!BL
	SUBROUTINE swap_z(a,b)
	COMPLEX(kind=8), INTENT(INOUT) :: a,b
	COMPLEX(kind=8)                :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_z
!BL
	SUBROUTINE swap_zv(a,b)
	COMPLEX(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
	COMPLEX(kind=8), DIMENSION(SIZE(a))          :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_zv
!BL
	SUBROUTINE swap_zm(a,b)
	COMPLEX(kind=8), DIMENSION(:,:), INTENT(INOUT)  :: a,b
	COMPLEX(kind=8), DIMENSION(size(a,1),size(a,2)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_zm
!BL
        SUBROUTINE masked_swap_rs(a,b,mask)
        REAL(kind=8), INTENT(INOUT) :: a,b
        LOGICAL, INTENT(IN)         :: mask
        REAL(kind=8)                :: swp
        if (mask) then
                swp=a
                a=b
                b=swp
        end if
        END SUBROUTINE masked_swap_rs
!BL
        SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL, DIMENSION(:), INTENT(IN)         :: mask
        REAL(kind=8), DIMENSION(size(a))          :: swp
        where (mask)
                swp=a
                a=b
                b=swp
        end where
        END SUBROUTINE masked_swap_rv
!BL
        SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(kind=8), DIMENSION(:,:), INTENT(INOUT)  :: a,b
        LOGICAL, DIMENSION(:,:), INTENT(IN)          :: mask
        REAL(kind=8), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
                swp=a
                a=b
                b=swp
        end where
        END SUBROUTINE masked_swap_rm
!BL
!BL
	FUNCTION imaxloc_r(arr)
	REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
	INTEGER                                :: imaxloc_r
	INTEGER, DIMENSION(1)                  :: imax
	imax=maxloc(arr(:))
	imaxloc_r=imax(1)
	END FUNCTION imaxloc_r
!BL
	FUNCTION imaxloc_i(iarr)
	INTEGER, DIMENSION(:), INTENT(IN) :: iarr
	INTEGER, DIMENSION(1)             :: imax
	INTEGER                           :: imaxloc_i
	imax=maxloc(iarr(:))
	imaxloc_i=imax(1)
	END FUNCTION imaxloc_i
!BL
!BL
	FUNCTION assert_eq2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN)          :: n1,n2
	INTEGER                      :: assert_eq2
	if (n1 == n2) then
		assert_eq2=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq2'
	end if
	END FUNCTION assert_eq2
!BL
	FUNCTION assert_eq3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN)          :: n1,n2,n3
	INTEGER                      :: assert_eq3
	if (n1 == n2 .and. n2 == n3) then
		assert_eq3=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq3'
	end if
	END FUNCTION assert_eq3
!BL
	FUNCTION assert_eq4(n1,n2,n3,n4,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN)          :: n1,n2,n3,n4
	INTEGER                      :: assert_eq4
	if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
		assert_eq4=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq4'
	end if
	END FUNCTION assert_eq4
!BL
	FUNCTION assert_eqn(nn,string)
	CHARACTER(LEN=*), INTENT(IN)      :: string
	INTEGER, DIMENSION(:), INTENT(IN) :: nn
	INTEGER                           :: assert_eqn
	if (all(nn(2:) == nn(1))) then
		assert_eqn=nn(1)
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eqn'
	end if
	END FUNCTION assert_eqn
!BL
!BL
	FUNCTION outerprod_r(a,b)
	REAL(kind=4), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=4), DIMENSION(size(a),size(b)) :: outerprod_r
	outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerprod_r
!BL
	FUNCTION outerprod_d(a,b)
	REAL(kind=8), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=8), DIMENSION(size(a),size(b)) :: outerprod_d
	outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerprod_d
!BL
	FUNCTION outerdiv(a,b)
	REAL(kind=8), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=8), DIMENSION(size(a),size(b)) :: outerdiv
	outerdiv = spread(a,dim=2,ncopies=size(b)) / &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiv
!BL
	FUNCTION outersum(a,b)
	REAL(kind=8), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=8), DIMENSION(size(a),size(b)) :: outersum
	outersum = spread(a,dim=2,ncopies=size(b)) + &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outersum
!BL
	FUNCTION outerdiff_r(a,b)
	REAL(kind=4), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=4), DIMENSION(size(a),size(b)) :: outerdiff_r
	outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_r
!BL
	FUNCTION outerdiff_d(a,b)
	REAL(kind=8), DIMENSION(:), INTENT(IN)   :: a,b
	REAL(kind=8), DIMENSION(size(a),size(b)) :: outerdiff_d
	outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_d
!BL
	FUNCTION outerdiff_i(a,b)
	INTEGER, DIMENSION(:), INTENT(IN)   :: a,b
	INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
	outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerdiff_i
!BL
	FUNCTION outerand(a,b)
	LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b
	LOGICAL, DIMENSION(size(a),size(b)) :: outerand
	outerand = spread(a,dim=2,ncopies=size(b)) .and. &
		spread(b,dim=1,ncopies=size(a))
	END FUNCTION outerand
!BL
END MODULE mynrutil
