!=====================================================================
!=====================================================================
!       SUBROUTINE SIMPLEX
!=====================================================================
!=====================================================================
! This routine takes the straight line raypath from subroutine
! rayweb and uses a simplex method of minimization to determine
! the exact 3-D raypath and travel time for a TELESEISMIC ray.
! This subroutine was written by Jess Taylor for the three dimensional
! raytracing work presented in Prothero et al. (1987).  It was 
! subsequently modified for Teleseismic rays, rather than local rays, 
! by Lee Steck (Steck and Prothero, (1990)).
! 
! Called by FORW
! calls TTIME2(VEL3), SWITCH, SORT(SWITCH)
!=====================================================================

    SUBROUTINE SIMPLEX(nrp,fstime,rp,velco)

      USE MOD_unit
      USE MOD_size
      USE MOD_vdata
      USE MOD_iloc
      USE MOD_layer

      USE MOD_switch
      USE MOD_sort
      USE MOD_ttime2
         
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of SIMPLEX
!=====================================================================
      integer,intent(in)                        :: nrp

      real(kind=8),intent(inout)                :: fstime
      real(kind=8),DIMENSION(:,:),intent(inout) :: rp
      real(kind=8),DIMENSION(:),intent(in)      :: velco

!=====================================================================
! Declaration of the dummy arguments of SIMPLEX
!=====================================================================
      integer                :: ni,i,j,m,l,har,k,n,inrp,count

      real(kind=8)           :: tt1,tt,tt2,tt3,ttr,tte,ttc,pi
      real(kind=8)           :: xj,x,y,xy,xx,yy,zz,xyz,diff
      real(kind=8)           :: cosphi,sinphi,ar1,amp1,xhar,angle,ctf
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: sine,cosine
      real(kind=8),DIMENSION(:),ALLOCATABLE  :: vdl,hdl
      real(kind=8),DIMENSION(:,:),ALLOCATABLE:: worst,second,center,&
                                                reflect,expanded,contract
!=====================================================================
! Some variable initializations
!=====================================================================
      if(ioutext .eq. 1) then
         write(inout,*)'Starting the simplex algorithm in SIMPLEX &
              &subroutine'
      end if

      tt1=fstime
      tt=fstime
      ni=nrp-1
      pi = 4.d0 * datan(1.d0)

      ALLOCATE (sine(iarsize))
      ALLOCATE (cosine(iarsize))
      ALLOCATE (vdl(iarsize))
      ALLOCATE (hdl(iarsize))
      ALLOCATE (worst(3,iarsize))
      ALLOCATE (second(3,iarsize))
      ALLOCATE (center(3,iarsize))
      ALLOCATE (reflect(3,iarsize))
      ALLOCATE (expanded(3,iarsize))
      ALLOCATE (contract(3,iarsize))
!=====================================================================
! calculate the direction cosines which will rotate
! the horizontal distortion so that it is perpendicular
! to the plane of the raypath.
!=====================================================================
      x=rp(1,nrp)-rp(1,1)
      y=rp(2,nrp)-rp(2,1)
      xy=dsqrt(x*x + y*y)
      if(xy.eq.0.0) then
         cosphi=0.
         sinphi=0.
      else
         cosphi=x/xy
         sinphi=y/xy
      end if
!=====================================================================
! calculate the vertical distortion direction cosines which
! will perturb the raypath in the direction normal to each
! segment along the raypath
!=====================================================================
      do i=2,ni
         m=i+1
         l=i-1
         zz=rp(3,m)-rp(3,l)
         yy=rp(2,m)-rp(2,l)
         xx=rp(1,m)-rp(1,l)
         xy=dsqrt(xx*xx+yy*yy)
         xyz=dsqrt(xx*xx + yy*yy + zz*zz)
         if(xyz .eq. 0.) then
            sine(i)=0.
            cosine(i)=0.
         else
            sine(i)=zz/xyz
            cosine(i)=xy/xyz
         end if
      end do
      sine(ni+1) = sine(ni)
      cosine(ni+1) = cosine(ni)
!=====================================================================
! compute proper number of harmonic to be used.  need at least
! 2N+1 points for the Nth harmonic to avoid aliasing. hh is the
! maximum allowed set in the control file.
!=====================================================================
      har=(nrp-1)/2
      if(har.gt.h) then
         har=h
      end if
      if(ioutext.eq.1) then
         write(inout,*) har,'    harmonics to deal with in SIMPLEX'
      end if
      amp1=amp
      ar1=ar
!=====================================================================
! compute starting simplex
!
! k = number of sweeps
! n = harmonic number
!=====================================================================
      vdl(1) = 0.0d0
      hdl(1) = 0.0d0
      do k=1,nsweep
         if(k.gt.1) then
            amp1=amp1/4
            ar1=ar1/4
         end if
!=====================================================================
! Starting loops over the harmonics
!=====================================================================
         do n=1,har
            xhar = dble(n)
            angle = (xhar - 0.5d0) * pi / dble(ni)
!=====================================================================
! compute the harmonic distortions
! changed loops bounds because sin(0) = 0 as defined above
!=====================================================================
            do i = 2, nrp
               xj = dble(i-1)
               vdl(i) = (amp1/xhar) * dsin (angle * xj)
               hdl(i) = (ar1/xhar)  * dsin (angle * xj)
            end do
!=====================================================================
! calculate the 'horizontally' perturbed raypath
!=====================================================================
            worst(:,1) = rp(:,1)
            worst(:,nrp) = rp(:,nrp)

            do i=1,nrp
               inrp = nrp - i+1
               worst(1,(inrp)) = rp(1,(inrp)) - sinphi*hdl(i)
               worst(2,(inrp)) = rp(2,(inrp)) + cosphi*hdl(i)
               worst(3,(inrp)) = rp(3,(inrp))
            end do
!=====================================================================
! calculate the 'vertically' perturbed raypath
!=====================================================================
            second(:,1) = rp(:,1)
            second(:,nrp) = rp(:,nrp)

            do i=1,nrp
               inrp = nrp - i+1
               second(1,(inrp)) = rp(1,(inrp)) - vdl(i)*sine(i)*cosphi
               second(2,(inrp)) = rp(2,(inrp)) - vdl(i)*sine(i)*sinphi
               second(3,(inrp)) = rp(3,(inrp)) + vdl(i)*cosine(i)
            end do
!=====================================================================
! Compute travel times
!=====================================================================
            CALL TTIME2(nrp,second,tt2,velco)
            CALL TTIME2(nrp,worst,tt3,velco)
              
!=====================================================================
! To be sure to pass one time in this loop, ctf=1.
! then, test if cut-off was reach, stop if more than 80 loops, but
! this part can be changed
!=====================================================================
            ctf = 1.
            count=0
            DO WHILE (ctf.ge.cf)
               count=count+1
               if(count.ge.80) then
                  write(*,*)count,'passes in simplex and still not &
                       &reach the cutoff (now= ',ctf,' and cutoff= ',cf,')'
                  STOP 'in SIMPLEX!'
               end if
!=====================================================================
! Sort out travel times  (rp=best vertex, then second and worst)
!=====================================================================
               CALL SORT(rp,second,worst,tt1,tt2,tt3,nrp)
!=====================================================================
! calculate the centroid
!=====================================================================
               do i=1,3
                  do j=1,nrp
                     center(i,j) = (rp(i,j)+second(i,j))/2.0d0
                  end do
               end do
!=====================================================================
! Compute reflected vertex (reflect thru avg ray)
!=====================================================================
               do i=1,3
                  do j=1,nrp
                     reflect(i,j) = 2.0d0 * center(i,j)-worst(i,j)
                  end do
               end do
               CALL TTIME2(nrp,reflect,ttr,velco)
!=====================================================================
! is the reflected vertex better than the best?
!=====================================================================
               if(ttr.le.tt1) then
!=====================================================================
! find expanded vertex
!=====================================================================
                  do i=1,3
                     do j=1,nrp
                        expanded(i,j)=3.0d0*center(i,j)-2.0d0*worst(i,j)
                     end do
                  end do
                  CALL TTIME2(nrp,expanded,tte,velco)
!=====================================================================
! is expanded the best?
!=====================================================================
                  if(tte.le.tt1) then
!=====================================================================
! keep expanded and reject the worst
!=====================================================================
                     CALL SWITCH(rp,expanded,tt1,tte,nrp)
                     CALL SWITCH(second,expanded,tt2,tte,nrp)
                     CALL SWITCH(worst,expanded,tt3,tte,nrp)
                  else
!=====================================================================
! keep reflected and reject the worst
!=====================================================================
                     CALL SWITCH(rp,reflect,tt1,ttr,nrp)
                     CALL SWITCH(second,reflect,tt2,ttr,nrp)
                     CALL SWITCH(worst,reflect,tt3,ttr,nrp)
                  end if
!=====================================================================
! if reflected fails, then begin here
!=====================================================================
               else
!=====================================================================
! is reflected better than worst?
!=====================================================================
                  if(ttr.lt.tt3) then
!=====================================================================
! accept reflected and reject worst
!=====================================================================
                     CALL SWITCH(worst,reflect,tt3,ttr,nrp)
!=====================================================================
! sort travel times
!=====================================================================
                     CALL SORT(rp,second,worst,tt1,tt2,tt3,nrp)
                  else
!=====================================================================
! find contracted vertex
!=====================================================================
                     do i=1,3
                        do j=1,nrp
                           contract(i,j)=0.5d0 * (center(i,j)+worst(i,j))
                        end do
                     end do
                     CALL TTIME2(nrp,contract,ttc,velco)
!=====================================================================
! is the contracted better than worst
!=====================================================================
                     if(ttc.lt.tt3) then
!=====================================================================
! accept contracted and reject worst
!=====================================================================
                        CALL SWITCH(worst,contract,tt3,ttc,nrp)
!=====================================================================
! sort travel times
!=====================================================================
                        CALL SORT(rp,second,worst,tt1,tt2,tt3,nrp)
                     else
!=====================================================================
! shrink all vertex half way towards the best.
!=====================================================================
                        do i=1,3
                           do j=1,nrp
                              second(i,j) = 0.5d0 * (second(i,j)+rp(i,j))
                              worst(i,j)  = 0.5d0 * (worst(i,j)+rp(i,j))
                           end do
                        end do
                        CALL TTIME2(nrp,second,tt2,velco)
                        CALL TTIME2(nrp,worst,tt3,velco)
!=====================================================================
! sort travel times
!=====================================================================
                        CALL SORT(rp,second,worst,tt1,tt2,tt3,nrp)
                     end if
                  end if
               end if
!=====================================================================
! determine if cutoff has been reached, end of DO WHILE loop
!=====================================================================
               ctf = tt3-tt1
            end do
!=====================================================================
! End of loops over the harmonics
!=====================================================================
         end do
!=====================================================================
! decide if another sweep is necessary
!=====================================================================
         diff=tt-tt1
         if(diff.lt.tmin) exit
         tt=tt1
!=====================================================================
! End of loops over the sweeps
!=====================================================================
      end do
      fstime=tt1
      if (ioutext .eq. 1) then
         write(inout,'(''    Sweep exited at k= '',i2,'' when time &
              &difference is less than '',f8.5)') k,tmin
      end if
!=====================================================================
! De-allocation of the memory for unused array
!=====================================================================

      DEALLOCATE (sine)
      DEALLOCATE (cosine)
      DEALLOCATE (vdl)
      DEALLOCATE (hdl)
      DEALLOCATE (worst)
      DEALLOCATE (second)
      DEALLOCATE (center)
      DEALLOCATE (reflect)
      DEALLOCATE (expanded)
      DEALLOCATE (contract)

    END SUBROUTINE SIMPLEX

!=====================================================================
!=====================================================================
!       SUBROUTINE SWITCH
!=====================================================================
!=====================================================================
!
! this subroutine switches array a with b and
! travel times ttw with ttb
!
! Called by SIMPLEX
! Calls none
!=====================================================================

    SUBROUTINE SWITCH(a,b,ttb,ttw,nrp)


      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of SWITCH
!=====================================================================
      integer,intent(in)                           :: nrp
      
      real(kind=8),intent(inout)                   :: ttb,ttw
      real(kind=8),DIMENSION(:,:),intent(inout)    :: a,b
!=====================================================================
! Declaration of the dummy arguments of SWITCH
!=====================================================================
      integer         :: i,j
      real(kind=8)    :: tp,temp
!=====================================================================
! Checking the array size. If wrong, stops
!=====================================================================
      if(size(a,1).gt.3 .or. size(b,1).gt.3) then
         write(*,*)'Error of array size in SWITCH subroutine:'
         write(*,*)'size(A,1) is ',size(a,1),' and size(B,1) is ',&
              size(b,1),' instead of 3. Stooooop in SIMPLEX!'
         STOP
      end if
      if(size(a,2).lt.nrp .or. size(b,2).lt.nrp) then
         write(*,*)'Error of array size in SWITCH subroutine:'
         write(*,*)'size(A,2) is ',size(a,2),' and size(B,2) is ',&
              size(b,2),' whereas nrp is ',nrp,'. Stooooop in SIMPLEX!'
         STOP
      end if
!=====================================================================
! switching the 2 arrays a and b
!=====================================================================
      do i=1,3
         do j=1,nrp
            tp=a(i,j)
            a(i,j)=b(i,j)
            b(i,j)=tp
         end do
      end do
!=====================================================================
! switching the 2 times ttb and ttw
!=====================================================================
      temp=ttb
      ttb=ttw
      ttw=temp

    END SUBROUTINE SWITCH


!=====================================================================
!=====================================================================
!       SUBROUTINE SORT
!=====================================================================
!=====================================================================
!
! this routine sorts the travel times and paths via
! the bubble method.
!
! Called by SIMPLEX
! Calls SWITCH
!
!=====================================================================
    SUBROUTINE SORT(a,b,c,tt1,tt2,tt3,nrp)

      USE MOD_switch

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of SORT
!=====================================================================
      integer,intent(in)                        :: nrp
      
      real(kind=8),intent(inout)                :: tt1,tt2,tt3
      real(kind=8),DIMENSION(:,:),intent(inout) :: a,b,c
!=====================================================================
! Beginning of the selection
!=====================================================================

      if(tt1.gt.tt2) then
         CALL SWITCH(a,b,tt1,tt2,nrp)
      end if

      if(tt2.gt.tt3) then
         CALL SWITCH(b,c,tt2,tt3,nrp)
      end if

      if(tt1.gt.tt2) then
         CALL SWITCH(a,b,tt1,tt2,nrp)
      end if

    END SUBROUTINE SORT

!=====================================================================
!=====================================================================
!        SUBROUTINE TTIME2
!
! This routine calculates the travel time of the raypath
! generated by the simplex routine.  It differs from 
! subroutine ttime only in that it works with a raypath
! dimensioned pth(3) rather than pth(3*iarsize) as in ttime.
!
! Called by SIMPLEX
! Calls VEL3
!=====================================================================
!=====================================================================

    SUBROUTINE TTIME2(nrp,pth,tt,velco)

      USE MOD_size
      USE MOD_iloc
      USE MOD_layer

      USE MOD_vel3

      IMPLICIT NONE
!=====================================================================
! Declaration of the in/out arguments of TTIME2
!=====================================================================
      integer,intent(in)                     :: nrp

      real(kind=8),intent(out)               :: tt
      real(kind=8),DIMENSION(:,:),intent(in) :: pth
      real(kind=8),DIMENSION(:),intent(in)   :: velco
!=====================================================================
! Declaration of the dummy arguments of TTIME2
!=====================================================================
      integer                               :: i,i1,ip,jp,kp

      real(kind=8)                          :: dt,x,y,z,v,psep
      real(kind=8)                          :: dz2,dy2,dx2
      real(kind=8),DIMENSION(:),ALLOCATABLE :: v1
!=====================================================================
! Initialization of some variables and V1 allocation
!=====================================================================
      tt=0.0d0
      dt=0.0d0
      ALLOCATE (v1(iarsize))
!=====================================================================
! Find velocities along the raypath.
!=====================================================================
      do i=1,nrp
         x=pth(1,i)
         y=pth(2,i)
         z=pth(3,i)
         CALL VEL3(x,y,z,v,velco,ip,jp,kp)
         v1(i)=v
      end do
!=====================================================================
! Caculate the travel time.
!
! modified 13-oct-92 to only accumulate travel time above a certain depth
!=====================================================================
      do i=2,nrp
         i1=i-1
         dz2 = pth(3,i)-pth(3,i1)
         dx2 = pth(1,i)-pth(1,i1)
         dy2 = pth(2,i)-pth(2,i1)
         psep = dsqrt( dx2*dx2 + dy2*dy2 + dz2*dz2)
         if(v1(i).eq.0 .or. v1(i1).eq.0) then
            write(*,*)'Division by 0 in TTIME2, stooooop in TTIME2!'
            write(*,*)'i,i1 = ',i,i1,' v1(i),v1(i1) = ', v1(i),v1(i1)
            STOP
         end if
         dt = psep * ( 1.0d0/v1(i) + 1.0d0/v1(i1) )
         tt = tt+dt
      end do
      tt = 0.5d0 * tt

      DEALLOCATE (v1)
      
    END SUBROUTINE TTIME2
