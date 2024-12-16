!=====================================================================
!=====================================================================
!     SUBROUTINE CALDT
!=====================================================================
!=====================================================================
!> This subroutine computes :
!!
!!   - an array with the parameter indice of the hitten nodes by 3D 
!!     raytracing (IBOVE)
!!   - the frechet matrix (ADERI) according to the linear algebra of the
!!     method
!!   - The calculated delay times CALDATA (zero average) with the 
!!     new parameters PAR
! 
! Called by MAIN
! calls none
!=====================================================================

      INCLUDE 'INTERF/MOD_dgemv.f'

    SUBROUTINE CALDT(ibove,invnod,caldata,vels,ieq,aderi,par,&
                     val_ad_s,columnAD,rowAD,idx_sp_A,nnz,ndat,npar)

      USE MOD_unit
      USE MOD_delim
      USE MOD_layer

      USE MOD_dgemv
      
      IMPLICIT NONE

!=====================================================================
! Declaration of the in/out arguments of CALDT
!=====================================================================
      integer,intent(in)                         :: invnod
      integer,DIMENSION(:),intent(inout)         :: ibove
      integer,DIMENSION(:),pointer               :: ieq

      real(kind=8),DIMENSION(:),intent(in)       :: par
      real(kind=8),DIMENSION(:),intent(inout)    :: caldata
      real(kind=8),DIMENSION(:,:),intent(inout)  :: aderi
      real(kind=8),DIMENSION(:,:,:,:),intent(in) :: vels

      real(kind=8),DIMENSION(:),intent(inout) 	 :: val_ad_s
	   integer, DIMENSION(:),intent(inout)   :: columnAD
	   integer, DIMENSION(:), intent(inout)  :: rowAD
      integer, intent(inout)                 :: idx_sp_A
      integer,intent(in)					 :: nnz
	   integer,intent(in)					 :: ndat
	   integer, intent(in)					:: npar

!=====================================================================
! Declaration of the dummy arguments of CALDT
!=====================================================================
      integer                               :: i,ii,j,jj,ier,ib,jb,m,nray
      integer                               :: no_event,nnod,num 
      integer                               :: k,x,y,maxi
      integer,DIMENSION(:),ALLOCATABLE      :: count

      !index for Sparse Aderi density
      integer                                :: ind_sp_v
      integer                                :: check
      !tets if sparse matix*vec give the same (compare with caldata)
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: cald_sp, cald_blas
      real(kind=8)                            :: one_a, zero_a

      real(kind=8)                          :: der_slow,vinit,vpert
      real(kind=8),DIMENSION(:),ALLOCATABLE :: mean
!=====================================================================
! Initialization of ibove and caldata arrays to 0.
!=====================================================================
      ibove(:) = 0
      caldata(jbegin(2):jend(2)) = 0.d0

   !constant for sparse blas multiplication
      one_a = 1.0d0
      zero_a = 0.0d0

      write(*,*)''
      write(*,*)'FORWARD CALCULATION OF SYNTHETIC DELAY TIMES'
      write(inout,*)''
      write(inout,*)'FORWARD CALCULATION OF SYNTHETIC DELAY TIMES'
!=====================================================================
! Opening the file inmat (synthe.fre)
!=====================================================================
      open(inmat,file='synth.fre',status='old',iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in opening the file SYNTHE.FRE, logical unit ',&
              inmat,'. Stooooop in CALDT!'
         STOP
      end if
!=====================================================================
! Check if reading parameters are non null
!=====================================================================
      if(ibegin(2).eq.0) then
         write(*,*)''
         write(*,*)'Inconsistency in CALDT (ibegin(2)=0). Stooooop!'
         STOP
      end if
!=====================================================================
! Check counter of data (jb) to starting value
! Every measured delay time increases this counter by 1
!=====================================================================
      if(jbegin(2).ne.0) then
         jb = jbegin(2)-1
      else
         jb = 0
      end if
!=====================================================================
! Starting reading the informations in the synthe.fre file (INMAT)
! INMAT is organized in blocks, with header line of:
!       - the number of node in the velocity model (m)
!       - the number of rays passing through this node (nray)
!       - the new number of nodes to be inverted (nnod)
! Then for each block:
!       - the number of the event(no_event)
!       - the partial derivative(der_slo)
!=====================================================================
      do i=1,invnod
         read(inmat,*,end=100) m,nray,nnod 
!=====================================================================
! Correspondance between node number and parameter number is stored 
! in ibove
! ibove(no_of_the_node) = no_of_the_parameter in PAR array if node is hit
!                       = 0                                if node is not hit
!=====================================================================
         ib = ibegin(2)+m-1
         ibove(m) = ib
!=====================================================================aderi
! Find the initial value of the velocity for this node and store it
!=====================================================================
         num = 0
 exter:  do k=1,nznode-1
            do y=1,nynode
               do x=1,nxnode
                  num = num+1
                  if(num.eq.m) then
                     vinit = vels(x,y,k,1)
                     exit exter
                  end if
               end do
            end do
         end do exter
!=====================================================================
! For each ray passing through the node, calculate the new perturbed
! velocity (vpert) and the corresponding delay time (caldata)
! In terms of dimension (L=length, T=time):
! der_slo = [L]
! par = vpert = vinit = [L/T]
! aderi = [T**2/L]
! caldata = [T]
!=====================================================================
         do j=1,nray
            read(inmat,*) no_event,der_slow
            jj=jb+no_event
            vpert = vinit + par(ib)
            aderi(jj,ib) = -der_slow/(vpert*vpert)
            ind_sp_v = idx_sp_A + (m-1)*(jend(2)-jb) + no_event                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
            val_ad_s(ind_sp_v) =  -der_slow/(vpert*vpert)
            rowAD(ind_sp_v) = jj
            columnAD(ind_sp_v) = ib
            caldata(jj) = caldata(jj) + aderi(jj,ib)*par(ib) !replace by sparse matrix*vector
         end do
      end do
!=====================================================================
! Closing the file inmat (synthe.fre)
!=====================================================================
 100  close(inmat,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file SYNTHE.FRE, logical unit ',&
              inmat,'. Stooooop in CALDT!'
         STOP
      end if
!==================================================

! the following multiplication was equivalent for a vector equal to 0,
      !check if caldata full give the same


   !CALL mkl_dcoomv('N',ndat, npar,one_a,'G**F',val_ad_s,rowAD,columnAD,&
   !                  nnz,par,zero_a,caldata)
    !===========================================================  


! test caldata and sparce matrix*vec
    !  write(*,*) 'fin loop, begin mkl'
     ! check =0


      !allocate(cald_sp(ndat))
      !cald_sp(:) = 0.0d0
      !one_a = 1.0d0
      !zero_a = 0.0d0
      !call mkl_dcsrmv('N', m, n, alpha, 'G**F', val_a, row_a_s, , x, beta, y)
      !CALL mkl_dcoomv('N',ndat, npar,one_a,'G**F',val_ad_s,rowAD,columnAD,&
       !              nnz,par,zero_a,cald_sp)


      !do j=jbegin(2),jend(2)
      !   if (caldata(j).ne.cald_sp(j)) then
      !      check = check +1
      !      print *, 'no sparsze result', caldata(j), 'sparse ', cald_sp(j)
      !      !EXIT
      !   end if
      !end do

      !write(*,*)'total diff ', check , &
       !       ' total ', jend(2) - jbegin(2)
      

!=====================================================================
! The average delay times for every event is zero, thus the
! synthetic delay times should be normalized the same way.
! The mean for each event is subtracted to delay times
! ieq(ii) is the event number whose ray ii belongs to (ii=1,n_data)
!=====================================================================
      maxi = MAXVAL(ieq)
      ALLOCATE(mean(maxi))
      ALLOCATE(count(maxi))

      mean(:) = 0.d0
      count(:) = 0
      do i=jbegin(2),jend(2)
         ii=i-jbegin(2)+1
         j=ieq(ii)
         mean(j) = mean(j) + caldata(i)
         count(j) = count(j) + 1
      end do

      do i=1,maxi
         if(count(i).ne.0) then
            mean(i) = mean(i)/count(i)
         else
            write(*,*)''
            write(*,*)'WARNING: Event No ',i,' has no rays. Check &
                 &data file.'
         end if
      end do

      do i=jbegin(2),jend(2)
         ii=i-jbegin(2)+1
         j=ieq(ii)
         caldata(i) = caldata(i) - mean(j)
      end do

      ! uodate idx_sp_A
      idx_sp_A = idx_sp_A + ind_sp_v
      DEALLOCATE(mean)
      DEALLOCATE(count)

    END SUBROUTINE CALDT
