!=====================================================================
! Module for explicite interface of subroutine BLDMAT
!=====================================================================
     !INCLUDE '../sblas/mkl_spblas.f90'
	module MOD_BLDMAT

      IMPLICIT NONE
        contains

        SUBROUTINE BLDMAT(iiter,bderi,punvar,varpar,h1,diff,npar,&
            npar1,ndat,smooth,iside,jside,ivside,jvside,&
            xb,yb,vxnodes,vynodes,par,ibove,ismooth,ilay,&
            ddvr,nbod,val_ad_s, rowAD, columnAD, nnz ,&
            AtCA, B_sp)

            USE MOD_delim
USE MOD_unit
USE MOD_inv

USE MOD_dsmooth
USE MOD_vsmooth
USE MOD_denvel
USE MOD_dgemm
USE MOD_dgemv

USE MKL_SPBLAS
!USE blas
!USE lapack

use ISO_Fortran_env, only: i4=>int32, i8=>int64



!=====================================================================
! Declaration of the in/out arguments of BLDMAT
!=====================================================================
integer,intent(in)                         :: iiter,npar,npar1,ndat
integer,intent(in)                         :: ismooth
integer, intent(in)                        :: nbod
integer,DIMENSION(:),intent(in)            :: iside,ivside
integer, DIMENSION(:),intent(in)           :: ibove,ilay
integer,DIMENSION(:,:),intent(in)          :: jside,jvside

real(kind=8),DIMENSION(:),intent(in)       :: diff,smooth,par,varpar
real(kind=8),DIMENSION(:),intent(in)       :: ddvr
real(kind=8),DIMENSION(:),intent(in)       :: vxnodes,vynodes
real(kind=8),DIMENSION(:,:),intent(in)     :: xb,yb
!real(kind=8),DIMENSION(:,:),intent(in)     :: aderi
real(kind=8),DIMENSION(:),intent(inout)    :: h1
real(kind=8),DIMENSION(:,:),intent(inout)  :: bderi
real(kind=8),DIMENSION(:),pointer          :: punvar

real(kind=8), DIMENSION(:), intent(in)     :: val_ad_s
integer, DIMENSION(:), intent(in)          :: rowAD
integer, DIMENSION(:), intent(in)          :: columnAD
integer, intent(in)                        :: nnz

type(SPARSE_MATRIX_T), intent(out)         :: B_sp
type(SPARSE_MATRIX_T), intent(out)         :: AtCA 
!=====================================================================
! Declaration of the dummy arguments of BLDMAT
!=====================================================================
integer                                    :: i,ii,j,k,l,dcont,vcont
integer,DIMENSION(8)                       :: time_miter1,time_miter2, time_multi
integer,DIMENSION(8)                       :: t_h11,t_h12, time_coo, time_fin_coo 

real(kind=8)                               :: dsum,vsum
real(kind=8)                               :: tracetot,vectrtot
real(kind=8),DIMENSION(4)                  :: traces,vectrs
real(kind=8),DIMENSION(3)                  :: itrace
real(kind=8),DIMENSION(3,4)                :: vectra,trace

!real(kind=8),DIMENSION(npar,ndat)          :: matinter !to compute aderi*punvar
real(kind=8),PARAMETER                     :: one = 1.0d0
real(kind=8),PARAMETER                     :: zero = 0.0d0
real(kind=8), DIMENSION(nnz)               :: weighted_values
type(SPARSE_MATRIX_T)                      ::Aderi_csr
type(SPARSE_MATRIX_T)                       :: Aderi_coo_t

type(SPARSE_MATRIX_T)                      :: Aderi_w_csr
type(SPARSE_MATRIX_T)                        :: Aderi_coo
integer                                    :: status
type(MATRIX_DESCR)                           :: descr

real(kind=8), dimension(npar)                :: diag_Bderi
real(kind=8), dimension(npar,npar)           :: Bd_test
!! Id sparse matrix in COO
! type(SPARSE_MATRIX_T)                        :: ID_coo
! type(SPARSE_MATRIX_T)                        :: ID_crs
! integer,DIMENSION(npar)                       :: col_id, row_id
! real(kind=8), DIMENSION(npar)                 ::  val_id
! real(kind=8), dimension(npar,npar)           :: B_Id
!=====================================================================
! Initialization of some arrays:
! VECTRA = sum of the H1 terms for each parameter type (dens., vel...)
!          2: sum of At*Cd**(-1)*DIFF
!          3: sum of At*Cd**(-1)*DIFF + smoothing
!          4: sum of At*Cd**(-1)*DIFF + smoothing + linear relation
! TRACE  = trace of the BDERI array for each parameter type
!          1: trace of At*Cd**(-1)*A
!          2: trace of At*Cd**(-1)*A + Cp**(-1)
!          3: trace of At*Cd**(-1)*A + Cp**(-1) + Dp*Cp**(-1)
!          4: trace of At*Cd**(-1)*A + Cp**(-1) + Ds*Cs**(-1) + Db*Cb**(-1)
! ITRACE = number of diagonal elts of BDERI for each parameter type
! TRACETOT = trace of the final BDERI matrix/npar
! VECTRTOT = total sum of At*Cd**(-1)*DIFF + smoothing + linear relation/npar
!=====================================================================
if(iiter.eq.1) then
trace(:,:)=0.d0
vectra(:,:)=0.d0
itrace(:)=0.d0

tracetot=0.d0
vectrtot=0.d0
traces(:)=0.d0
vectrs(:)=0.d0
end if

h1(:)=0.d0 
Bd_test(:,:) = 0.d0


write(inout,*)''
write(inout,*)'COMPUTING THE MATRIX TO BE INVERTED (BLDMAT)'
write(inout,*)'    Iteration no ',iiter

write(*,*)''
write(*,*)'COMPUTING THE MATRIX TO BE INVERTED (BLDMAT)'
write(*,*)'    Iteration no ',iiter
write(*,*)''
!=====================================================================
! Calculation of At*A (of dimension nparXnpar)
! and store the result in the array BDERI 
!=====================================================================
! write(inout,*)'    Computing the partial derivative second part'
! write(*,*)'    Computing the partial derivative second part'
!=====================================================================
! Time consumming with BLAS subroutine is much better than with MATMUL!!
!=====================================================================
! do i=1,npar
!   do j=1,i
!    do l=1,ndat
!      bderi(i,j)=bderi(i,j)+aderi(l,i)*aderi(l,j)*punvar(l)
!  end do
! if(j.ne.i) bderi(j,i)=bderi(i,j)
! end do
! end do


!============================================================
!          handel COO
!          create crs
!=============================================================

! Step 1: Prepare weighted COO values
!do i = 1, nnz
!  weighted_values(i) = sqrt(punvar(row_indices(i))) * values(i)
!end do
CALL DATE_AND_TIME(VALUES=time_coo)

!===================================================================
! Creation of hadel COO and covertion in CRS
! 
! do j=1,ndat
     ! do i=1,npar
     ! matinter(i,j)=aderi(j,i)*punvar(j)
     ! end do
! end do
! Step 1: Prepare weighted COO values
do i = 1, nnz
weighted_values(i) =(punvar(rowAD(i))) * val_ad_s(i)
end do


status = mkl_sparse_d_create_coo(Aderi_coo_t, SPARSE_INDEX_BASE_ONE,&
ndat, npar, nnz, rowAD, columnAD, weighted_values)

print *, 'COO w', status

status = mkl_sparse_d_create_coo(Aderi_coo, SPARSE_INDEX_BASE_ONE,&
ndat, npar, nnz, rowAD, columnAD, val_ad_s)

print *, 'Coo ', status

!status =0
!!! notice that in order to obtain a convertion, the finteger must be 8
!!!====================================================================
!  Notice that ADER_w_crs = A*Cd**(-1) in sparse representation
status = MKL_SPARSE_CONVERT_CSR(Aderi_coo_t,SPARSE_OPERATION_TRANSPOSE, Aderi_w_csr)

print *, 'CRS w', status

status =0
status = MKL_SPARSE_CONVERT_CSR(Aderi_coo,SPARSE_OPERATION_NON_TRANSPOSE, Aderi_csr)

print *, 'CRS ', status
status =0

status = mkl_sparse_destroy(Aderi_coo)
status = mkl_sparse_destroy(Aderi_coo_t)


!=============================================================================

write(inout,*)'    Computing the partial derivative second part'
write(*,*)'    Computing the partial derivative second part'


! !=====================================================================
! ! Calculation of At*Cd**(-1)*A.
!=======================================================================
!stat = mkl_sparse_spmm (operation, A, B, C)

!CALL DGEMM ('N','N',npar,npar,ndat,one,matinter,npar,aderi,ndat,&
!      zero,bderi, npar)


!status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, &
 !    Aderi_w_csr, Aderi_csr, AtCA) 
  !   print *, 'result  mm sparse ', status

!stat = mkl_sparse_d_sp2md (transA, descrA, A, transB, descrB, B, alpha, beta, C, layout, ldc )
!stat = mkl_sparse_d_spmmd (operation, A, B, layout, C, ldc)

status = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, Aderi_w_csr, Aderi_csr,SPARSE_LAYOUT_ROW_MAJOR, bderi, npar)
print *, 'result  mm dense bderi', status



! !=====================================================================
! ! Calculation of H1=At*Cd**(-1)*DIFF for each parameter
! !=====================================================================
 write(inout,*)'    Computing the second end member'
 write(*,*)'    Computing the second end member'

 !  CALL DGEMV('N',npar,ndat,one,matinter,npar,diff,1,zero,h1,1)
!stat = mkl_sparse_d_mv (operation, alpha, A, descr, x, beta, y)


descr%type = SPARSE_MATRIX_TYPE_GENERAL
status = mkl_sparse_d_mv (SPARSE_OPERATION_NON_TRANSPOSE, one, Aderi_w_csr, &
            descr,diff, zero, h1)

print *, 'result  mv', status
!CALL DGEMV('N',npar,ndat,one,matinter,npar,diff,1,zero,h1,1)
CALL DATE_AND_TIME(VALUES=time_fin_coo)
CALL TIMECAL(time_coo,time_fin_coo)



if(iiter.eq.1) then
     do i=1,npar
          j=1
          do while(iend(j).lt.i)
               j=j+1
          end do
          vectra(j,2)=vectra(j,2)+h1(i)
     end do
end if
write(inout,*)'    ....OK'
write(*,*)'    ....OK'

! !=====================================================================
! ! Calculation of At*Cd**(-1)*A + Cp**(-1)
! ! If (as supposed) all parameters are independent of each other, only
! ! diagonal elements of Cp**(-1) has to be added. Thus, we just have to
! ! add VARPAR elements.
! !=====================================================================
 write(inout,*)'    Computing the parameter variance part'
 write(*,*)'    Computing the parameter variance part'

!call SP_Diag(AtCA,diag_Bderi, npar)

do i=1,npar
     if(iiter.eq.1) then
        j=1
        do while(iend(j).lt.i)
             j=j+1
        end do
        trace(j,1)=trace(j,1)+bderi(i,i)
        trace(j,1)=trace(j,1)+diag_Bderi(i)
        itrace(j)=itrace(j)+1
     end if
     bderi(i,i)=bderi(i,i)+varpar(i)
     !diag_Bderi(i) = diag_Bderi(i) + varpar(i)

     if(iiter.eq.1) then
          j=1
          do while(iend(j).lt.i)
               j=j+1
          end do
          trace(j,2)=trace(j,2)+diag_Bderi(i)
     end if
  end do
write(inout,*)'    ....OK'
write(*,*)'    ....OK'
! !=====================================================================
! ! If smoothing constraint is on:
! ! Calculation of At*Cd**(-1)*A + Cp**(-1) - Ds*Cs**(-1)
! ! i.e. we add smoothing terms on the diagonal terms of BDERI
! !=====================================================================
if(ismooth.ne.0) then
   write(inout,*)'    Computing the smoothing part'
write(*,*)'    Computing the smoothing part'
if(INVD.or.INVGR) then
     CALL DSMOOTH(smooth,iside,jside,xb,yb,par,dsum,&
             dcont,bderi,h1,ismooth)
end if
if(INVV) then
     CALL VSMOOTH(smooth,ivside,jvside,par,ibove,vxnodes,&
             vynodes,vsum,vcont,bderi,h1,ismooth)
end if
write(inout,*)'    ....OK'
write(*,*)'    ....OK'
end if
! 
if(iiter.eq.1) then
   do i=1,npar
        j=1
        do while (iend(j).lt.i)
             j=j+1
        end do
        trace(j,3)=trace(j,3)+bderi(i,i)
        vectra(j,3)=vectra(j,3)+h1(i)
   end do
end if
! !=====================================================================
! ! Calculation of At*Cd**(-1)*A + Cp**(-1) - Ds*Cs**(-1) - Db*Cb**(-1)
! ! Adding the matrix connecting densities and velocities
! !=====================================================================
if((INVD.and.INVV).or.(INVGR.and.INVV)) then
   if(ibegin(3).eq.0) then
        write(*,*)''
        write(*,*)'Inconsistency between INVV, INVD (or INVGR) and ibegin(3)'
        write(*,*)'ibegin(3) should not be 0!'
        STOP 'in BLDMAT'
   end if
  write(inout,*)'    Computing the linear relation part'
         write(*,*)'    Computing the linear relation part'
         CALL DENVEL(bderi,h1,ilay,ibove,par,ddvr,nbod) !MP27 aucune correction nécessaire
         write(inout,*)'    ....OK'
         write(*,*)'    ....OK'
      end if

if(iiter.eq.1) then

do i=1,npar
j=1
do while (iend(j).lt.i)
   j=j+1
end do
trace(j,4)=trace(j,4)+bderi(i,i)
vectra(j,4)=vectra(j,4)+h1(i)
tracetot=tracetot+bderi(i,i)
vectrtot=vectrtot+h1(i)
end do

tracetot=tracetot/npar
vectrtot=vectrtot/npar

do i=1,3
if(itrace(i).gt.0) then
   do j=4,1,-1
      if(j.gt.1) then
         trace(i,j)=trace(i,j)-trace(i,j-1)
         traces(j)=traces(j)+trace(i,j)
         trace(i,j)=trace(i,j)/itrace(i)
         vectra(i,j)=vectra(i,j)-vectra(i,j-1)
         vectrs(j)=vectrs(j)+vectra(i,j)
         vectra(i,j)=vectra(i,j)/itrace(i)
      else
         traces(1)=traces(1)+trace(i,1)
         trace(i,1)=trace(i,1)/itrace(i)
      end if
   end do
end if
end do

traces(:)=traces(:)/npar
vectrs(:)=vectrs(:)/npar



!=====================================================================
! Writting some outputs
!=====================================================================
write(inout,*)''
write(inout,*)'AFTER FIRST ITERATION OF INVERSION:'
write(inout,*)'    Partial traces: Data, Par, Smooth, B:'
write(inout,'(4x,4E12.4)') traces(:)
write(inout,*)'    Total traces, matrices and vectors:'
write(inout,'(4x,2E12.4)') tracetot,vectrtot

write(*,*)''
write(*,*)'AFTER FIRST ITERATION OF INVERSION:'
write(*,*)'    Partial traces: Data, Par, Smooth, B:'
write(*,'(4x,4E12.4)') traces(:)
write(*,*)'    Total traces, matrices and vectors:'
write(*,'(4x,2E12.4)') tracetot,vectrtot

end if


status = mkl_sparse_destroy(Aderi_w_csr)
status = mkl_sparse_destroy(Aderi_csr)

     END SUBROUTINE BLDMAT

end module MOD_BLDMAT

!==================================
!    SUBROUTINE SPARSE Diagonal
!         
!    computes the diagonal of a matrix 
!    in mkl format
!===================================

     SUBROUTINE SP_Diag(A_handle, diagonal, n)

          USE MKL_SPBLAS

          IMPLICIT NONE

          integer, intent(in)                     :: n ! size of matrix 
          type(SPARSE_MATRIX_T), intent(in)       :: A_handle
          real(8), intent(out)                    :: diagonal(n)! Output diagonal
          

          !========================================================
          ! Dummy variables

          type(MATRIX_DESCR)                           :: descr
          integer                                      :: info,i
          real(8), allocatable :: e_i(:)              ! Unit vector
          real(8), allocatable :: result(:)           ! Result vector
          !=====================================================
          descr%type = SPARSE_MATRIX_TYPE_GENERAL

          ! Allocate working vectors
          allocate(e_i(n), result(n))

          do i = 1, n
          ! Create unit vector e_i
               e_i = 0.0d0
               e_i(i) = 1.0d0
          
               ! Perform matrix-vector multiplication: result = A * e_i
               info = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0d0, &
                               A_handle, descr, e_i, 0.0d0, result)
          
               if (info /= SPARSE_STATUS_SUCCESS) then
                    print *, "Error in matrix-vector multiplication"
                    return
               end if

               ! Add diagonal element to trace
               diagonal(i) = result(i)
          end do

          deallocate(e_i, result)

          end subroutine 