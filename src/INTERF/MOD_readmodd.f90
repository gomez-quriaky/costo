!======================================================================
! Module for explicite interface of subroutine READ_MODD
!======================================================================
	module MOD_readmodd
	 interface
	  SUBROUTINE READ_MODD(zroref,roexp,ilay,xb,yb,zb,par,varpar,npar,nbod,&
                           x_block,y_block,ddep,dxcor,dycor)
	   integer,intent(in)                         :: roexp
           integer,intent(inout)                      :: npar
           integer,intent(out)                        :: nbod
           integer,DIMENSION(:),intent(inout)         :: x_block,y_block
           integer,DIMENSION(:),intent(inout)         :: ilay
           real (kind=8),DIMENSION(:),intent(inout)   :: ddep
           real (kind=8),DIMENSION(:,:),intent(inout) :: dxcor
           real (kind=8),DIMENSION(:,:),intent(inout) :: dycor
           real(kind=8),intent(in)                  :: zroref
           real(kind=8),DIMENSION(:,:),intent(inout):: xb,yb,zb
           real(kind=8),DIMENSION(:),intent(inout)  :: par,varpar
	  END SUBROUTINE READ_MODD
	 end interface
	end module MOD_readmodd
