!=====================================================================
! Module for explicite interface of subroutine INTMAP_VAR
!=====================================================================
	module MOD_intmapvar
	 interface
	  SUBROUTINE INTMAP_VAR(x,y,z,ip,jp,kp,vdep,vxnodes,vynodes,&
                            x_node,y_node)
	   integer,intent(out)                   :: ip,jp,kp
           integer,DIMENSION(:),intent(in)       :: x_node,y_node
           real(kind=8),intent(in)               :: x,y,z
           real(kind=8),DIMENSION(:),intent(in)  :: vdep
           real(kind=8),DIMENSION(:,:),intent(in):: vxnodes,vynodes
	  END SUBROUTINE INTMAP_VAR
	 end interface
	end module MOD_intmapvar
