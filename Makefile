############################################################
#
#  makefile for fortran programm joint_inv
#
############################################################

###########################################################
# defining compilation options 
# BE CAREFULL: the options can chage with different compilers!!
############
# Main options for the f90 compiler on a SUN/UNIX environment:
# -C Check array references for out of range subscripts.
#    This option will increase the size of the executable
#    file and degrade execution performance.  It should only
#    be used while debugging.
# -e extends line source to 132 caracters
# -O optimization of the compilation.
#    Cannot be used with -C
# -free assume free-form source (f90 only)
# -fixed assume fixed form source (f90+f77)
# -fast optimized the calculations for your computer
# -ftrap explicit errors of what is following this option
#
# -Xlist check across routine for consistency. Put notes in
#        .lst file
############
# Main options for the Absoft f90/95 compiler on Mac/Darwin environment
# -v verbose mode
# -f free assume free-form source (f90/95 only)
#-----------------------------------------------------------
F77=gfortran
#OPTION = -ffree-form
#OPTION = -fdefault-real-8 -O0 -g  -fbounds-check -Wall -Wtabs -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -ftrapv -fimplicit-none 
MKLROOT = /opt/intel/oneapi/mkl/2024.2
#FFLAGS = -I$(HOME)/include/mkl/intel64/lp64 -fdefault-integer-8  -I$(MKLROOT)/include
FFLAGS = -I$(HOME)/include/mkl/intel64/lp64   -I$(MKLROOT)/include
#LIBS = -L$(HOME)/lib/intel64/libmkl_blas_ilp64.a
#LDFLAGS =   -m64 -Wl,--start-group ${MKLROOT}/lib/libmkl_gf_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
LDFLAGS =  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#F77=ifort
# option de compil anne
# =  -O3 -xhost -ipo -fp-model precise -heap-arrays
USER_OBJ_MAIN = \
		joint_inv.f90

USER_OBJ_ALL = \
		read_layer.f90 \
		read_modd.f90 \
		read_modv.f90 \
		readbco.f90 \
		statcoord.f90 \
		rdata.f90 \
		info3d.f90 \
		adnoise.f90 \
		bldmap.f90 \
		vel3.f90 \
		ttmder.f90 \
		simplex.f90 \
		rayweb.f90 \
		forw.f90 \
		frechet.f90 \
		caldt.f90 \
		findnod.f90 \
		zeroav.f90 \
		dsmooth.f90 \
		vsmooth.f90 \
		ddiffer.f90 \
		vdiffer.f90 \
		bdiffer.f90 \
		calgra.f90 \
		contbod.f90 \
		contnod.f90 \
		bcalc.f90 \
		denvel.f90 \
		dgemm.f90 \
		dgemv.f90 \
		bldmat.f90 \
		ludcmp.f90 \
		lubksb.f90 \
		invermat.f90 \
		resol.f90 \
		perturb.f90 \
		output.f90 \
		timecal.f90 \
		grread.f90 \
		ddiffergr.f90 \
		calgradio.f90


all :: costo doc

costo :: 

# compilation des programmes en C pour le calcul de tesseroid
	gcc -c lib_C/constants.c -Ilib_C
	gcc -c lib_C/geometry.c -Ilib_C
	gcc -c lib_C/grav_prism.c -Ilib_C
	gcc -c tesseroid.c -Ilib_C

	$(F77) $(FFLAGS) $(LIBS) -o costo $(USER_OBJ_MAIN) $(USER_OBJ_ALL)  $(LDFLAGS) tesseroid.o constants.o geometry.o grav_prism.o 

doc ::
	doxygen Doxyfile_costo

clean ::
	rm mod_*.mod
	rm *.o
	rm costo



#USER_OBJ_ALL = \
#		read_layer.f90 \
#		read_modd.f90 \
#		read_modv.f90 \
#		readbco.f90 \
#		statcoord.f90 \
#		rdata.f90 \
#		info3d.f90 \
#		adnoise.f90 \
#		bldmap.f90 \
#		vel3.f90 \
#		ttmder.f90 \
#		simplex.f90 \
#		rayweb.f90 \
#		forw.f90 \
#		frechet.f90 \
#		caldt.f90 \
#		findnod.f90 \
#		zeroav.f90 \
#		dsmooth.f90 \
#		vsmooth.f90 \
#		ddiffer.f90 \
#		vdiffer.f90 \
#		bdiffer.f90 \
#		calgra.f90 \
#		contbod.f90 \
#		contnod.f90 \
#		bcalc.f90 \
#		denvel.f90 \
#		dgemm.f90 \
#		dgemv.f90 \
#		bldmat.f90 \
#		ludcmp.f90 \
#		lubksb.f90 \
#		invermat.f90 \
#		resol.f90 \
#		perturb.f90 \
#		output.f90 \
#		timecal.f90
#joint_inv :: 
#	$(F77) $(OPTION) -o joint_inv $(USER_OBJ_MAIN) $(USER_OBJ_ALL)






#USER_OBJ_JOINT = \
#		read_layer.f90 \
#		readbco.f90 \
#		rdata.f90 \
#		findnod.f90 \
#		zeroav.f90 \
#		bdiffer.f90 \
#		bcalc.f90 \
#		denvel.f90 \
#		dgemm.f90 \
#		dgemv.f90 \
#		bldmat.f90 \
#		ludcmp.f90 \
#		lubksb.f90 \
#		invermat.f90 \
#		resol.f90 \
#		perturb.f90 \
#		output.f90 \
#		timecal.f90
#
#USER_OBJ_SISMO = \
#		read_modv.f90 \
#		info3d.f90 \
#		bldmap.f90 \
#		vel3.f90 \
#		ttmder.f90 \
#		simplex.f90 \
#		rayweb.f90 \
#		forw.f90 \
#		frechet.f90 \
#		caldt.f90 \
#		vsmooth.f90 \
#		vdiffer.f90 \
#		contnod.f90
#
#USER_OBJ_GRAVI = \
#		read_modd.f90 \
#		statcoord.f90 \
#		adnoise.f90 \
#		dsmooth.f90 \
#		ddiffer.f90 \
#		calgra.f90 \
#		contbod.f90
#
#joint_inv :: 
#	$(F77) $(OPTION) -o joint_inv $(USER_OBJ_MAIN) $(USER_OBJ_JOINT) $(USER_OBJ_SISMO) $(USER_OBJ_GRAVI)

