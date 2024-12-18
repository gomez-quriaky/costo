! MP
!> \mainpage
!! The algorithm is based on the `joint_inv` algo from C. Tiberi. We will briefly detail how to use the code here and the main informations. 
!!
!! See Tiberi et al. (2003, 2018) for more details about the joint inversion between tomography and gravity, 
!! 	  Plasman et al (2020) for more details about the gravity gradients inversion and 
!!    Plasman et al (in prep) for a joint inversion of the three methods. 
!! 
!! See the `joint_inv.f90` program for more details on the several improvements and corrections
!! realized
!! 
!! \section aa 1. General descriptions
!!
!! The overall method is based on the theory of Zeyen & Achauer (1997), <em>Joint inversion of teleseismic
!!    delay times and gravity anomaly data for regional structures
!!    in Upper mantle heterogeneities from active and passive seismology.</em>
!!    K.Fuchs. Then, it has been add :
!!
!! * a 3D ray tracing from <em>Steck & Prothero (1995)</em> for the teleseismic part.
!!
!! * a GRAVITY computation from BLAKELY subroutine (rigth rectangular prism)
!!    (see Li & Chouteau, 1998). Three dimensional gravity modeling in
!!     all space, Surveys in Geophysics, 19, 339-368,for a methods overview).
!!
!! * a Gravity gradient computation for rectangular prism, with the TESSEROID algorithm from Uieda et al. 2016.
!!
!! The inversion of separate dataset can be resumed to the minimization of the difference between the observed and the synthetic data 
!!    to retrieved the best 3-D distribution of rectangular density prisms (gravity or gravity gradient) or velocity nodes (tomography).
!!    The minimization is realized with a Bayesian optimization (Menke 1984) to include a priori and a posteriori information. 
!!    For the gravity gradient part, we can considere all the gradient components, or a combination of them, or a single one.
!!
!! We suppose a Gaussian probability density function for both data and parameters distribution, and use covariance matrices 
!!    to weight the data and the parameters distribution.
!!
!! At each iteration, a new model (density and/or velocity) is computed with respect to the previous density/velocity distribution p0. 
!!    The convergence of the inversion is verified 
!!    through the RMS decrease. The process is stopped either after a determined number of iterations, or when a given threshold (ε) is reached 
!!    [(dobs − dcalc) ≤ ε] or when the difference between two consecutive iterations is small enough (≤ε/100, Zeyen & Achauer 1997).
!!
!! At the end, the resolution matrix is also computed (Zeyen & Achauer 1997; Tiberi et al. 2003). Its diagonal terms give 
!!    the resolution of each separate parameter. And the off-diagonal terms quantify the interdependency between each parameter. 
!!    This information is more complete but much more difficult to analyse because the more parameters, the biggest the matrix. 
!!    That is why, we will hereafter only use the diagonal terms. Classicaly in this type of gravity inversion, especially with synthetic 
!!    inversion of homogeneous data, the resolution is quite homogeneous inside a same layer and decrease with depth. The detail of the parameters 
!!    equation can be found in Zeyen & Achauer (1997) and Tiberi et al. (2003).
!!
!!
!!
!! Few points :
!! * The depth axis (z) is POSITIVE DOWN!
!!
!! * The densities are in g/cm3
!!
!! * The velocities are in km/s
!!
!! * The b coefficients, needed for the joint inversion between density and tomography are in (km/s)/(g/cm3)
!! 
!!
!! \section cc 2.2 Input files
!!
!! \subsection bb 2.1 File parameter.inp
!!
!! This file contains all the informations about the inversion parameters (smoothing, nature of the data, etc.)
!!    as well as the informations about the model(s), the dataset(s), inversion type, etc.
!! 
!! It is composed as follow :
!!
!! <b>1 - TITLE</b> <em>(A80)</em> \n
!!   title of the inversion
!!
!! <b>2 - INV,INVD,INVV 	</b>			        <em>3(l2,5x),l2</em> \n
!!   inv  =T inversion, 
!!        =F direct problem \n
!!   invd =T density inverted (or if INV=F, direct problem) \n
!!	 if =F    "    not inverted \n
!!   invv =T velocity inverted (or if INV=F, direct problem) \n
!!   if	=F    "     not inverted
!!
!! <b>3 - ITER,EPSI	</b>	 				<em>free format</em> \n
!!   iter = number of iterations \n
!!   epsi = maximum allowed difference between model and data
!!	   for stopping the iterations before iter. \n
!!	  epsil = maximum difference between two consecutive
!!                  iterations for stopping the iterations before
!!                  MAXITE.
!!
!! <b>3b - CORREL THRESHOLD		</b>			<em>free format</em> \n
!!    thrhold = threshold for the B value that relates the velocity
!!              and the density. If the correlation calculated within the BCALC subroutine is
!!		greater than thrhold, then the new B-value is taken. Else, the
!!		old  one is conserved.
!!
!! <b>4 - ISMOOTH,SMOOTH(1),SMOOTH(2)	</b>			<em>free format</em> \n
!!   ismooth=0 no smoothest model algorithm \n
!!	   =1 smoothest model for all blocks \n
!!   smooth(1) smoothing factor for density (variance) \n
!!   smooth(2) smoothing factor for velocity (variance)
!! 
!! <b>5 - REGULARIZATION,LAMBDA FACTOR	</b>		<em>free fromat</em> \n
!!   regul=0 no regularization \n
!!        =1 regularization by lambda factor \n
!!   lambda: multiply the diagonal elts of BDERI with this factor
!!
!! <b>6 - NOISED,NOISEV,SIGNOISD,SIGNOISV</b>
!!   Always read, but only taken into account for direct problem 
!!   (i.e. INV=.FALSE.) \n
!!   noised = 1 Adding random gaussian noise to synthetic gravity \n
!!            0 No gaussian noise added \n
!!   noisev = 1 Adding random gaussian noise to synthetic delay times \n
!!            0 No gaussian noise added \n
!!   noisegr = 1 Adding random gaussian noise to synthetic gradio \n
!!            0 No gaussian noise added \n
!!   signoisd = Sigma for the gaussian curve (density) \n
!!   signoisv = Sigma for the gaussian curve (velocity) \n
!!   signoisgr = Sigma for the gaussian curve (gradio) \n
!!
!! <b>7 - ZROREF,ROEXP              </b>          		<em>free format</em> \n
!!   Parameters for calculating a depth dependent covariances of
!!   densities: Var=VAR0*(Z/ZROREF)**ROEXP \n
!!   with ROEXP>0, it favours the density contrast in deep layers. \n
!!   with ROEXP<0, it favours the density contrast in shallow layers. \n
!!   zroref= Reference depth \n
!!   roexp= Exponent \n
!!
!! <b>8 - NLAYER				</b>		<em>Free Format</em> \n
!!   nlayer = Number of density layers to be used in the inversion
!!
!! <b>9 - NLR,BINIT,VARB,DDVR,DRMAX   </b>  <em>Free Format</em> \n
!!   NLR   : Number of (density) layer \n
!!   BINIT = Starting value of dV/dRho (B-coeff) in (km/s)/(g/cm3) \n
!!   VARB  = Standard deviation of B-coefficient \n
!!   DDVR  = Variance of velocity-density relation \n
!!   DRMAX = Maximum variation of density in this layer (allowed RHO is between +/-DRMAX) (g/cm3) \n
!!   The calculation is done with relative density contrasts.
!!
!!<b>10 - DENSMOD</b> \n
!!   Name of the density model file                      <em>Free Format</em> \n
!!   Open it only if needed
!!
!!<b>11 - VELMOD </b>\n
!!   Name of the velocity model file                    <em>Free Format</em> \n
!!   Open it only if needed
!!
!!<b>12 - NFIL</b>						<em>Free Format</em> \n
!!   nfil : number of input data files
!!
!!-Following blocks (12 to 14) are repeated for any number of files nfil \n
!! <b> WARNING: IF GRAVITY AND VELOCITY FILES ARE INPUT, GRAVITY SHOULD BE
!!          THE FIRST... Because of reading procedure in
!!          rdata subroutine.</b>
!!
!!<b>13 - CODATA		</b>				<em>Free Format</em> \n
!!   codata : code for the data format \n
!!          =1 irregular gravity data \n
!!          =2 regular gravity data \n
!!          =3 delay times \n
!!          =4 gravity gradient data \n
!!
!!<b>14 - DFILE		</b>				<em>(A80)</em> \n
!!   Name of the data file to be read
!!
!!<b>15 - MODVAR,RTVAR  </b>              			<em>Free Format</em> \n
!!   modvar =0 Constant covariance, sigma is RTVAR \n
!!	      Standard deviation = RTVAR, \n
!!             Covariance = 1/Sigma**2 \n
!!	   =1 Used only with the delaytime data. \n
!!             Weighting the data with quality factors.
!!             There are 4 quality factors for 4 sigma values
!!             In this case, RTVAR is read but not taken into account
!! 
!! \subsection cc 2.2 Input models
!! \subsubsection ss Density model
!! 
!! The density model is the same for gravity and/or gravity gradient inversions
!!
!! * DUMMY <em>(a80)</em>
!!
!! * nxbloc,nybloc : maximum nb of x,y blocks              <em>(Free Format)</em>
!!
!! * DUMMY <em>(a80)</em>
!!
!! * DDEP(NLAYER+1): depth of the interfaces of the model  <em>(Free Format)</em>
!!
!! * DUMMY <em>(a80)</em>
!!
!! * D0(NLAYER) : reference density for each layer         <em>(Free Format)</em>
!!
!!  for each layer then:
!!
!! * DUMMY, x_block,y_block  : number of EW(x) and NS(y) blocks
!!
!! * dxcor(x_block+1)        : coordinates of the x corners of the blocks
!!
!! * dycor(y_block+1)        : coordinates of the y corners of the blocks
!!
!! * density(x_block,y_block,ddep)  : each line is EW, each column is NS
!!
!! * dvari(x_block,y_block,ddep)  : (std.dev.)each line is EW, each column is NS
!!
!! \subsubsection sss Velocity model
!!
!! * DUMMY <em>(a80)</em>
!!
!! * nxnode,nynode,nznode : nb of x,y,z nodes              <em>(Free Format)</em>
!!
!! * DUMMY <em>(a80)</em>
!!
!! * VXNODES(nxnode) : nodes coordinates in x-direction    <em>(Free Format)</em>
!!
!! * VYNODES(nynode) : nodes coordinates in y-direction    <em>(Free Format)</em>
!!
!! * VZNODES(nznode) : depth of the layers of the model    <em>(Free Format)</em>
!!
!! for each layer then :
!!
!! * DUMMY
!!
!! * vels(x_node,y_node,z_node,1)   : each line is EW, each column is NS
!!
!! * vvari(x_node,y_node,z_node,1)  : each line is EW, each column is NS
!!
!!
!! ------
!!
!! __WARNING__ 
!! 1. The number of layers in the velocity model should be greater than NLAYER, 
!! so that the 3D raytracing could take into account the last layer.
!! 2. Both models should have the same limits in X and Y axes. Otherwise, the B-coeff 
!! calculation does not know how to deal with. (error from VEL3 used in BLOCKMEAN subroutine).
!!
!! \subsection ca 2.3 Data files
!! Data files can be  :
!! * for the gravity data : regular or irregular data set
!! * for the tomo : delay times data
!! * for the gravity gradient : regular or irregular data set
!!
!! \subsubsection sq Regular gravity data file
!!
!! <b>DUMMY</b>\n
!! <b>x0,y0,z0,NC,NL,dx,dy,dval,amask	</b>		<em>Free Format</em> \n
!!    general boundaries of the data
!!     * X0,Y0: Coordinate of the NW corner of the area (Y=North,
!!           X=East) (km)
!!     * Z0   : Elevation (km).  
!!           if Z0 = 99, read elevation for each data points after data 
!!           if z0/= 99, take Z0 has the common elevation value for data 
!!     * NC   : Number of N-S profiles 
!!     * NL   : Number of points on every profile (inEW direction) 
!!     * DX,DY: Distance between the points in E-W / N-S direction. 
!!     * DVAL : =0 Mean value of the measured data is substracted  \n
!!          /=0 DVAL is added to the measured data in order to
!!           get rid of a constant offset. 
!!     * AMASK: Value for points with no data \n
!!
!!
!! <b>DUMMY</b>\n
!! <b>x2,y2,z2,NC2,NL2,dx2,dy2  </b>                    <em>Free Format</em> \n
!!    used boundaries for the inversion, should be multiples of
!!    x0,y0,z0,NL,NC,dx,dy \n
!!    With same rules for Z2 as those for Z0. \n
!! <b>DATA</b> \n
!!    Each column is a NS profile, each line a WE profile
!! if Z0=99 or Z2=99: \n
!! <b>DUMMY</b> \n
!! <b>ELEVATIONS </b>                                  <em>Free Format</em> \n
!!    Each column is a NS profile, each line a WE profile
!!
!!
!! \subsubsection sq Irregular gravity data file
!! <b>DUMMY</b>\n
!! NDATA : Number of data in this file                    <em>Free Format</em> \n
!! <b>DUMMY</b>\n
!! XMIN,XMAX,YMIN,YMAX: Delimitation of the used area     <em>Free Format</em> \n
!! <b>DUMMY</b>\n
!! DVAL AMASK:                                            <em>Free Format</em>
!!    * DVAL : =0 Mean value of the measured data is substracted  \n
!!          /=0 DVAL is added to the measured data in order to
!!           get rid of a constant offset.
!!    * AMASK: Value for points with no data
!! <b>DUMMY</b>\n
!! X Y Z DATA: Coordinates and value of the data point    <em>Free Format</em> \n
!!    X,Y,Z in km \n
!!    DATA in mGal
!!
!!
!! \subsubsection sz Delay times data file
!! <b>DUMMY</b>\n
!! <b>ORLON ORLAT</b> : Longitude and latitude of the middle point
!!               used only if LATMOD=1                        <em>Free Format</em> \n
!! <b>DUMMY</b>\n
!! <b>LATMOD THA</b>
!!    * latmod = 1 : input stations coordinates in longitude/latitude\n
!!           = 0 : input stations coordinates in km
!!    * tha = rotating angle (counterclockwise from North) to aligned the
!!          station coordinates with nodes. Used only if latmod=1
!!
!! <b>DUMMY</b>\n
!! <b>NSTAT</b> = number of stations\n
!! <b>STN STC(1:3) OFFS</b> : Information for the teleseismic stations
!!    * STN      = name of the station                                 <em>(A4)</em> \n
!!    * STC(1:3) = x,y,z coordinates of the station             <em>Free Format</em> \n
!!    * OFFS     = time offset to substract to data if ioff=1   <em>Free Format</em> \n
!!
!! <b>DUMMY</b>\n
!! <b>IOFF </b>= 1 : substract OFFS to each data                     <em>Free Format</em> \n
!!      = 0 : no time correction made
!! <b>DUMMY</b>\n
!! <b>N_DATA NEQ</b>                                                 <em>Free Format</em> \n
!!    * n_data : number of data (rays)
!!    * neq    : number of earthquakes
!!
!! <b>DUMMY</b>\n
!! <b>QUAL(1:4)</b> : Quality of the data (Std.Dev.)                 <em>Free Format</em> \n
!! <b>DUMMY</b>\n
!! <b>IEQ IST RAYPARM BAZIN TOBS TPRED TDIFF QUAL</b>               <em>Free Format</em> \n
!!    * ieq : number of the earthquake
!!    * ist : number of the station
!!    * rayparm : ray parameter of the ray
!!    * bazin : backazimuth of the ray
!!    * tobs : observed arrival time for the ray (should be = tdiff)
!!    * tpred : predicted arrival time for the ray (NOT USED AT PRESENT = 0)
!!    * tdiff : tobs-tpred, DATA READ IN RDATA SUBROUTINE
!!    * qual : quality of the data according to QUAL(1:4)
!!
!!
!! \subsubsection st Gravity gradient data file
!! <b>DUMMY</b>\n
!! NDATA : Number of data in this file                    <em>Free Format</em> \n
!! One data is composed of the six tensor components
!! <b>DUMMY</b>\n
!! XMIN,XMAX,YMIN,YMAX: Delimitation of the used area     <em>Free Format</em> \n
!! <b>DUMMY</b>\n
!! <b>DATA (X Y Z XX XY XZ YY YZ ZZ)</b>    <em>Free Format</em> \n
!! Coordinates of the data point and gravity gradient value in E
!!
!!
!! __WARNING__ 
!! With this approach, the second derivatives of gz are not defined for points located inside or at the surface of a prism.
!!
!! \section XX 3. Sorting programs
!!
!! \subsection oo 3.1 For tomography 
!! * `read_modv.f90` : reads the velocity body information in files velmod
!! * \ref vread : Read an irregular data file (velocity)
!! * `info3d.f90` : Read the last information for the 3D ray tracing velocity model
!! * `bldmap.f90` : points to node indices for any x,y,z for the 3D raytracing velocity model
!! * \ref vcoefs : computes the coefficients for the linear velocity interpolation
!! * `frechet.f90` : calculates the frechet matrix and stored the number of rays passing through each node
!! * `forw.f90` which traces ray for a single event to single station.
!! Thus it is called foreach station per event (calls several subroutines)
!! * \ref makemap : looks for the points of each layer where rays pass through
!! * `caldt.f90` : computes 1) the parameter indice of the hitten nodes and 2) calculated delay times CALDATA
!! (zero average) with the new parameters PAR
!! * `vdiffer.f90` : calculates the difference between observed and calculated
!! delaytime data and uses it to compute the RMS
!! * `vsmooth.f90` : computes the smoothing constrained for BDERI matrix
!! for the velocity bodies
!! * `contnod.f90` :looks for contiguous velocity nodes, in order to store
!! the parameters of the smoothing matrix
!!
!! \subsection oi 3.2 For gravity 
!! * `read_modd.f90` : read of density model information in file densmod
!! * \ref data_reg :  Read a regular data file (gravity)
!! * \ref data_irr :  Read a irregular data file (gravity)
!! * `ddiffer.f90` : calculates the difference between observed and calculated
!! density data and uses it to compute the RMS
!! * `dsmooth.f90` : computes the smoothing constrained for BDERI matrix
!! for the density bodies
!! * `calgra.f90` : compute the frechet matrix and the response model
!! * \ref gbox : computes the vertical attraction of a rectangular prism for one 
!! observation point.
!! * `contbod.f90` : looks for contiguous bodies, in order to store
!! the parameters of the smoothing matrix
!!
!! \subsection op 3.3 For gravity gradient 
!! * `calgradio.f90` : compute the gravity gradient response of the density model
!! * `calgra_gz.f90` : compute the gravity effect but with the tesseroid code
!! * `ddiffgr.f90` : compute the differece between syn and obs 
!! * `grread.f90` : read the gravity gradient file data
!! * + all the MOD_* files associated
!! * + all the tesseroid files
!! * + bldmat.f90, ddiffer.f90, joint_inv.f90, output.f90, perturb.f90, rdata.f90, read_layer.f90, resol.f90, zeroav.f90, vdiffer.f90 have been modified to integrate the gravity gradient
!!
!! \subsection ol 3.4 For joint inversion 
!! * `read_layer.f90` : Read the maximum number of layer, x, y, blocks or nodes
!! in files dens.mod and / or vel.mod
!! * `readbco.f90` : reads in file parameter.inp the initial values of B-coeff
!! * `rdata.f90` : read the field point location and the data from the data files.
!! * `findnod.f90` : creates file blocnod where for each density block the corresponding number 
!! of the velocity nodes including information are stored.
!! * `bdiffer.f90` : calculates the difference between observed and calculated
!! density data and uses it to compute the standar deviation.
!! * \ref blockmean : computes the average velocity for a density block.
!! * `bcalc.f90` : computes the B-coefficient of the linear relation
!! between density and velocity.
!! * `bldmat.f90` : calculates the matrix to invert.
!! * `denvel.f90` : prepares the linear set of equations for
!! connection of density and velocity by a linear correlation factor.
!! * `invermat.f90` : inverts the matrix BDERI + \ref ludcmp.f90 +  \ref lubksb.f90 + \ref dsinv.
!! * `perturb.f90` : calculates the perturbation of the parameters from the inversion of the matrix (BDERI).
!! * `resol.f90` : calculates the resolution terms for BDERI matrix.
!!
!!\subsection om 3.5 Others
!! * `statcoord.f90` : called by \ref vread and converts latitude and longitude coordinates
!!  of the stations into kilometric ones.
!! * \ref delaz : called by `statcoord.f90` and calculates the geocentric distance and azimuths.
!! * \ref cvrtop : called by \ref delaz and converts from rectangular to polar coordinates.
!! * `adnoise.f90` : called by `frechet.f90` and `calgra.f90` and finds gaussian random numbers
!!  and adds it to __synthetic__ datum.
!! * \ref ran1 : called by `adnoise.f90` subroutine returns a random number between 0 and 1. 
!! * `zeroav.f90` : reduces the mean of velocity and density variations for each layer to 0.
!! * `output.f90` : write several output files.
!!
!! \section Compilation
!!  - make costo (or make all) : compile all the programs and create the executable costo.
!!  - make clean : remove all mod_*.mod and executable costo.
!!


!
!=====================================================================
!========================OTHER FILE NEEDED============================
!=====================================================================
!
! RAYTRAC.PARAM  : Input parameters for 3D raytracing in velocity model
! -------------    Read in INFO3D
!
! DUMMY:   step length along raypath 
!          ex: 3
! DUMMY:   step length for velocity model partial derivatives
!          ex: 3
! DUMMY:   nb of harmonics, amplitude, amplitude ratio, cutoff, 
!          nb of sweep, min. time diff
!          ex: 9  1  1  0.0001 20 0.001
! DUMMY:   Three dimensional raytracing? [yes=1,no=0]
!          ex: 1
! DUMMY:   Extented output for 3D raytracing? [yes=1,no=0]
!          ex: 0
! DUMMY:   Minimum number of rays to take into account of the nodes
!          ex: 5
!=====================================================================
!=====================================================================
!
! PROGRAM FOR THE JOINT INVERSION OF GRAVITY, VELOCITY AND GRAVITY GRADIENT DATA

!=====================================================================
!=====================================================================
!  INCLUDE OF MODULE FOR VARIABLES DECLARATION
!=====================================================================
      INCLUDE 'DECLA/MOD_unit.f'
      INCLUDE 'DECLA/MOD_layer.f'
      INCLUDE 'DECLA/MOD_inv.f'
      INCLUDE 'DECLA/MOD_delim.f'
      INCLUDE 'DECLA/MOD_vdata.f'
      INCLUDE 'DECLA/MOD_iloc.f'

!=====================================================================
     !INCLUDE 'sblas/mkl_spblas.f90'
!=====================================================================
!  INCLUDE OF MODULE FOR EXPLICIT INTERFACE (SUBROUTINES)
!=====================================================================
      INCLUDE 'INTERF/MOD_readmodd.f'
      INCLUDE 'INTERF/MOD_readmodv.f'
      INCLUDE 'INTERF/MOD_readbco.f'
      INCLUDE 'INTERF/MOD_rdata.f'
      INCLUDE 'INTERF/MOD_bldmap.f'
      INCLUDE 'INTERF/MOD_vcoefs.f'
      INCLUDE 'INTERF/MOD_frechet.f'
      INCLUDE 'INTERF/MOD_caldt.f'
      INCLUDE 'INTERF/MOD_findnod.f'
      INCLUDE 'INTERF/MOD_zeroav.f'
      INCLUDE 'INTERF/MOD_ddiffer.f'
      INCLUDE 'INTERF/MOD_ddiffergr.f'      
      INCLUDE 'INTERF/MOD_vdiffer.f'
      INCLUDE 'INTERF/MOD_bdiffer.f'
      INCLUDE 'INTERF/MOD_calgra.f'
      INCLUDE 'INTERF/MOD_calgradio.f'  
      INCLUDE 'INTERF/MOD_contbod.f'
      INCLUDE 'INTERF/MOD_contnod.f'
      INCLUDE 'INTERF/MOD_bcalc.f'
      INCLUDE 'INTERF/MOD_dsmooth.f'
      INCLUDE 'INTERF/MOD_vsmooth.f'
      INCLUDE 'INTERF/MOD_bldmat.f90'
      INCLUDE 'INTERF/MOD_bldmat_modf.f90'
      INCLUDE 'INTERF/MOD_invermat.f'
      INCLUDE 'INTERF/MOD_resol.f'
      INCLUDE 'INTERF/MOD_perturb.f'
      INCLUDE 'INTERF/MOD_output.f'
      INCLUDE 'INTERF/MOD_timecal.f'
!=====================================================================
!  BEGINNING OF MAIN PROGRAM
!=====================================================================
       
      PROGRAM joint_inv

!=====================================================================
!  DECLARATION OF COMMUN MODULES and EXPLICIT INTERFACES
!=====================================================================


      USE MOD_unit
      USE MOD_layer
      USE MOD_size
      USE MOD_inv
      USE MOD_delim
      USE MOD_vdata
      USE MOD_iloc

      USE MOD_readmodd
      USE MOD_readmodv
      USE MOD_readbco
      USE MOD_rdata
      USE MOD_bldmap
      USE MOD_vcoefs
      USE MOD_frechet
      USE MOD_caldt
      USE MOD_findnod
      USE MOD_zeroav
      USE MOD_ddiffer
      USE MOD_vdiffer
      USE MOD_bdiffer
      USE MOD_ddiffergr  
      USE MOD_calgra
      USE MOD_calgradio 
      USE MOD_contbod
      USE MOD_contnod
      USE MOD_bcalc
      USE MOD_dsmooth
      USE MOD_vsmooth
      USE MOD_bldmat
      !USE MOD_bldmat_modf
      USE MOD_invermat
      USE MOD_resol
      USE MOD_perturb
      USE MOD_output
      USE MOD_timecal


      use ISO_Fortran_env
      USE MKL_SPBLAS

!=====================================================================
!  VARIABLES DECLARATION
!=====================================================================
      IMPLICIT NONE

      !INCLUDE 'mkl.fi'

      logical                                 :: temp1
	
      character(len=1)                           :: answer
      character(len=80)                          :: dummy,TITLE
      character(len=80),DIMENSION(:),ALLOCATABLE :: dfile

      integer                                 :: ier,i,j,npar1,nfil,temp
      integer                                 :: iter,iiter,ismooth,roexp
      integer                                 :: noised,noisev,noisegr,regul
      integer                                 :: nnod,npar,nbod,ndat,maxpar
      integer                                 :: nrp,invnod,NL,NC
      integer                                 :: constrain,m,nunode,velnod

      integer,DIMENSION(:),ALLOCATABLE        :: nlr
      integer,DIMENSION(:),ALLOCATABLE        :: codata,modvar
      integer,DIMENSION(:),ALLOCATABLE        :: x_block,y_block
      integer,DIMENSION(:),ALLOCATABLE        :: npts,ibove,ilay
      integer,DIMENSION(:),ALLOCATABLE        :: iside,ivside
      integer,DIMENSION(:,:),ALLOCATABLE      :: jside,jvside
      integer,DIMENSION(8)                    :: time_ori,time_blt1,time_blt2
      integer,DIMENSION(8)                    :: time_inv,time_fin
        
      integer,DIMENSION(:),pointer            :: ieq,ist

      real(kind=8)                            :: epsi,epsil,dtot0,dtot,lambda
      real(kind=8)                            :: thrhold
      real(kind=8)                            :: zroref
      real(kind=8)                            :: signoisd,signoisv,signoisgr
      real(kind=8)                            :: varv,varg
      real(kind=8),DIMENSION(6) 			  :: vargr
      real(kind=8),DIMENSION(3)               :: smooth
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: rmsg,rmst,errclt,errclg
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: rmsgr,errclgr
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: binit,varb,ddvr,drmax
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: ddtot,vdtot,bdtot
      real(kind=8),DIMENSION(:),ALLOCATABLE :: grdtot
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: rovar,velvar
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: rtvar
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: amask,dval
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: xb,yb,zb
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: par,par0,varpar
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: ddep
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: dxcor,dycor
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: vxnodes,vynodes,&
                                                 vznodes,v0
        
      real(kind=8),DIMENSION(:),ALLOCATABLE   :: velco,diff,h1
      real(kind=8),DIMENSION(:,:),ALLOCATABLE :: G,aderi,bderi
      real(kind=8),DIMENSION(:,:,:,:),ALLOCATABLE :: vels,vvari

      real(kind=8),DIMENSION(:),pointer       :: FX,FY,FZ,punvar
      real(kind=8),DIMENSION(:),pointer       :: usedata,caldata
      real(kind=8),DIMENSION(:),pointer       :: rayparm,bazin,weight
      real(kind=8),DIMENSION(:,:),pointer     :: stc,ttt

      integer                                 :: ii,ipar,k,l  
      real(kind=8)                            :: x,y,z  
      character*3 						          :: scount     
      character*20 						          :: file_name  
      real(kind=8)                            :: dvel     
      integer                                 :: ncomp     
      
   !  Aderi sparse variables
      
      real(kind=8), DIMENSION(:), ALLOCATABLE  :: val_ad_s
      integer, DIMENSION(:), ALLOCATABLE       :: columnAD, rowAD
      integer                                  :: nnz,idx_sp_A

      type(sparse_matrix_t)                     :: AtCA
      type(sparse_matrix_t)                     :: b_sp
      integer                                   :: status_mkl

      integer                                   :: memory_usage_bytes_Aderi
      integer                                   :: memory_usage_bytes_A_s
      integer                                   :: memory_usage_bytes_col_s
      integer                                   :: memory_usage_bytes_row_s
      

!=====================================================================
! Calculation of the date and time of the beginning of the run
!=====================================================================
      CALL DATE_AND_TIME(VALUES=time_ori)
!=====================================================================
! initialization of some arrays
!=====================================================================
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '    **************************************************************'
      write(*,*) '    ********** JOINT INVERSION OF TOMO / GRAVI / GRADIO **********'
      write(*,*) '    **************************************************************'
      write(*,*) ' '
      write(*,*) ' '


      npar1 = 0
      dtot0 = 0.d0
!=====================================================================
! Read the input files of the inversion parameters
!=====================================================================
      write(*,*)''
      write(*,*)'********************************************************'
      write(*,*)'             READING THE PARAMETER.INP FILE'
      write(*,*)'********************************************************'

      open(inpar,file='parameter.inp',status='OLD',iostat=ier)
      if(ier.ne.0) then
         write(*,*)''
         write(*,*)'Error in opening file: PARAMETER.INP, logical unit ', &
              inpar
         STOP 'in MAIN!'
      end if

      read(inpar,'(a80)') dummy
      read(inpar,'(a80)') TITLE

      read(inpar,'(/a80)') dummy
      read(inpar,'(l2,5x,l2,5x,l2,5x,l2)') INV,INVD,INVV,INVGR

      if(INVD.or.INVGR) then
        npar1=npar1+1
      end if
      if(INVV) then
        npar1=npar1+1
      end if
      if(npar1.gt.1) then
        npar1=npar1+1
      end if

      if((.not.INV) .and. ((.not.INVD) .and. (.not.INVV) .and. (.not.INVGR))) then
         write(*,*)''
         write(*,*)'Error 1 in consistency between inversion modes...'
         STOP 'in MAIN. Error in PARAMETER.INP file'
      end if
      if(INV.and.((.not.INVD).and.(.not.INVV).and.(.not.INVGR))) then
         write(*,*)''
         write(*,*)'Error 2 in consistency between inversion modes...'
         STOP 'in MAIN. Error in PARAMETER.INP file'
      end if
      
      read(inpar,'(/a80)') dummy
      read(inpar,*) iter,epsi
      epsil=epsi/100.

      if((.not.INV) .and. (iter.ne.1)) then
         write(*,*)''
         write(*,*)'WARNING: Direct Problem is asked, but there is more '
         write(*,*)'than 1 iteration in PARAMETER.INP.'
         write(*,*)''
         write(*,*)'Would you like to go on however? [Y/N]'
!         read(*,*) answer
         answer='Y'
         if(answer.eq.'N' .or. answer.eq.'n') then
            STOP 'You have to change ITER in PARAMETER.INP file'
         else if(answer.eq.'Y' .or. answer.eq.'y') then
            write(*,*)'Ok, we go on...'
         end if
      end if
!=====================================================================
! Allocation of arrays of iter dimension
!=====================================================================
      ALLOCATE (ddtot(iter+1,4))
      ALLOCATE (vdtot(iter+1,4))
      ALLOCATE (grdtot(iter+1))
      ALLOCATE (bdtot(iter+1,3))
      ALLOCATE (rovar(iter+1,3))
      ALLOCATE (velvar(iter+1,2))
      ALLOCATE (rmsg(iter+1))
      ALLOCATE (rmst(iter+1))
      ALLOCATE (rmsgr(iter+1,7))
      ALLOCATE (errclgr(iter+1,7))
      ALLOCATE (errclt(iter+1))
      ALLOCATE (errclg(iter+1))

      rmsg(:) = 0.d0
      rmst(:) = 0.d0
      rmsgr(:,:) = 0.d0

      read(inpar,'(/a80)') dummy
      read(inpar,*) thrhold
      if(thrhold.gt.1.) then
         write(*,*)''
         write(*,*)'THRHOLD value should be less than 1'
         STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if

      read(inpar,'(/a80)') dummy
      read(inpar,*) ismooth,smooth(1),smooth(2)
      smooth(3)=0.D0

      read(inpar,'(/a80)') dummy
      read(inpar,*) regul,lambda
      if(lambda.lt.0) then
        write(*,*)''
        write(*,*)'LAMBDA value should be greater than 1'
        STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if
      if(regul.ne.0 .and. regul.ne.1) then
        write(*,*)''
        write(*,*)'REGUL value should be 0 or 1'
        STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if

      read(inpar,'(/a80)') dummy
      read(inpar,*) noised,noisev,noisegr,signoisd,signoisv,signoisgr
      if(noised.ne.0 .and. noised.ne.1) then
         write(*,*)''
         write(*,*)'NOISED value should be 0 or 1'
         STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if
      if(noisev.ne.0 .and. noisev.ne.1) then
         write(*,*)''
         write(*,*)'NOISEV value should be 0 or 1'
         STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if
      if(noisegr.ne.0 .and. noisegr.ne.1) then
         write(*,*)''
         write(*,*)'NOISEGR value should be 0 or 1'
         STOP 'in MAIN! Problem in PARAMETER.INP file'
      end if      

      read(inpar,'(/a80)') dummy
      read(inpar,*) zroref,roexp
      read(inpar,'(/a80)') dummy
      read(inpar,*) nlayer

!=====================================================================
! Allocation of arrays of nlayer dimension
!=====================================================================
      ALLOCATE (nlr(nlayer))
      ALLOCATE (binit(nlayer))
      ALLOCATE (varb(nlayer))
      ALLOCATE (ddvr(nlayer))
      ALLOCATE (drmax(nlayer))

      read(inpar,'(/a80)') dummy
      do i=1,nlayer
         read(inpar,*) nlr(i),binit(i),varb(i),ddvr(i),drmax(i)
      end do
        
      read(inpar,'(/a80)') dummy
      read(inpar,*) densmod
      read(inpar,'(/a80)') dummy
      read(inpar,*) velmod

      read(inpar,'(/a80)') dummy
      read(inpar,*) nfil

!=====================================================================
! Allocation of arrays of nfil dimension
!=====================================================================
      ALLOCATE (codata(nfil))
      ALLOCATE (dfile(nfil))
      ALLOCATE (modvar(nfil))
      ALLOCATE (dval(nfil))
      ALLOCATE (amask(nfil))
      ALLOCATE (npts(nfil))

	   if(INVGR) then
        ALLOCATE (rtvar(nfil+5))
      elseif(.not.INVGR) then
        ALLOCATE (rtvar(nfil))
  	   endif

      do i=1,nfil
         read(inpar,'(/a80)') dummy
         read(inpar,*) codata(i)
         read(inpar,'(/a80)') dummy
         read(inpar,'(a80)') dfile(i)
         read(inpar,'(/a80)') dummy
         if(INVGR) then
         	if ( i == nfil ) then
                read(inpar,*) modvar(i),rtvar(i),rtvar(i+1),rtvar(i+2),rtvar(i+3),&
                rtvar(i+4),rtvar(i+5)
         	else 
                read(inpar,*) modvar(i),rtvar(i)
         	endif
     	   else
            read(inpar,*) modvar(i),rtvar(i)
      	endif
         if(modvar(i).ne.0 .and. modvar(i).ne.1) then
            write(*,*)''
            write(*,*)'Allowed values for MODVAR is 0 or 1.'
            STOP 'in MAIN. Error in PARAMETER.INP file'
         end if
      end do

!=====================================================================
! Writting some output comments in inout file "parameter.out"
!=====================================================================

      open(inout,file='parameter.out',status='REPLACE',iostat=ier)
      if(ier.ne.0) then
         write(*,*) 'Error in opening file: parameter.out, &
              & logical unit ',inout
         STOP 'in MAIN.'
      end if

      write(inout,*)' '
      write(inout,*)' '
      write(inout,*)'*********************************************'
      write(inout,*)'********** JOINT INVERSION OF DATA **********'
      write(inout,*)'*********************************************'
      write(inout,*)' '
      write(inout,*)' '
      write(inout,*)' '
      write(inout,*)'CHECKING PARAMETERS OF THE INVERSION (MAIN)'
      write(inout,*)' '
      write(inout,'(4x,''Title:'',1x,a80/)') TITLE

! écritures des output si c'est juste probleme DIRECT    	
      if(.not.INV) then
         write(*,'(4x,''Direct Problem:'')')
         write(inout,'(4x,''Direct Problem:'')')
         if(INVV) then
            write(*,'(7x,''-Synthetic delay times calculated'')')
            write(inout,'(7x,''-Synthetic delay times calculated'')')
            if(noisev.eq.1) then
               write(*,'(7x,''-Adding gaussian noise to synthetic (sigma='',&
                    &f6.3,'')'')') signoisv
               write(inout,'(7x,''-Adding gaussian noise to synthetic &
                    &(sigma='',f6.3,'')'')') signoisv
            else
               write(*,'(7x,''-No gaussian random noise to synthetic'')')
               write(inout,'(7x,''-No gaussian random noise to synthetic'')')
            end if
         end if
         if(INVD) then
            write(*,'(7x,''-Synthetic gravity data calculated'')')
            write(inout,'(7x,''-Synthetic gravity data calculated'')')
            if(noised.eq.1) then
               write(*,'(7x,''-Adding gaussian noise to synthetic (sigma='',&
                    &f6.3,'')'')') signoisd
               write(inout,'(7x,''-Adding gaussian noise to synthetic &
                    &(sigma='',f6.3,'')'')') signoisd
            else
               write(*,'(7x,''-No gaussian random noise to synthetic'')')
               write(inout,'(7x,''-No gaussian random noise to synthetic'')')
            end if
         end if
         if(INVGR) then
            write(*,'(7x,''-Synthetic gravity gradients calculated'')')
            write(inout,'(7x,''-Synthetic gravity gradients calculated'')')
            if(noised.eq.1) then
               write(*,'(7x,''-Adding gaussian noise to synthetic (sigma='',&
                    &f6.3,'')'')') signoisgr
               write(inout,'(7x,''-Adding gaussian noise to synthetic &
                    &(sigma='',f6.3,'')'')') signoisgr
            else
               write(*,'(7x,''-No gaussian random noise to synthetic'')')
               write(inout,'(7x,''-No gaussian random noise to synthetic'')')
            end if
         end if

      else

! écritures des output si IL Y A INVERSION      	
         write(*,'(4x,''Inversion Problem:'')')
         write(inout,'(4x,''Inversion Problem:'')')
         if(INVV.and.(.not.INVD).and.(.not.INVGR)) then
            write(*,'(7x,''-Inversion of teleseismic delay times'')')
            write(inout,'(7x,''-Inversion of teleseismic delay times'')')
         end if
         if(INVD.and.(.not.INVV).and.(.not.INVGR)) then
            write(*,'(7x,''-Inversion of gravity'')')
            write(inout,'(7x,''-Inversion of gravity'')')
         end if
         if(INVV.and.INVD.and.(.not.INVGR)) then
            write(*,'(7x,''-Joint inversion of velocity/gravity'')')
            write(inout,'(7x,''-Joint Inversion of velocity/gravity'')')
         end if
         if(INVGR.and.(.not.INVV) .and.(.not.INVD)) then
            write(*,'(7x,''-Inversion of gravity gradients'')')
            write(inout,'(7x,''-Inversion of gravity gradients'')')
         end if
         if(INVGR.and.INVV.and.(.not.INVD)) then
            write(*,'(7x,''-Joint inversion of gravity gradients/velocity'')')
            write(inout,'(7x,''-Joint inversion of gravity gradients/velocity'')')
         end if
         if(INVGR.and.INVD.and.(.not.INVV)) then
            write(*,'(7x,''-Joint inversion of gravity gradients/gravity'')')
            write(inout,'(7x,''-Joint inversion of gravity gradients/gravity'')')
         end if
         if(INVGR.and.INVD.and.INVV) then
            write(*,'(7x,''-Joint inversion of gravity gradients/gravity/velocity'')')
            write(inout,'(7x,''-Joint inversion of gravity gradients/gravity/velocity'')')
         end if

      end if

      write(inout,*)''
      write(inout,'(4x,''Number of parameters to be inverted &
           &per block:'',1x,i5/)') npar1
           
      select case (ismooth)
             case (0)
                write(inout,'(4x,''No smoothest model'')')
             case (1)
                write(inout,'(4X,''Smoothest model for all blocks'')')
                write(inout,'(7x,''-Smoothing factors:'',1p,3E10.2)')&
                     (smooth(i),i=1,3)
             case default
                write(*,*)''
                write(*,*)'Ismooth value must be 0 or 1'
                STOP 'in MAIN! Error in PARAMETER.INP file.'
      end select

      select case (regul)
	     case (0)
            write(inout,'(4x,''No regularization of the matrix'')')
	     case (1)
            write(inout,'(4X,''Regularization of the matrix'')')
            write(inout,'(7x,''-regularization factor:'',1p,E10.2)')&
               lambda
	     case default
            write(*,*)''
            write(*,*)'Regul value must be 0 or 1'
            STOP 'in MAIN! Error in PARAMETER.INP file.'
      end select

      write(inout,'(4x,''Description of the'',i4,'' layers'')') nlayer
      write(inout,'(4x,''NL   B-init    VAR.B0    VAR.VR    DROMAX'')')
      write(inout,'(4x,43(''-''))')
      do i=1,nlayer
        write(inout,'(i5,1x,4e10.2)') nlr(i),binit(i),varb(i),&
         ddvr(i),drmax(i)
      end do

!=====================================================================
! Closing the input parameter file parameter.inp
!=====================================================================
      close(inpar,iostat=ier)
      if(ier.ne.0) then
         write(*,*)'Error in closing the file PARAMETER.INP, logical &
                   &unit ',inpar
         STOP 'in MAIN.'
      end if

!=====================================================================
! Read the maximum size of the array model (density and velocity)
! and allocate the memory needed for those arrays
!=====================================================================

      CALL READ_LAYER

      if(INV.and.INVV.and.INVD.and.INVGR) then
         maxpar = nxnode*nynode*nznode + nlayer*nxbloc*nybloc + nlayer
      elseif(INV.and.INVD.and.INVGR) then
         maxpar = nlayer*nxbloc*nybloc
      elseif(INV.and.INVV.and.INVGR) then
         maxpar = nxnode*nynode*nznode + nlayer*nxbloc*nybloc + nlayer
      elseif(INV.and.INVV.and.INVD) then
         maxpar = nxnode*nynode*nznode + nlayer*nxbloc*nybloc + nlayer
      elseif(INV.and.INVGR) then
         maxpar = nlayer*nxbloc*nybloc
      elseif(INV.and.INVD) then
         maxpar = nlayer*nxbloc*nybloc
      elseif(INV.and.INVV) then
         maxpar = nxnode*nynode*nznode

      elseif((.not.INV).and.INVGR) then
         maxpar = nlayer*nxbloc*nybloc
      elseif((.not.INV).and.INVD) then
         maxpar = nlayer*nxbloc*nybloc
      elseif((.not.INV).and.INVV) then
         maxpar = nxnode*nynode*nznode
      end if

      write(*,*)''
      write(*,*)'Maximum of parameters for memory allocation:',maxpar
      write(inout,*)''
      write(inout,*)'Maximum of parameters for memory allocation:',maxpar
      ALLOCATE(par(maxpar))
      ALLOCATE(par0(maxpar))
      ALLOCATE(varpar(maxpar))
      par(:)=0.d0
      par0(:)=0.d0
      varpar(:)=0.d0

      if(INVD.or.INVGR) then
         ALLOCATE(x_block(nlayer))
         ALLOCATE(y_block(nlayer))
         ALLOCATE(ddep(nlayer+1))
         ALLOCATE(dxcor(nxbloc+1,nlayer))
         ALLOCATE(dycor(nybloc+1,nlayer))
         ALLOCATE(ilay(nlayer*nxbloc*nybloc))
         ALLOCATE(xb(nlayer*nxbloc*nybloc,2))
         ALLOCATE(yb(nlayer*nxbloc*nybloc,2))
         ALLOCATE(zb(nlayer*nxbloc*nybloc,2))
      end if

      if(INVV) then
         ALLOCATE(vxnodes(nxnode))
         ALLOCATE(vynodes(nynode))
         ALLOCATE(vznodes(nznode))
         ALLOCATE(v0(nznode))
         ALLOCATE(vels(nxnode,nynode,nznode,iter+1))
         ALLOCATE(vvari(nxnode,nynode,nznode,iter))
      end if

!=====================================================================
! Read the input files for the model parameters of density and
! velocity and initial values of B-coeff
!=====================================================================
      ibegin(:)=0
      iend(:)=0

      npar=0
      nbod=0
      nnod=0
!=====================================================================
! If inversion of data, need to know the model geometry (parameters cald_sp(iend(2)))
! location).
! If direct problem, need to know the model geometry to pass through
! IBEGIN and IEND arrays are initialized here.
!=====================================================================
      write(*,*)''
      write(*,*)'********************************************************'
      write(*,*)'             READING THE MODEL FILE(S)'
      write(*,*)'********************************************************'

      if(INVD.or.INVGR) then
         print*,"début lecture"
         CALL READ_MODD(zroref,roexp,ilay,xb,yb,zb,par,varpar,npar,nbod,&
                        x_block,y_block,ddep,dxcor,dycor)
         print*,"fin lecture"
      else
         write(inout,*)''
         write(inout,*)'NO DENSITY MODEL READ'
         write(inout,*)'    ibegin(1)=0'
         write(inout,*)'    iend(1)=0'
      end if

      if(INVV) then
         CALL READ_MODV(v0,par,varpar,nnod,nbod,npar,vels,vvari,vxnodes,&
                        vynodes,vznodes)
      else
         write(inout,*)''
         write(inout,*)'NO VELOCITY MODEL READ'
         write(inout,*)'    ibegin(2)=0'
         write(inout,*)'    iend(2)=0'
      end if


      !if((INV.and.INVV.and.INVD).or.(INV.and.INVV.and.INVGR)) then
      if(INV.and.INVV.and.(INVD .or. INVGR)) then
         CALL READBCO(binit,varb,par,varpar,npar)
      else
         write(inout,*)''
         write(inout,*)'NO LINEAR RELATION BETWEEN VELOCITY AND DENSITY'
         write(inout,*)'    ibegin(3)=0'
         write(inout,*)'    iend(3)=0'
      end if

!=====================================================================
! Checking the good storage of the model parameters
!=====================================================================

      write(inout,*)' '
      write(inout,*)'TOTAL NUMBER OF BODIES:',npar
      write(inout,*)' '
     
!=====================================================================
! If inversion of data, read the datapoint location and data from the
! data files, and substract the mean value of the data
! If Direct Problem, we need the fieldpoint location to know where
! to calculate synthetic data
! JBEGIN, JEND arrays are initialized here
!=====================================================================
      CALL RDATA(nfil,codata,dfile,modvar,rtvar,&
                 FX,FY,FZ,DVAL,AMASK,usedata,npts,rayparm,bazin,&
                 weight,ndat,ist,stc,ieq,ttt,punvar,varv,varg,vargr,ncomp)

!=====================================================================
! Read information about the velocity 3D ray tracing model
!=====================================================================
      if(INVV) then
         CALL INFO3D
!=====================================================================
! Building a Map with index of the nodes for 3D raytracing and initial
! velocity models
!=====================================================================

         CALL BLDMAP(vxnodes,vynodes,vznodes)
      end if

!=====================================================================
! END OF INPUTS MODEL AND DATA
! Checking the correct storage of the data
!=====================================================================
      write(inout,*)' '
      write(inout,*)'CHECKING THE GOOD STORAGE OF DATA AND PARAMETERS'
      write(inout,*)' '
      write(inout,*)'DENSITIES   FROM COL.',ibegin(1),' TO ',iend(1)
      write(inout,*)'VELOCITIES  FROM COL.',ibegin(2),' TO ',iend(2)
      write(inout,*)'dV/dRHO, V0 FROM COL.',ibegin(3),' TO ',iend(3)
      write(inout,*)'NUMBER OF PARAMETERS PER BODY:............',npar1
      write(inout,*)'THEORIC NUMBER OF PARAMETERS TO INVERT:...',npar
      write(inout,*)'GRAVITY   DATA FROM ROW ',jbegin(1),' TO ',jend(1)
      write(inout,*)'DELAYTIME DATA FROM ROW ',jbegin(2),' TO ',jend(2)
      write(inout,*)'FTG DATA FROM ROW ',jbegin(3),' TO ',jend(3)     
      write(inout,*)'B-COEFF DATA FROM ROW ',jbegin(4),' TO ',jend(4) 
      write(inout,*)'TOTAL NUMBER OF DATA:.....................',ndat   

!=========================
      write(*,*)' '
      write(*,*)'CHECKING THE GOOD STORAGE OF DATA AND PARAMETERS'
      write(*,*)' '
      write(*,*)'DENSITIES   FROM COL.',ibegin(1),' TO ',iend(1)
      write(*,*)'VELOCITIES  FROM COL.',ibegin(2),' TO ',iend(2)
      write(*,*)'dV/dRHO, V0 FROM COL.',ibegin(3),' TO ',iend(3)
      write(*,*)'NUMBER OF PARAMETERS PER BODY:............',npar1
      write(*,*)'THEORIC NUMBER OF PARAMETERS TO INVERT:...',npar
      write(*,*)'GRAVITY   DATA FROM ROW ',jbegin(1),' TO ',jend(1)
      write(*,*)'DELAYTIME DATA FROM ROW ',jbegin(2),' TO ',jend(2)
      write(*,*)'FTG DATA FROM ROW ',jbegin(3),' TO ',jend(3)     
      write(*,*)'B-COEFF DATA FROM ROW ',jbegin(4),' TO ',jend(4) 
      write(*,*)'TOTAL NUMBER OF DATA:.....................',ndat        

!=====================================================================
! Allocate memory for the ADERI array(nbdata,nbparam) and initialization
!=====================================================================
      !Aderi blocks dim
      !A_rho_dim = iend(1)*jend(1) 
      !A_v_dim = (iend(2)-ibegin(2))*(jend(2)-jbegin(2)) ->some elements migth be Zero
      !A_B_dim = (iend(3)-ibegin(3))*ndat
      
      nnz = iend(1)*jend(1) + (iend(2)-ibegin(2))*(jend(2)-jbegin(2)) +&
            (iend(1)*(jend(3)-jbegin(3)))

      write(*,*) 'value of nnz ', nnz, 'rho ', iend(1)*jend(1) ,&
                  'v ',(iend(2)-ibegin(2))*(jend(2)-jbegin(2)) , &
                  'gra ',  (iend(1)*(jend(3)-jbegin(3)))

      ALLOCATE (val_ad_s(nnz))
      ALLOCATE (columnAD(nnz))
      ALLOCATE (rowAD(nnz))

      ALLOCATE (aderi(ndat,npar))
      ALLOCATE (bderi(npar,npar))
      ALLOCATE (h1(npar))
      ALLOCATE (caldata(ndat))
      ALLOCATE (diff(ndat))
      ALLOCATE (ibove(npar))
      ALLOCATE (iside(nbod))
      ALLOCATE (jside(nbod,4))
      ALLOCATE (ivside(nxnode*nynode*(nznode-1)))
      ALLOCATE (jvside(nxnode*nynode*(nznode-1),4))
      caldata(:) = 0.d0
      aderi(:,:) = 0.d0
      bderi(:,:) = 0.d0
      diff(:) =0.d0

      val_ad_s(:) = 0.d0
      columnAD(:) = 0
      rowAD(:)    = 0
      idx_sp_A = iend(1)*jend(1) 

      memory_usage_bytes_Aderi = npar*ndat/(1.d-9)! Size in bytes
      write (*,*) "Approx. memory usage for 'ADERI': ", memory_usage_bytes_Aderi, " bytes"
      write(inout,*) "Approx. memory usage for 'ADERI': ", memory_usage_bytes_Aderi, " bytes"


      memory_usage_bytes_A_s= size(val_ad_s) * storage_size(val_ad_s) / 8 ! Size in bytes
      write (*,*) "Approx. memory usage for 'ADERI sparse' array val: ", memory_usage_bytes_A_s , " bytes"
      write(inout,*) "Approx. memory usage for 'ADERI sparse' array val: ", memory_usage_bytes_A_s , " bytes"

      write(*,*) 'nnz ', nnz ,  'size As', size(val_ad_s) ,'size col_a ', size(columnAD) , &
               'storage_size col_s ' , storage_size(columnAD), 'COO = ', nnz*16/1.d-9
         
      write(inout,*)  'nnz ', nnz ,  'size As', size(val_ad_s) ,'size col_a ', size(columnAD) , &
      'storage_size col_s ' , storage_size(columnAD)

      memory_usage_bytes_col_s = size(columnAD) * storage_size(columnAD) / 8 ! Size in bytes
      write (*,*) "Approx. memory usage for 'ADERI sparse' array col: ", memory_usage_bytes_col_s, " bytes"




     ! memory_usage_bytes_col_s = size(columnAD) * storage_size(columnAD) 
     ! write (*,*) "Approx. memory usage for 'ADERI sparse' col array: ", memory_usage_bytes_col_s , " bits"


      !memory_usage_bytes_row_s = size(rowAD) * storage_size(rowAD) 
      !write (*,*) "Approx. memory usage for 'ADERI sparse' row array: ", memory_usage_bytes_row_s , " bits"
  

      if(INVV) then
         ALLOCATE (velco(8*nxnode*nynode*nznode))
         ALLOCATE (G(n_data,nnod))
         velco(:) = 0.d0
         G(:,:) = 0.d0
      end if
        
!=====================================================================
! If inversion of delay-times or calculation of delay-times, we have
! to compute parameters for the 3D raytracing and the density nodes
! that are constrained by velocity nodes
!=====================================================================
      iiter = 1
      if(INVV) then
!=====================================================================
! compute the coefficients for the linear velocity interpolation of 3D RT
!=====================================================================
         CALL VCOEFS(velco,iiter,vels,vxnodes,vynodes,vznodes)
        
!=====================================================================
! Making the frechet matrix to see which nodes are hitten by the rays
! invnod = number of inverted velocity nodes
!
! In case of direct problem, return the travel times in TIME.PRED file
! and stop the program here.
!=====================================================================
         CALL FRECHET(vxnodes,vynodes,vznodes,bazin,rayparm,&
                      ieq,stc,ist,G,nrp,velco,invnod,signoisv,noisev)
!=====================================================================
! Construct the file with hit nodes number IBOVE (CALDT)
! and find the density blocks constrained (stored in file bloc.nod, FINDNOD)
!=====================================================================
         !CALL CALDT(ibove,invnod,caldata,vels,ieq,aderi,par)
         CALL CALDT(ibove,invnod,caldata,vels,ieq,aderi,par,val_ad_s,&
               columnAD,rowAD, idx_sp_A,nnz,ndat,npar)
         if(INVD.or.INVGR) then
            CALL FINDNOD(xb,yb,zb,vxnodes,vynodes,vznodes,nbod,ibove)

            write(inout,*)''
            write(inout,*)'AFTER CALL OF FINDNOD'

            open(blocnod,file='bloc.nod',status='OLD',iostat=ier)
            if(ier.ne.0) then
               write(*,*)'Error in opening the file BLOC.NOD, &
                    &logical unit ',blocnod
               STOP 'in MAIN.'
            end if

            constrain=0
            temp = 0
            do i=1,nbod
               read(blocnod,*) m,nunode
               if(nunode.ne.0) then
                  constrain=constrain+1
                  temp1=.false.
                  do j=1,nunode
                     read(blocnod,*) velnod
                     write(inout,*)'    Dens. block no ',i,' contains &
                          &vel. node no ',velnod
                     if(ibove(velnod).ne.0) then
                        temp1=.true.
                     end if
                  end do
                  if(temp1) temp=temp+1
               end if
            end do
            
            write(inout,*)''
            write(inout,*)'    Number of density blocks containing &
                 &velocity node(s): ',constrain
            write(inout,*)'    Number of real node constained density &
                 &blocks: ',temp

            close(blocnod,iostat=ier)
            if(ier.ne.0) then
               write(*,*)'Error in closing the file BLOC.NOD, logical &
                    &unit ',blocnod
               STOP 'in MAIN.'
            end if

         end if

      end if
!=====================================================================
! Put the initial parameters in array PAR0
! And get the first standard deviation by calling DIFFER (sdevg0,sdevt0)
! But does it really matter?????
!=====================================================================

      par0(:) = par(:)
      rovar(:,:) = 0.d0
      velvar(:,:) = 0.d0

      if(INV.and.INVD) then
         CALL DDIFFER(usedata,caldata,rmsg,iiter,par,par0,punvar,rovar,&
                     varpar,xb,yb,zb,ismooth,smooth,iside,jside,diff,ddtot,&
                     h1,bderi,errclg)
      end if
      if(INV.and.INVV) then
         CALL VDIFFER(usedata,caldata,rmst,iiter,par,par0,punvar,&
                     varpar,vxnodes,vynodes,ismooth,smooth,ivside,&
                     jvside,ibove,diff,vdtot,velvar,h1,bderi,errclt)
      end if

      if(INV.and.INVGR) then
         CALL DDIFFERGR(usedata,caldata,rmsgr,iiter,par,par0,punvar,rovar,&
                     varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
                     diff,ddtot,grdtot,h1,bderi,npts,errclgr,rtvar,ncomp)
      end if

      if(INV.and.INVV.and.INVD.and.(.not.INVGR)) then
        CALL BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,velco,&
                v0,nbod,bdtot,ibove)
      elseif(INV.and.INVV.and.INVGR.and.(.not.INVD)) then
        CALL BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,velco,&
                v0,nbod,bdtot,ibove)
      elseif(INV.and.INVV.and.INVD.and.INVGR) then
        CALL BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,velco,&
                v0,nbod,bdtot,ibove)
      end if

!=====================================================================
! Beginning of the calculations.......
! Either the inverse problem, or the forward calculation
!=====================================================================
      if(INV) then
         write(*,*)''
         write(*,*)'********************************************************'
         write(*,*)'   Starting the ',iter,' iterations for the inversion'
         write(*,*)'********************************************************'

         write(inout,*)''
         write(inout,*)'*************************************************&
              &*******'
         write(inout,*)'   Starting the ',iter,' iterations for the inversion'
         write(inout,*)'*************************************************&
              &*******'
         write(inout,*)''
      else
         write(*,*)''
         write(*,*)'********************************************************'
         write(*,*)'   Starting the forward calculation of synthetic data'
         write(*,*)'********************************************************'

         write(inout,*)''
         write(inout,*)'****************************************************&
              &****'
         write(inout,*)'   Starting the forward calculation of synthetic data'
         write(inout,*)'****************************************************&
              &****'
         write(inout,*)''
      end if

!=====================================================================
! Beginning of the loops over the iterations
! We go until iiter=iter+1, because we need to pass through the
! subroutine of zero average and differences calculation for the last
! iteration. It does not invert for iter+1, it exits the loop before.
!=====================================================================
      if((INV.and.INVV.and.INVD.and.(.not.INVGR)).or.(INV.and.INVV.and.INVGR.and.(.not.INVD)).or.&
          (INV.and.INVV.and.INVGR.and.INVD)) then
         open(correl,file='correl.out',status='REPLACE',position='APPEND',&
              iostat=ier)
         if(ier.ne.0) then
            write(*,*)''
            write(*,*)'Error in opening the file CORREL.OUT'
            STOP 'in MAIN'
         end if
         write(correl,*)'iter layer  density   velocity'

         open(bcoeff,file='bcoeff.res',status='REPLACE',position='APPEND',&
              iostat=ier)
         if(ier.ne.0) then
            write(*,*)''
            write(*,*)'Error in opening the file BCOEFF.OUT'
            STOP 'in MAIN'
         end if
      end if

      do iiter=1,iter+1

         write(*,*)''
         write(*,*)'***************'
         write(*,*)'ITERATION No ',iiter
         write(*,*)'***************'
         
         write(inout,*)''
         write(inout,*)'***************'
         write(inout,*)'ITERATION No ',iiter
         write(inout,*)'***************'


!=====================================================================
! Compute zero average of the model parameters for inversion hypothesis
!=====================================================================
         if(INV) then
!         	print*,"#######################################"
!         	print*,"#######################################"
!         	print*," ATTENTION ZEROV N EST PAS CALCULE CAR LE CODE BEUGGUAIT &
!         		SINON, SUREMENT UN PROBLEME AVEC LES FICHIERS INPUT"
!         	print*,"#######################################"
!         	print*,"#######################################"

            CALL ZEROAV(x_block,y_block,par,drmax,ibove)
         end if
!=====================================================================
! Forward calculation of gravity anomalies (CALGRA) with zero average.
! Synthetic data are stored in CALDATA array, frechet matrix in ADERI.
! If direct problem the program is stopped in CALGRA.
!
! BE CARREFUL...
! there is a problem with aderi(:,:)=0 and the first lines in CALGRA.
! I have to suppress aderi=0 or to suppress the first part of IF in
! CALGRA calculation... Have to test this.
! After test, it gives the same results...
!=====================================================================
         aderi(:,:) = 0.d0
         bderi(:,:) = 0.d0


         val_ad_s(:) = 0.d0
         columnAD(:) = 0
         rowAD(:)    = 0
         idx_sp_A = 0
         
         if(INVD) then
            CALL CALGRA(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
                        noised,signoisd,val_ad_s,columnAD,rowAD, idx_sp_A)
         end if

         if(INVGR) then
            CALL CALGRADIO(iiter,aderi,caldata,par,FX,FY,FZ,nbod,xb,yb,zb,&
                             noisegr,signoisgr,rtvar,ncomp,val_ad_s,columnAD,rowAD, idx_sp_A)
         end if

!=====================================================================
! Forward calculation of delay-times (CALDT) using a 3D raytracing
! synthetic data are stored in CALDATA array, frechet matrix in ADERI
! If direct problem the program is stopped in FRECHET.
!=====================================================================
         if(INVV) then
            CALL VCOEFS(velco,iiter,vels,vxnodes,vynodes,vznodes)
            CALL FRECHET(vxnodes,vynodes,vznodes,bazin,rayparm,ieq,&
                         stc,ist,G,nrp,velco,invnod,signoisv,noisev)
!=====================================================================
! Here should be the subroutine for the ray density calculation
! and the new nodes coordinates calculation:
!        if(move.eq.1) then
!           CALL RAYDENS()
!        end if
!=====================================================================
            CALL CALDT(ibove,invnod,caldata,vels,ieq,aderi,par,&
                        val_ad_s,columnAD,rowAD, idx_sp_A,nnz,ndat,npar)
         end if



!=====================================================================
! From here now, it only concerns inverse problem of data as the
! forward problems has stopped before in CALGRA and CALDT subroutines
!=====================================================================

!=====================================================================
! Calcul of smoothing constraint (SMODEL,VSMODEL)
!=====================================================================
         if(INV .and. ismooth.eq.1 .and. iiter.eq.1) then
            write(inout,*)''
            write(inout,*)'DETERMINATION OF CONTIGUOUS BODY, &
                 &ISMOOTH = 1, FIRST ITERATION'

            if(INVD.or.INVGR) CALL CONTBOD(xb,yb,nbod,ilay,iside,jside)
            if(INVV) CALL CONTNOD(ibove,ivside,jvside)
         end if
         
!=====================================================================
! In case of joint inversion, calculate the B-coeff and write the
! densities/velocities in correl.out file
!=====================================================================
         if((INV.and.INVD.and.INVV).or.(INV.and.INVGR.and.INVV)) then
!=====================================================================
!        if(move.eq.1) then
!           CALL FINDNOD()
!        end if
!=====================================================================
            CALL BCALC(nbod,ilay,par,ibove,iiter,thrhold)
         end if

!=====================================================================
! difference between obs. and calc. data, and actual and previous parameter
!=====================================================================
         if(INV.and.INVD) then
            CALL DDIFFER(usedata,caldata,rmsg,iiter,par,par0,punvar,rovar,&
                        varpar,xb,yb,zb,ismooth,smooth,iside,jside,diff,&
                        ddtot,h1,bderi,errclg)
         end if

         if(INV.and.INVGR) then
            CALL DDIFFERGR(usedata,caldata,rmsgr,iiter,par,par0,punvar,rovar,&
                     varpar,xb,yb,zb,ismooth,smooth,iside,jside,&
                     diff,ddtot,grdtot,h1,bderi,npts,errclgr,rtvar,ncomp)
         end if

         if(INV.and.INVV) then
            CALL VDIFFER(usedata,caldata,rmst,iiter,par,par0,punvar,&
                        varpar,vxnodes,vynodes,ismooth,smooth,ivside,&
                        jvside,ibove,diff,vdtot,velvar,h1,bderi,errclt)
         end if

         if((INV.and.INVD.and.INVV).or.(INV.and.INVGR.and.INVV)) then
            CALL BDIFFER(iiter,par,par0,varpar,ilay,xb,yb,zb,varb,&
                         velco,v0,nbod,bdtot,ibove)
         end if

         dtot = ddtot(iiter,4) + vdtot(iiter,4) + bdtot(iiter,3) + grdtot(iiter)

         ! enregistrement modèle densité
         if(INVD.or.INVGR) then
            WRITE(scount,'(i3)') iiter
            file_name  = scount//"_iter_density.res"
            open(unit = 444,file = file_name,status ='replace') 
                     ii = 0
                     do i=ibegin(1),iend(1)
                        ii = ii + 1
                        x=(xb(ii,1)+xb(ii,2))*0.5
                        y=(yb(ii,1)+yb(ii,2))*0.5
                        z=(zb(ii,1)+zb(ii,2))*0.5
                        write(444,'(4(f10.3,2x))') x,y,z,par(i)
                     end do
            close(444)
         endif

         ! enregistrement modèle vitesse
         if(INVV) then
            WRITE(scount,'(i3)') iiter
            file_name  = scount//"_iter_velocity.res"
            open(unit = 444,file = file_name,status ='replace') 
                     ipar = 0
                     do k=1,nznode-1
                        do j=1,nynode
                           do i=1,nxnode
                              ipar = ipar + 1
                              if(ibove(ipar).ne.0) then
                                 dvel=par(ibove(ipar))/vels(i,j,k,1)*100
                              else
                                 dvel=0.d0
                              end if
                              write(velout,'(4(f10.4,2x),10f10.3)') vxnodes(i),vynodes(j),&
                                    vznodes(k),dvel,(vels(i,j,k,l),l=1,iter+1)
                           end do
                        end do
                     end do
            close(444)
         endif

! enregistrement de la réponse du modèle

!=====================================================================
! If the difference between 2 iterations is small enough, stop the
! process and write the results
! If not, check if dtot is lower than epsi to stop the process.
! if the number of iteration is reached, exit the loop and write the
! results.
! dtot0 = 0 at first iteration
!=====================================================================
         if(DABS(dtot0-dtot).lt.epsil) then
            write(*,*)''
            write(*,'(''Mean difference between two iterations '',&
                 &f6.4,'' is < epsi/100.'')') DABS(dtot0-dtot)
            write(*,*)'I quit the iterations and write the results'

            write(inout,*)''
            write(inout,'(''Mean difference between two iterations '',&
                 &f6.4,'' is < epsi/100.'')') DABS(dtot0-dtot)
            write(inout,*)'I quit the iterations and write the results'
            exit
         else
            dtot0=dtot
            if(DABS(dtot).lt.epsi) then
               write(*,*)''
               write(*,'(''Mean difference between model and data '',&
                    &f6.4,'' is < epsi: '')') DABS(dtot)
               write(*,*)'I stop the iteration process and write the results'

               write(inout,*)''
               write(inout,'(''Mean difference between model and data '',&
                    &f6.4,'' is < epsi: '')') DABS(dtot)
               write(inout,*)'I stop the iteration process and write &
                    &the results'
               exit
            end if
         end if
         if(iiter.eq.(iter+1)) then
            write(*,*)''
            write(*,*)'End of the iterations (iter=',iter,')'
            write(*,*)'I stop the calculations and write the results'

            write(inout,*)''
            write(inout,*)'End of the iterations (iter=',iter,')'
            write(inout,*)'I stop the calculations and write the results'
            exit
         end if
!=====================================================================
! Computation of the matrix (At.Cd**-1.A + Cp**-1 - Cd**-1.Db - Cs**-1.Ds)
! to invert
! and calulation of the computing time
!=====================================================================
         CALL DATE_AND_TIME(VALUES=time_blt1)

         CALL BLDMAT(iiter,aderi,bderi,punvar,varpar,h1,diff,npar,&
                     npar1,ndat,smooth,iside,jside,ivside,jvside,&
                     xb,yb,vxnodes,vynodes,par,ibove,ismooth,ilay,&
                     ddvr,nbod,val_ad_s, rowAD, columnAD, nnz,&
                     AtCA,B_sp)
                     
         CALL DATE_AND_TIME(VALUES=time_blt2)
         CALL TIMECAL(time_blt1,time_blt2)
!=====================================================================
! inversion of the matrix to obtain modelparameters estimation
! and calulation of the computing time
!=====================================================================/..
         CALL INVERMAT(bderi,npar,iiter,regul,lambda)

         CALL DATE_AND_TIME(VALUES=time_inv)
         CALL TIMECAL(time_blt2,time_inv)
!=====================================================================
! Calcul of the new parameters in PAR and VARPAR array
! Calcul of the standard deviation
!=====================================================================
         CALL PERTURB(bderi,h1,npar,drmax,ilay,par,par0,vels,iiter,&
                      ibove,xb,yb,zb)

!=====================================================================
! calcul of the resolution... It seems correct
!=====================================================================
         if (iiter.gt.1) then
            write(*,*)''
            write(*,*)'RESOLUTION CALCULATION'

            write(inout,*)''
            write(inout,*)'RESOLUTION CALCULATION'

            CALL RESOL(bderi,aderi,npar,ndat,xb,yb,ilay,vxnodes,vynodes,ibove,punvar)
         end if
        
!=====================================================================
! End of loops over the iterations
!=====================================================================
      end do
!=====================================================================
! End of Calculation
! Write the outputs on files 
!=====================================================================

      CALL OUTPUT(iter,ddtot,grdtot,vdtot,bdtot,velvar,rovar,par,ibove,xb,yb,&
                  zb,vxnodes,vynodes,vznodes,vels,varv,varg,vargr,epsi,epsil,&
                  rmsg,rmst,rmsgr,errclt,errclgr,errclg,rtvar,ncomp)

      CALL DATE_AND_TIME(VALUES=time_fin)
      CALL TIMECAL(time_ori,time_fin)

      memory_usage_bytes_Aderi = size(ADERI) * storage_size(ADERI) / 8 ! Size in bytes
      write (*,*) "Approx. memory usage for 'ADERI': ", memory_usage_bytes_Aderi, " bytes"
      write(inout,*) "Approx. memory usage for 'ADERI': ", memory_usage_bytes_Aderi, " bytes"


      memory_usage_bytes_A_s= size(val_ad_s) * storage_size(val_ad_s) / 8 ! Size in bytes
      write (*,*) "Approx. memory usage for 'ADERI sparse' array val: ", memory_usage_bytes_A_s , " bytes"
      write(inout,*) "Approx. memory usage for 'ADERI sparse'array val: ", memory_usage_bytes_A_s, " bytes"



      write(*,*) 'nnz ', nnz ,  'size As', size(val_ad_s) ,'size col_a ', size(columnAD) , &
               'storage_size col_s ' , storage_size(columnAD) 

     !memory_usage_bytes_col_s = size(columnAD) * storage_size(columnAD) 
     ! write (*,*) "Approx. memory usage for 'ADERI sparse' col array: ", memory_usage_bytes_col_s , " bits"
      !write(inout,*) "Approx. memory usage for 'ADERI sparse' col array : ", memory_usage_bytes_col_s, " bits"

      !memory_usage_bytes_row_s = size(rowAD) * storage_size(rowAD) 
      !write (*,*) "Approx. memory usage for 'ADERI sparse' row array: ", memory_usage_bytes_row_s , " bits"
      !write(inout,*) "Approx. memory usage for 'ADERI sparse' row array : ", memory_usage_bytes_row_s, " bits"


      write(*,*)''
      write(*,*)'********************************************************'
      write(*,*)'                      GREETINGS!!'
      write(*,*)'              GOOD END OF JOINT INVERSION'
      write(*,*)'********************************************************'

      write(inout,*)''
      write(inout,*)'********************************************************'
      write(inout,*)'                      GREETINGS!!'
      write(inout,*)'              GOOD END OF JOINT INVERSION'
      write(inout,*)'********************************************************'
!=====================================================================
! Closing the files
!=====================================================================

      close(inout,iostat=ier)
      if(ier.ne.0) then
         write(*,*)''
         write(*,'(''Error in closing the file PARAMETER.OUT, logical &
              &unit '',i4)') inout
         STOP 'in MAIN'
      end if

      if((INV.and.INVD.and.INVV).or.(INV.and.INVGR.and.INVV)) then
         close(correl,iostat=ier)
         if(ier.ne.0) then
            write(*,*)''
            write(*,'(''Error in closing the file CORREL.OUT, logical &
              &unit '',i4)') correl
            STOP 'in MAIN'
         end if

         close(bcoeff,iostat=ier)
         if(ier.ne.0) then
            write(*,*)''
            write(*,'(''Error in closing the file BCOEFF.OUT, logical &
              &unit '',i4)') bcoeff
            STOP 'in MAIN'
         end if
      end if

!=====================================================================
! Deallocation of the allocatable arrays
!=====================================================================
      !Aderi sparse
      DEALLOCATE (val_ad_s)
      DEALLOCATE(columnAD)
      DEALLOCATE(rowAD)

      !status_mkl = mkl_sparse_destroy(Aderi_crs)

      DEALLOCATE (ddtot)
      DEALLOCATE (vdtot)
      DEALLOCATE (grdtot)
      DEALLOCATE (bdtot)
      DEALLOCATE (rovar)
      DEALLOCATE (velvar)
      DEALLOCATE (nlr)
      DEALLOCATE (binit)
      DEALLOCATE (varb)
      DEALLOCATE (ddvr)
      DEALLOCATE (drmax)
      DEALLOCATE (aderi)
      DEALLOCATE (bderi)
      DEALLOCATE (h1)
      DEALLOCATE (caldata)
      DEALLOCATE (diff)
      DEALLOCATE (rmsg)
      DEALLOCATE (rmst)

      if(INVD) then
         DEALLOCATE(x_block)
         DEALLOCATE(y_block)
         DEALLOCATE(ddep)
         DEALLOCATE(dxcor)
         DEALLOCATE(dycor)
      end if

      if(INVV) then
         DEALLOCATE(vxnodes)
         DEALLOCATE(vynodes)
         DEALLOCATE(vznodes)
         DEALLOCATE(v0)
         DEALLOCATE(vels)
         DEALLOCATE(vvari)
         DEALLOCATE(velco)
         DEALLOCATE(G)
      end if
!=====================================================================
! End of the main
!=====================================================================
    END PROGRAM joint_inv