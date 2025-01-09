#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "geometry.h"
#include "grav_prism.h"



tesseroid_(
	x1,x2,y1,y2,z1,z2,density,
	xb,yb,zb,
	resxx,resxy,resxz,resyy,resyz,reszz,
	res2xx,res2xy,res2xz,res2yy,res2yz,res2zz
	)

double *x1,*x2,*y1,*y2,*z1,*z2,*density;
double *xb,*yb,*zb;
double *resxx,*resxy,*resxz,*resyy,*resyz,*reszz;
double *res2xx,*res2xy,*res2xz,*res2yy,*res2yz,*res2zz;

{

		int iii;
	    /**printf("Test de l'implémention de tesseroid dans fortran\n");
	    printf("Les infos sur le prisme sont x1, y1, z1 %.15f %.15f %.15f  \n",*x1,*y1,*z1);
	    printf("Les infos sur le prisme sont x2, y2, z2 %.15f %.15f %.15f  \n",*x2,*y2,*z2);
	    printf("Les infos sur le point sont xb, yb, zb %.15f %.15f %.15f  \n",*xb,*yb,*zb);*/

		PRISM modtest;

	    modtest.x1=*x1;
		modtest.x2=*x2;
		modtest.y1=*y1;
	  	modtest.y2=*y2;
	  	modtest.z1=*z1;
	  	modtest.z2=*z2;
	  	modtest.density=*density;
		modtest.lon=0;
		modtest.lat=0;
		modtest.r=0;

		/**printf("Début du calcul des gradients \n ");*/


/** 1er appel on veut récupérer la vraie réponse*/
		iii=1;

	    *resxx = prism_gxx(modtest,*xb,*yb,-*zb,iii);
	    *resxy = prism_gxy(modtest,*xb,*yb,-*zb,iii);
	    *resxz = prism_gxz(modtest,*xb,*yb,-*zb,iii);
	    *resyy = prism_gyy(modtest,*xb,*yb,-*zb,iii);
	    *resyz = prism_gyz(modtest,*xb,*yb,-*zb,iii);
	    *reszz = prism_gzz(modtest,*xb,*yb,-*zb,iii);


/** 2ème appel on veut récupérer la réponse sans la densité*/
		iii=2;	    
	    *res2xx = prism_gxx(modtest,*xb,*yb,-*zb,iii);
		*res2xy = prism_gxy(modtest,*xb,*yb,-*zb,iii);
		*res2xz = prism_gxz(modtest,*xb,*yb,-*zb,iii);
		*res2yy = prism_gyy(modtest,*xb,*yb,-*zb,iii);
		*res2yz = prism_gyz(modtest,*xb,*yb,-*zb,iii);
	    *res2zz = prism_gzz(modtest,*xb,*yb,-*zb,iii);


        /**printf("On sort de tesseroid \n");*/

        return *resxx,*resxy,*resxz,*resyy,*resyz,*reszz,*res2xx,
        *res2xy,*res2xz,*res2yy,*res2yz,*res2zz;


}	

int pause()
{
	printf("Press enter to continue...\n");
	getchar();
} 
