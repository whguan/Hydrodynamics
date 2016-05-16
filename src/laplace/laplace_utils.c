/*
  laplace_utils.c: utility routines used in adaptive fmm-laplace package
  Copyright (c) 2012 Bo Zhang, Jingfang Huang, Nikos P. Pitsianis, Xiaobai Sun

  This program is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public Licenses as published by 
  the Free Software Foundation, either version 3 of the Licenses, or any 
  later version. 

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "fmm_utils.h"
#include "fmm_ds.h"
#include "fmm_global.h"
#include "fmm_laplace.h"

void adap_fmm_init(const int accuracy, const int nparts)
{
  // The init function completes three tasks: (1) Based on the accuracy 
  // requirement, it determines the length of the multipole/local expansion; 
  // (2) It computes the coefficient invariants in the multipole-to-multipole,
  // local-to-local, and multipole-exponential-local operators. 

  if ( accuracy == 3 ) {
    PTERMS = 9;
    NLAMBS = 9;
  } else if ( accuracy == 6) {
    PTERMS = 18;
    NLAMBS = 18;
  }else if(accuracy == 9){
	PTERMS = 27;
    NLAMBS = 27;
  }else if(accuracy == 12){
	PTERMS = 36;
    NLAMBS = 36;
  }else {
    printf("Error: wrong accuracy input\n");
    exit(-1);
  }

  NPARTS = nparts;
  PGSZ = pow(PTERMS+1,2);
  PTERMS2 = PTERMS*PTERMS;

  // Allocate memory to hold the coefficient invariants in all sorts of operators.
  NUMPHYS = (int *)calloc(NLAMBS, sizeof(int));
  NUMFOUR = (int *)calloc(NLAMBS, sizeof(int));
  WHTS = (double *)calloc(NLAMBS, sizeof(double));
  RLAMS = (double *)calloc(NLAMBS, sizeof(double));
  RDPLUS = (double *)calloc(PGSZ*(2*PTERMS+1), sizeof(double));
  RDMINUS = (double *)calloc(PGSZ*(2*PTERMS+1), sizeof(double));
  RDSQ3 = (double *)calloc(PGSZ*(2*PTERMS+1), sizeof(double));
  RDMSQ3 = (double *)calloc(PGSZ*(2*PTERMS+1), sizeof(double));
  DC = (double *)calloc(pow(2*PTERMS+1,2), sizeof(double));
  YTOPC = (double *)calloc(3721, sizeof(double));
  YTOPCS = (double *)calloc(3721, sizeof(double));
  YTOPCSINV = (double *)calloc(3721, sizeof(double));
  RLSC = (double *)calloc(PGSZ*NLAMBS, sizeof(double));
  CARRAY = (double *)calloc(pow(4*PTERMS+1,2), sizeof(double));

  int allocFailure = (NUMPHYS==0) || (NUMFOUR==0) || (WHTS==0) || 
    (RLAMS==0) || (RDPLUS==0) || (RDMINUS==0) || (RDSQ3==0) || 
    (RDMSQ3==0) || (DC==0) || (YTOPC==0) || (YTOPCS==0) || (YTOPCSINV==0) ||
    (RLSC==0) || (CARRAY==0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Generate coefficient invariants
  frmini(YTOPC, YTOPCS, YTOPCSINV);
  rotgen(PTERMS, CARRAY, RDPLUS, RDMINUS, RDSQ3, RDMSQ3, DC);
  vwts(NLAMBS, RLAMS, WHTS);
  numthetahalf(NLAMBS, NUMFOUR);
  numthetafour(NLAMBS, NUMPHYS);
  rlscini(NLAMBS, PTERMS, RLAMS, RLSC);
  
  NEXPTOT = 0; 
  NTHMAX = 0; 
  NEXPTOTP = 0; 
  int i;
  for ( i = 1; i <= NLAMBS; i++ ) {
    NEXPTOT += NUMFOUR[i-1];
    if ( NUMFOUR[i-1] > NTHMAX ) 
      NTHMAX = NUMFOUR[i-1];
    NEXPTOTP += NUMPHYS[i-1];
  }
  NEXPTOTP /= 2.0;
  NEXPMAX = MAX(NEXPTOT, NEXPTOTP) + 1;

  XS = (dcomplex *)calloc(NEXPMAX*3, sizeof(dcomplex));
  YS = (dcomplex *)calloc(NEXPMAX*3, sizeof(dcomplex));
  ZS = (double *)calloc(NEXPMAX*3, sizeof(double));
  // wenhua Oct. 24
  FEXPE = (dcomplex *)calloc(30000, sizeof(dcomplex));
  FEXPO = (dcomplex *)calloc(30000, sizeof(dcomplex));
  FEXPBACK = (dcomplex *)calloc(30000, sizeof(dcomplex));
  
  allocFailure = (XS==0) || (YS==0 ) || (ZS==0) ||
    (FEXPE==0) || (FEXPO==0) || (FEXPBACK==0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }    
  mkfexp(NLAMBS, NUMFOUR, NUMPHYS, FEXPE, FEXPO, FEXPBACK);
  mkexps(RLAMS, NLAMBS, NUMPHYS, NEXPTOTP, XS, YS, ZS);

  FMMLOC    = (double *)calloc(3*NPARTS, sizeof(double));
  FMMCHARGE = (double *)calloc(NPARTS, sizeof(double));
  //march 3th
  RPYCHARGE = (double *)calloc(3*NPARTS, sizeof(double));
  RPYPOT = (double *)calloc(3*NPARTS, sizeof(double));

  FMMPOT   = (double *)calloc(NPARTS, sizeof(double));
  FMMFIELD = (double *)calloc(3*NPARTS, sizeof(double));
  
  FMMPOTN   = (double *)calloc(NPARTS, sizeof(double));
  FMMFIELDN = (double *)calloc(3*NPARTS, sizeof(double));
  
  FMMDXX = (double *)calloc(NPARTS, sizeof(double));
  FMMDYY = (double *)calloc(NPARTS, sizeof(double));
  FMMDZZ = (double *)calloc(NPARTS, sizeof(double));
  FMMDXY = (double *)calloc(NPARTS, sizeof(double));
  FMMDXZ = (double *)calloc(NPARTS, sizeof(double));
  FMMDYZ = (double *)calloc(NPARTS, sizeof(double));
  
  FMMDXXN = (double *)calloc(NPARTS, sizeof(double));
  FMMDYYN = (double *)calloc(NPARTS, sizeof(double));
  FMMDZZN = (double *)calloc(NPARTS, sizeof(double));
  FMMDXYN = (double *)calloc(NPARTS, sizeof(double));
  FMMDXZN = (double *)calloc(NPARTS, sizeof(double));
  FMMDYZN = (double *)calloc(NPARTS, sizeof(double));

  allocFailure = (FMMLOC==0) || (FMMFIELD==0) || (FMMPOT==0) ||
    (FMMFIELD==0)|| (FMMPOTN==0) || (FMMFIELDN==0) || (FMMDXX==0) || (FMMDYY==0)||
    (FMMDZZ==0) || (FMMDXY==0) || (FMMDXZ==0) || (FMMDYZ==0) || (FMMDXXN==0) || 
    (FMMDYYN==0)|| (FMMDZZN==0) || (FMMDXYN==0) || (FMMDXZN==0) || (FMMDYZN==0);
  if ( allocFailure ) {   
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }


  IFL_UP[0] = 3; IFL_UP[1] = 4; IFL_UP[2] = 2; IFL_UP[3] = 1;
  IFL_UP[4] = 3; IFL_UP[5] = 4; IFL_UP[6] = 2; IFL_UP[7] = 1;
  IFL_DN[0] = 1; IFL_DN[1] = 2; IFL_DN[2] = 4; IFL_DN[3] = 3;
  IFL_DN[4] = 1; IFL_DN[5] = 2; IFL_DN[6] = 4; IFL_DN[7] = 3;
}
// March 3 wenhua
//void adap_fmm_graph (const int nparts, const int s, 
//		     const double *ploc, const double *pcharge)
void adap_fmm_graph(const int nparts, const int s, const double *ploc)
{
  // The graph function partition and relocate data into boxes and build graph
  // and all sorts of associated lists. 
  PERM = (int *)calloc(NPARTS+1, sizeof(int));
  
  int allocFailure = (PERM==0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  if (create_all_box(ploc, nparts, s, PERM, &NBOXES, &NLEV, &SIZE, &BOXES, &CONTENT)) {
    printf("Error in %s, line %d: unable to generate adaptive tree\n", __FILE__, __LINE__);
    exit(-1);
  }
  
  LIST = (fmmlist *)calloc(1+NBOXES, sizeof(fmmlist));
  allocFailure = (LIST==0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  if (create_colleague(BOXES, CONTENT, NLEV, SIZE, LIST)) {
    printf("Error in %s, line %d: unable to generate colleague list\n", __FILE__, __LINE__);
    exit(-1);
  }

  if (create_list134(BOXES, CONTENT, NLEV, NBOXES, SIZE, LIST)) {
    printf("Error in %s, line %d: unable to generate lists 1, 3, and 4\n", __FILE__, __LINE__);
    exit(-1);
  }
  // Relocate input data into boxes
  int i;
  //2015
  //cilk_for ( i = 0; i < NPARTS; i++ ) {
  for ( i = 0; i < NPARTS; i++ ) {
    int j = PERM[i+1]-1;
    FMMLOC[3*i] = ploc[3*j];
    FMMLOC[3*i+1] = ploc[3*j+1];
    FMMLOC[3*i+2] = ploc[3*j+2];
    //FMMCHARGE[i] = pcharge[j];
  }
  // Compute scaling factor
  SCALE = (double *)calloc(1+NLEV, sizeof(double));
  if ( SCALE==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }
  SCALE[0] = 1/SIZE;
  for (i = 1; i <= NLEV; i++ ) 
    SCALE[i] = 2*SCALE[i-1];

}
// March 3
// To assign the charge
// if ifdirect == 0, assign pcharge to FMMCHARGE --- far field caculation using laplace fmm
// if ifdirect != 0, assign rpycharge to RPYCHARGE -- near field caultion using rpy kernel 
void adap_fmm_charge(int ifdirect, const double *charge, const double *rpycharge){
	int allocFailure;
	// Allocate memory to hold multipole, local, and six directional exponential expansions
	MPOLE = (dcomplex *)calloc((1+NBOXES)*PGSZ, sizeof(dcomplex));
	LOCAL = (dcomplex *)calloc((1+NBOXES)*PGSZ, sizeof(dcomplex));
	LEXPU = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	LEXPD = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	LEXPN = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	LEXPS = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	LEXPE = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	LEXPW = (dcomplex *)calloc((1+NBOXES)*NEXPMAX, sizeof(dcomplex));
	allocFailure = (MPOLE==0) || (LOCAL==0) || (LEXPU==0) || (LEXPD==0) ||
    (LEXPN==0) || (LEXPS==0) || (LEXPE==0) || (LEXPW==0);
	if ( allocFailure ) {
		printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
		exit(-1);
	}
	int i;
	for ( i = 0; i < NPARTS; i++ ) {
	        int j = PERM[i+1]-1;
	// ifdirect is deciding whether pcharge or rpycharge has been assigned.
		if(ifdirect==0){
			FMMCHARGE[i] = charge[j];
			FMMPOT[i] = 0;
			FMMFIELD[3*i] =0;
			FMMFIELD[3*i+1] = 0;
			FMMFIELD[3*i+2] = 0;
			FMMDXX[i] = 0;
			FMMDYY[i] = 0;
			FMMDZZ[i] = 0;
			FMMDXY[i] = 0;
			FMMDXZ[i] = 0;
			FMMDYZ[i] = 0;
		
			FMMPOTN[i] = 0;
			FMMFIELDN[3*i] =0;
			FMMFIELDN[3*i+1] = 0;
			FMMFIELDN[3*i+2] = 0;
			FMMDXXN[i] = 0;
			FMMDYYN[i] = 0;
			FMMDZZN[i] = 0;
			FMMDXYN[i] = 0;
			FMMDXZN[i] = 0;
			FMMDYZN[i] = 0;
		}else{
		       // Near field calculation
            		RPYCHARGE[3*i] = rpycharge[3*j];
           		RPYCHARGE[3*i+1] = rpycharge[3*j+1];
          		RPYCHARGE[3*i+2] = rpycharge[3*j+2];
                        //note : here to avoid calloc new memory to for pot of RPY in near field.
			// FMMFIELDN is used to store the pot of RPY
		   	 FMMFIELDN[3*i] = 0;
			 FMMFIELDN[3*i+1] =0;
			 FMMFIELDN[3*i+2] =0;

		}

	}

}

void adap_fmm_post(int ifdirect, double *RPY, int ifpot, double *pot, int ifgrad, double *field, 
	int ifhess, double *dxx, double *dyy, double *dzz,double *dxy, double *dxz, double *dyz)
{
  // The post function completes two tasks: (1) It writes out data in original 
  // input order; (2) It frees all memory space allocated in the graph function.

  int i;
    if(ifdirect==0){
        //Laplace fmm  far field March 14th		
		for( i=0;i<NPARTS; i++){
			int j = PERM[i+1]-1;
			if(ifpot!=0)
			pot[j] = FMMPOT[i];//+FMMPOTN[i];
				
			if(ifgrad!=0){
				field[3*j]   = FMMFIELD[3*i];//+FMMFIELDN[3*i];
				field[3*j+1] = FMMFIELD[3*i+1];//+FMMFIELDN[3*i+1];
				field[3*j+2] = FMMFIELD[3*i+2];//+FMMFIELDN[3*i+2]; 
			}
		
			if(ifhess!=0){
				dxx[j] = FMMDXX[i];// + FMMDXXN[i];
				dyy[j] = FMMDYY[i];// + FMMDYYN[i];
				dzz[j] = FMMDZZ[i];// + FMMDZZN[i];
				dxy[j] = FMMDXY[i];// + FMMDXYN[i];
				dxz[j] = FMMDXZ[i];// + FMMDXZN[i];
				dyz[j] = FMMDYZ[i];// + FMMDYZN[i];
			}
	
		}
    } else {
	    /*printf("after RPYPOT\n");
	    for(i=0;i<NPARTS; i++){
		printf("%f %f %f \n", RPYPOT[3*i], RPYPOT[3*i+1], RPYPOT[3*i+2]);
	    } */
	    // Just the near field for RPY fmm
		for(i=0;i<NPARTS; i++){
			int j = PERM[i+1]-1;
			if(ifpot!=0){
				RPY[3*j]  += RPYPOT[3*i];
				RPY[3*j+1] += RPYPOT[3*i+1];
				RPY[3*j+2] += RPYPOT[3*i+2];
	
			}  
		}
    }
	free(MPOLE);
	free(LOCAL);
	free(LEXPU);
	free(LEXPD);
	free(LEXPN);
	free(LEXPS);
	free(LEXPE);
	free(LEXPW);
}

void adap_fmm_clean(void)
{
  // Free memory space allocated in the graph function
  free(PERM);
  free(BOXES);
  free(CONTENT);
  free_list(LIST, NBOXES);
  /*free(MPOLE);
  free(LOCAL);
  free(LEXPU);
  free(LEXPD);
  free(LEXPN);
  free(LEXPS);
  free(LEXPE);
  free(LEXPW);*/
  free(SCALE);
  // The clean function frees all the memory space allocated in the init function.
  free(NUMPHYS);
  free(NUMFOUR);
  free(WHTS);
  free(RLAMS);
  free(RDPLUS);
  free(RDMINUS);
  free(RDSQ3);
  free(RDMSQ3);
  free(DC);
  free(YTOPC);
  free(YTOPCS);
  free(YTOPCSINV);
  free(RLSC);
  free(CARRAY);
  free(XS);
  free(YS);
  free(ZS); 
  free(FEXPE); 
  free(FEXPO);
  free(FEXPBACK); 
  free(FMMLOC);
  free(FMMCHARGE);
  free(RPYCHARGE);

  free(RPYPOT);
  free(FMMPOT); 
  free(FMMFIELD); 
  free(FMMPOTN);
  free(FMMFIELDN);
  
  free(FMMDXX);
  free(FMMDYY);
  free(FMMDZZ);
  free(FMMDXY);
  free(FMMDXZ);
  free(FMMDYZ);
  
  free(FMMDXXN);
  free(FMMDYYN);
  free(FMMDZZN);
  free(FMMDXYN);
  free(FMMDXZN);
  free(FMMDYZN);
}

void rlscini(const int nlambs, const int pterms, const double *rlams, double *rlsc)
{
  double *factorial = (double *)calloc(2*pterms+1, sizeof(double));
  double *rlampow = (double *)calloc(1+pterms, sizeof(double));
  if ( factorial==0 || rlampow==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit (-1);
  }

  factorial[0] = 1;
  int i;
  for ( i = 1; i <= 2*pterms; i++)
    factorial[i] = factorial[i-1]*sqrt(i);

  int nell;
  for ( nell = 0; nell < nlambs; nell++ ) {
    double rmul = rlams[nell];
    rlampow[0] = 1;
    int j;
    for ( j = 1;  j <= pterms; j++)
      rlampow[j] = rlampow[j-1]*rmul;

    for ( j = 0; j <=pterms; j++) {
      int k;
      for ( k = 0; k <= j; k++ ) {
	rlsc[j+k*(pterms+1)+nell*PGSZ] = rlampow[j]/factorial[j-k]/factorial[j+k];
      }
    }
  }

  free (factorial);
  free (rlampow);
}

void mkfexp(const int nlambs, const int *numfour, const int *numphys,
	    dcomplex *fexpe, dcomplex *fexpo, dcomplex *fexpback)
{
  int i, j, nm, nexte, nexto, next, nalpha, nalpha2;
  double halpha, alpha;

  nexte = 0;
  nexto = 0;
  for ( i = 0; i < nlambs; i++ ) {
    nalpha = numphys[i];
    nalpha2 = nalpha/2;
    halpha = 2.0*M_PI/nalpha;
    for ( j = 1; j <= nalpha2; j++ ) {
      alpha = (j-1)*halpha;
      for ( nm = 2; nm <= numfour[i]; nm=nm+2 ) {
	fexpe[nexte] = cexp( (nm-1)*_Complex_I*alpha);
	nexte += 1;
      }
      for ( nm=3; nm <= numfour[i]; nm=nm+2) {
	fexpo[nexto] = cexp( (nm-1)*_Complex_I*alpha);
	nexto += 1;
      }
    }
  }

  next = 0;
  for ( i = 0; i < nlambs; i++) {
    nalpha = numphys[i];
    nalpha2 = nalpha/2;
    halpha = 2.0*M_PI/nalpha;
    for ( nm = 3; nm<=numfour[i]; nm=nm+2) {
      for ( j = 1; j <= nalpha2; j++) {
	alpha = (j-1)*halpha;
	fexpback[next] = cexp(-(nm-1)*_Complex_I*alpha);
	next += 1;
      }
    }
    for ( nm = 2; nm<=numfour[i]; nm=nm+2) {
      for ( j= 1; j <= nalpha2; j++) {
	alpha = (j-1)*halpha;
	fexpback[next] = cexp(-(nm-1)*_Complex_I*alpha);
	next += 1;
      }
    }
  }
}

void mkexps(const double *rlams, const int nlambs, const int *numphys,
	    const int nexptotp, dcomplex *xs, dcomplex *ys, double *zs)
{
  int ntot, nell, mth, ncurrent;
  double hu, u;

  ntot = 0;
  for ( nell = 0; nell < nlambs; nell++) {
    hu = 2.0*M_PI/numphys[nell];
    for ( mth = 0; mth < numphys[nell]/2; mth++ ) {
      u = mth*hu;
      ncurrent=3*(ntot+mth);
      zs[ncurrent] = exp( -rlams[nell] );
      zs[ncurrent+1] = zs[ncurrent]*zs[ncurrent];
      zs[ncurrent+2] = zs[ncurrent]*zs[ncurrent+1];
      xs[ncurrent] = cexp(_Complex_I*cos(u)*rlams[nell]);
      xs[ncurrent+1] = xs[ncurrent]*xs[ncurrent];
      xs[ncurrent+2] = xs[ncurrent]*xs[ncurrent+1];
      ys[ncurrent] = cexp(_Complex_I*sin(u)*rlams[nell]);
      ys[ncurrent+1] = ys[ncurrent]*ys[ncurrent];
      ys[ncurrent+2] = ys[ncurrent+1]*ys[ncurrent];
    }
    ntot += numphys[nell]/2;
  }
}

//Compute Legendre polynomial
void lgndr(const int nmax, const double x, double *y)
{
  // y is used as y(0:nmax,0:nmax) in JF's fortran routine
  int m, n, offset1, offset2, offset3;
  double u;

  n = (nmax+1)*(nmax+1);
  for (m=0; m<n; m++)
    y[m] = 0;

  offset1 = nmax+2;

  u = -sqrt(1-x*x);
  y[0] = 1;

  // m = 0 case
  y[1] = x*y[0];
  for ( n = 2; n<=nmax; n++ ) 
    y[n] = ((2*n-1)*x*y[n-1]-(n-1)*y[n-2])/n;
  
  // m = 1:nmax-1
  for ( m = 1; m<=nmax-1; m++ ) {
    offset2 = m*offset1;
    y[offset2] = y[offset2-offset1]*u*(2*m-1);
    y[offset2+1] = y[offset2]*x*(2*m+1);
    for ( n = m+2; n <=nmax; n++ ) {
      offset3 = n+m*(nmax+1);
      y[offset3] = ( (2.0*n-1)*x*y[offset3-1] - (n+m-1)*y[offset3-2])/(n-m);
    }
  }
  
  // m = nmax case
  y[nmax+nmax*(nmax+1)]  = y[nmax-1+(nmax-1)*(nmax+1)]*u*(2*nmax-1);

 return;
}

void fstrtn(const int pterms, double *d, const double *sqc, const double theta)
{
  // d (0:pterms,0:pterms, -pterms:pterms);
  // sqc(0:2*pterms, 0:2*pterms);
  int ij, im, imp;
  double ctheta, stheta, hsthta, cthtap, cthtan, precis, ww;
  precis = 1.0e-19;
  ww = sqrt(2)/2;

  ctheta = cos(theta);
  if ( fabs(ctheta) <= precis) 
    ctheta = 0.0;

  stheta = sin(-theta);
  if ( fabs(stheta) <= precis)
    stheta = 0.0;

  hsthta = ww*stheta;
  cthtap = ww*(1.0+ctheta);
  cthtan = -ww*(1.0-ctheta);

  d[pterms*PGSZ] = 1.0;

  for ( ij = 1; ij <= pterms; ij++ ) {
    // compute result for m'=0 case, use formula (1)
    for ( im = -ij; im <= -1; im++ ) {
      d[ij+(im+pterms)*PGSZ] = -sqc[ij-im+2*(1+2*pterms)]*d[ij-1 + (im+1+pterms)*PGSZ];

      if ( im > 1-ij ) {
	d[ij+(im+pterms)*PGSZ] += sqc[ij+im+2*(1+2*pterms)]*d[ij-1+(im-1+pterms)*PGSZ];
      }

      d[ij+(im+pterms)*PGSZ] *= hsthta;

      if ( im > -ij ) {
	d[ij+(im+pterms)*PGSZ] += d[ij-1+(im+pterms)*PGSZ]*ctheta*
	  sqc[ij+im+2*pterms+1]*sqc[ij-im + 2*pterms+1];
      }

      d[ij+(im+pterms)*PGSZ] /= ij;
    }

    d[ij+pterms*PGSZ] = d[ij-1+pterms*PGSZ]*ctheta;
    if ( ij > 1 ) {
      d[ij+pterms*PGSZ] = d[ij+pterms*PGSZ] + hsthta*sqc[ij+2*(1+2*pterms)]*(
						     d[ij-1+(-1+pterms)*PGSZ]+d[ij-1+(1+pterms)*PGSZ])/ij;
    }

    for ( im = 1; im <= ij; im++ ) {
      d[ij + (im+pterms)*PGSZ] = -sqc[ij+im+2*(1+2*pterms)]*d[ij-1+(im-1+pterms)*PGSZ];

      if ( im < ij-1 ) {
	d[ij+(im+pterms)*PGSZ] += sqc[ij-im+2*(1+2*pterms)]*d[ij-1+(im+1+pterms)*PGSZ];
      }

      d[ij+(im+pterms)*PGSZ] *= hsthta;

      if ( im < ij ) {
	d[ij+(im+pterms)*PGSZ] += d[ij-1+(im+pterms)*PGSZ]*ctheta*
	  sqc[ij+im+2*pterms+1]*sqc[ij-im+2*pterms+1];
      }

      d[ij+(im+pterms)*PGSZ] /= ij;
    }

    // compute the result for 0<m'<=j case, use formula (2)
    for ( imp = 1; imp <= ij; imp++ ) {
      for ( im = -ij; im <= -1; im++ ) {
	d[ij + imp*(pterms+1)+(im+pterms)*PGSZ] = 
	  d[ij-1 + (imp-1)*(pterms+1)+(im+1+pterms)*PGSZ]*cthtan*sqc[ij-im+2*(2*pterms+1)];

	if ( im > 1 - ij ) {
	  d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] -= 
	    d[ij-1+(imp-1)*(pterms+1)+(im-1+pterms)*PGSZ]*cthtap*sqc[ij+im+2*(pterms*2+1)];
	}

	if ( im > -ij ) {
	  d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] +=
	    d[ij-1+(imp-1)*(pterms+1)+(im+pterms)*PGSZ]*stheta*
	    sqc[ij+im+2*pterms+1]*sqc[ij-im+2*pterms+1];
	}

	d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] *= ww/sqc[ij+imp+2*(1+2*pterms)];
      }

      d[ij+imp*(pterms+1)+pterms*PGSZ] = ij*stheta*d[ij-1+(imp-1)*(pterms+1)+pterms*PGSZ];
      if ( ij > 1 ) {
	d[ij+imp*(pterms+1)+pterms*PGSZ] -= 
	  sqc[ij+2*(1+2*pterms)]*( d[ij-1+(imp-1)*(pterms+1)+(-1+pterms)*PGSZ]*cthtap
				   + d[ij-1+(imp-1)*(pterms+1)+(1+pterms)*PGSZ]*cthtan );
      }

      d[ij+imp*(pterms+1)+pterms*PGSZ] *= ww/sqc[ij+imp+2*(pterms*2+1)];

      for ( im = 1; im <= ij; im++ ) {
	d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] = 
	  d[ij-1+(imp-1)*(pterms+1)+(im-1+pterms)*PGSZ]*cthtap*sqc[ij+im+2*(1+2*pterms)];

	if ( im < ij - 1 ) {
	  d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] -= 
	    d[ij-1+(imp-1)*(pterms+1)+(im+1+pterms)*PGSZ]*cthtan*sqc[ij-im+2*(1+2*pterms)];
	}

	if ( im < ij ) {
	  d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] += 
	    d[ij-1+(imp-1)*(pterms+1)+(im+pterms)*PGSZ]*stheta*sqc[ij+im+1+2*pterms]
	    *sqc[ij-im+1+2*pterms];
	}

	d[ij+imp*(pterms+1)+(im+pterms)*PGSZ] *= ww/sqc[ij+imp+2*(1+2*pterms)];
      }
    }
  }
}

void bnlcft(double *c, const int pterms)
{
  int n, m;
  for ( n = 0; n <= pterms; n++ ) {
    c[n] = 1.0;
  }

  for ( m = 1; m <= pterms; m++ ) {
    int offset = m*(pterms+1);
    int offset1 = (m-1)*(pterms+1);
    c[m+offset] = 1.0;
    for ( n = m+1; n <= pterms; n++ ) {
      c[n+offset] = c[n-1+offset]+c[n-1+offset1];
    }
  }

  for ( m = 1; m <= pterms; m++ ) {
    int offset = m*(pterms+1);
    for ( n = m+1; n<=pterms; n++) {
      c[n+offset] = sqrt( c[n+offset] );
    }
  }
}

void rotgen(const int pterms, const double *carray, double *rdpi2, double *rdmpi2,
	    double *rdsq3, double *rdmsq3, double *dc)
{
  bnlcft(dc, 2*pterms);
  double theta = M_PI/2;
  fstrtn(pterms, rdpi2, dc, theta);
  theta = -theta;
  fstrtn(pterms, rdmpi2, dc, theta);
  theta = acos(sqrt(3)/3);
  fstrtn(pterms, rdsq3, dc, theta);
  theta = acos(-sqrt(3)/3);
  fstrtn(pterms, rdmsq3, dc, theta);
}
// Compute the coeffients
void frmini(double *c, double *cs, double *csinv)
{
  int ell, m;
  double d  = 1.0;

  double *factorial = (double *)calloc(121, sizeof(double));
  if ( factorial==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  factorial[0] = d;
  for ( ell = 1; ell <= 120; ell++ ) {
    d *= sqrt(ell);
    factorial[ell] = d;
  }

  cs[0] = 1.0;
  csinv[0] = 1.0;
  for ( ell = 1; ell <= 60; ell++ ) {
    for ( m = 0; m <= ell; m++ ) {
      c[ell+m*61] = factorial[ell-m]/factorial[ell+m];
      csinv[ell+m*61] = factorial[ell-m]*factorial[ell+m];
      cs[ell+m*61] = 1.0/csinv[ell+m*61];
    }
  }

  free(factorial);
}

void vecAdd(const dcomplex *b, dcomplex *a, const int nelements)
{
  int i;
  for ( i = 0; i < nelements; i++)
    a[i] += b[i];
}


void vwts(const int nlambs, double *rlams, double *whts)
/*
  vwts returns a set of Gaussian nodes and weights for integrating the functions
j0(rx)*exp(zx) with respect to x over the range x=0 to x=infty. The range of z is 
assumed to be [1,4], and the range of r is [0, 4*sqrt(2)]. 

nlambs is the number of weights and nodes in the quadature,  rlams is the nodes, 
and whts is the weights, nlambs should be an integer within [2, 39].

The approximate accuracy of the quadrature is as follows: 
	n	maximum error
          2     0.15318e+00
          3     0.76505e-01
          4     0.32149e-01
          5     0.15630e-01
          6     0.75110e-02
          7     0.35030e-02
          8     0.16243e-02
          9     0.72230e-03
        10     0.33074e-03
        11     0.15035e-03
        12     0.70952e-04
        13     0.31751e-04
        14     0.14589e-04
        15     0.64300e-05
        16     0.29477e-05
        17     0.13222e-05
        18     0.61488e-06
        19     0.27435e-06
        20     0.12534e-06
        21     0.55324e-07
        22     0.25257e-07
        23     0.11293e-07
        24     0.52063e-08
        25     0.23256e-08
        26     0.10580e-08
        27     0.46835e-09
        28     0.21286e-09
        29     0.95164e-10
        30     0.43599e-10
        31     0.19516e-10
        32     0.88491e-11
        33     0.39313e-11
        34     0.17821e-11
        35     0.79603e-12
        36     0.36460e-12
        37     0.16331e-12
        38     0.73497e-13
        39     0.31530e-13
 */
{
  if ( nlambs < 2 || nlambs > 39 ) {
    printf("Error in %s, line %d: invalid input in vwts()\n", __FILE__, __LINE__);
    exit (-1);
  }

  if ( nlambs == 9 ) {
    // 3-digit case
    rlams[0]=0.99273996739714473469540223504736787e-01;
    rlams[1]=0.47725674637049431137114652301534079e+00;
    rlams[2]=0.10553366138218296388373573790886439e+01;
    rlams[3]=0.17675934335400844688024335482623428e+01;
    rlams[4]=0.25734262935147067530294862081063911e+01;
    rlams[5]=0.34482433920158257478760788217186928e+01;
    rlams[6]=0.43768098355472631055818055756390095e+01;
    rlams[7]=0.53489575720546005399569367000367492e+01;
    rlams[8]=0.63576578531337464283978988532908261e+01;
    whts[0]=0.24776441819008371281185532097879332e+00;
    whts[1]=0.49188566500464336872511239562300034e+00;
    whts[2]=0.65378749137677805158830324216978624e+00;
    whts[3]=0.76433038408784093054038066838984378e+00;
    whts[4]=0.84376180565628111640563702167128213e+00;
    whts[5]=0.90445883985098263213586733400006779e+00;
    whts[6]=0.95378613136833456653818075210438110e+00;
    whts[7]=0.99670261613218547047665651916759089e+00;
    whts[8]=0.10429422730252668749528766056755558e+01;
  } else if ( nlambs == 18 ) {
    // 6-digit case
    rlams[0]=0.52788527661177607475107009804560221e-01;
    rlams[1]=0.26949859838931256028615734976483509e+00;
    rlams[2]=0.63220353174689392083962502510985360e+00;
    rlams[3]=0.11130756427760852833586113774799742e+01;
    rlams[4]=0.16893949614021379623807206371566281e+01;
    rlams[5]=0.23437620046953044905535534780938178e+01;
    rlams[6]=0.30626998290780611533534738555317745e+01;
    rlams[7]=0.38356294126529686394633245072327554e+01;
    rlams[8]=0.46542473432156272750148673367220908e+01;
    rlams[9]=0.55120938659358147404532246582675725e+01;
    rlams[10]=0.64042126837727888499784967279992998e+01;
    rlams[11]=0.73268800190617540124549122992902994e+01;
    rlams[12]=0.82774009925823861522076185792684555e+01;
    rlams[13]=0.92539718060248947750778825138695538e+01;
    rlams[14]=0.10255602723746401139237605093512684e+02;
    rlams[15]=0.11282088297877740146191172243561596e+02;
    rlams[16]=0.12334067909676926788620221486780792e+02;
    rlams[17]=0.13414920240172401477707353478763252e+02;
    whts[0]=0.13438265914335215112096477696468355e+00;
    whts[1]=0.29457752727395436487256574764614925e+00;
    whts[2]=0.42607819361148618897416895379137713e+00;
    whts[3]=0.53189220776549905878027857397682965e+00;
    whts[4]=0.61787306245538586857435348065337166e+00;
    whts[5]=0.68863156078905074508611505734734237e+00;
    whts[6]=0.74749099381426187260757387775811367e+00;
    whts[7]=0.79699192718599998208617307682288811e+00;
    whts[8]=0.83917454386997591964103548889397644e+00;
    whts[9]=0.87570092283745315508980411323136650e+00;
    whts[10]=0.90792943590067498593754180546966381e+00;
    whts[11]=0.93698393742461816291466902839601971e+00;
    whts[12]=0.96382546688788062194674921556725167e+00;
    whts[13]=0.98932985769673820186653756536543369e+00;
    whts[14]=0.10143828459791703888726033255807124e+01;
    whts[15]=0.10400365437416452252250564924906939e+01;
    whts[16]=0.10681548926956736522697610780596733e+01;
    whts[17]=0.11090758097553685690428437737864442e+01;
  }else if(nlambs ==27){
      rlams[0]=0.35834993247954083361861421508365311e-01;
	  rlams[1]=0.18605352495124757861155728733137948e+00;
	  rlams[2]=0.44642194330546725034025712375296280e+00;
	  rlams[3]=0.80414605853435716653621057048439980e+00;
	  rlams[4]=0.12462239924414626468518463298096322e+01;
	  rlams[5]=0.17611257461239471222569363817456178e+01;
	  rlams[6]=0.23391453785205094106913747964426875e+01;
	  rlams[7]=0.29722080942457780317056403873721138e+01;
	  rlams[8]=0.36535650947301854252202701900387183e+01;
	  rlams[9]=0.43775326951870896508012265258003026e+01;
	  rlams[10]=0.51393010711012010460763121955096722e+01;
	  rlams[11]=0.59347963207816585295972799940500408e+01;
	  rlams[12]=0.67605758319435445002909546019509435e+01;
	  rlams[13]=0.76137438021356000916739503736607730e+01;
	  rlams[14]=0.84918804165657331139982488821260631e+01;
	  rlams[15]=0.93929821745587691594892021385021508e+01;
	  rlams[16]=0.10315412711899460518338855763431638e+02;
	  rlams[17]=0.11257864229560141211550217121839523e+02;
	  rlams[18]=0.12219330080395609527954547957051545e+02;
	  rlams[19]=0.13199089643751765521528795943595469e+02;
	  rlams[20]=0.14196707628360135444722800457384437e+02;
	  rlams[21]=0.15212051628819105886236684455070645e+02;
	  rlams[22]=0.16245334166526408381514556822367013e+02;
	  rlams[23]=0.17297187705846923222452460322529078e+02;
	  rlams[24]=0.18368783996172027173088281415402889e+02;
	  rlams[25]=0.19462114300213894324542707181535661e+02;
	  rlams[26]=0.20582673852051787122263704077340662e+02;
	  whts[0]=0.91634452557743162337544617912499234e-01;
	  whts[1]=0.20725798653776011937210910218709614e+00;
	  whts[2]=0.31128301536359448409285732850548811e+00;
	  whts[3]=0.40197254484414085551335915624804329e+00;
	  whts[4]=0.48025776040426659818649568478576839e+00;
	  whts[5]=0.54793227641803410321585943165700883e+00;
	  whts[6]=0.60676658418922713167376059573143721e+00;
	  whts[7]=0.65823842224192508076896501734154299e+00;
	  whts[8]=0.70353143374483517735740178977721371e+00;
	  whts[9]=0.74360437862720618618084245099453256e+00;
	  whts[10]=0.77925510412986809249957786960294470e+00;
	  whts[11]=0.81116302338313539532776985652162693e+00;
	  whts[12]=0.83991454783642205406124503497267142e+00;
	  whts[13]=0.86601910131774351153666202662861906e+00;
	  whts[14]=0.88992107452719526605733335600234568e+00;
	  whts[15]=0.91201046747319902685546821885509416e+00;
	  whts[16]=0.93263337746385499560375365035724826e+00;
	  whts[17]=0.95210285562789409574691035231808200e+00;
	  whts[18]=0.97071063776060639849418976154993288e+00;
	  whts[19]=0.98874065527360133209811010601697490e+00;
	  whts[20]=0.10064860790234053666836189222522080e+01;
	  whts[21]=0.10242732841383590525197178067173809e+01;
	  whts[22]=0.10424992963077313046227345694205724e+01;
	  whts[23]=0.10616963916554709168593717549811117e+01;
	  whts[24]=0.10826735650889751649827985602314584e+01;
	  whts[25]=0.11072236156410122376314575376454741e+01;
	  whts[26]=0.11468295974452529240039666547090746e+01;
  
  } else if(nlambs == 36){
  	  rlams[0]=0.27142056328972239548358302840824763e-01;
	  rlams[1]=0.14178681907151427510349606109230081e+00;
	  rlams[2]=0.34342381977838998263763414797722362e+00;
	  rlams[3]=0.62543278953658920649161245819414034e+00;
	  rlams[4]=0.98021315858273505572384465267532505e+00;
	  rlams[5]=0.14002363417962260250249073578743264e+01;
	  rlams[6]=0.18785945055991226659841686341678724e+01;
	  rlams[7]=0.24091879673843461517890318646095693e+01;
	  rlams[8]=0.29867138306149492166241543600335717e+01;
	  rlams[9]=0.36065737319948958372606284683570266e+01;
	  rlams[10]=0.42647654793772025172415851557161659e+01;
	  rlams[11]=0.49577865215147900457282048591878265e+01;
	  rlams[12]=0.56825569402033835331167210824787617e+01;
	  rlams[13]=0.64363605312603828778605929983314127e+01;
	  rlams[14]=0.72167996452014815389475188567303121e+01;
	  rlams[15]=0.80217595113515507421197980875149369e+01;
	  rlams[16]=0.88493788614634816269699513213708997e+01;
	  rlams[17]=0.96980248504438417711526199127547443e+01;
	  rlams[18]=0.10566271182166131481494630861561745e+02;
	  rlams[19]=0.11452878934627758411579634412191808e+02;
	  rlams[20]=0.12356779911330106003219952981453389e+02;
	  rlams[21]=0.13277062520347925556052359752357006e+02;
	  rlams[22]=0.14212960280227379783468677487690002e+02;
	  rlams[23]=0.15163843131286117937861490645445883e+02;
	  rlams[24]=0.16129211834035846351298459921963513e+02;
	  rlams[25]=0.17108695899939579732063066330738366e+02;
	  rlams[26]=0.18102055768589796258538626716472208e+02;
	  rlams[27]=0.19109190390395511371934844646602869e+02;
	  rlams[28]=0.20130152154913822926118882605805993e+02;
	  rlams[29]=0.21165172450830002759403214440681040e+02;
	  rlams[30]=0.22214703800018561707929620752111077e+02;
	  rlams[31]=0.23279488597200625576988386455923319e+02;
	  rlams[32]=0.24360673016746289931688806973397732e+02;
	  rlams[33]=0.25460002587127451789683618699200451e+02;
	  rlams[34]=0.26580225791803091084375409991480410e+02;
	  rlams[35]=0.27728035052626292866762014455161989e+02;
	  whts[0]=0.69511944436498607213792411130270921e-01;
	  whts[1]=0.15909081406147243531457036169740604e+00;
	  whts[2]=0.24305240853533374711936687617708230e+00;
	  whts[3]=0.31968136778214911730700009684369434e+00;
	  whts[4]=0.38861583821018325091145584337937180e+00;
	  whts[5]=0.45027651490975162396068753878353164e+00;
	  whts[6]=0.50542366833399765546630533208372071e+00;
	  whts[7]=0.55488112295767544335234333630069159e+00;
	  whts[8]=0.59940607536473111682795433807768859e+00;
	  whts[9]=0.63964831820243539528547671579872258e+00;
	  whts[10]=0.67615269572807168430017554783262312e+00;
	  whts[11]=0.70937703945268570926430129475193098e+00;
	  whts[12]=0.73971183935845019608024131230195053e+00;
	  whts[13]=0.76749632153441160742346482948050834e+00;
	  whts[14]=0.79302992954234163835991466839914210e+00;
	  whts[15]=0.81658000692285714894325110435602255e+00;
	  whts[16]=0.83838691595559167168971725914161652e+00;
	  whts[17]=0.85866766530031546356127591934637167e+00;
	  whts[18]=0.87761878886172950409161330753704533e+00;
	  whts[19]=0.89541891941597706594535566182457842e+00;
	  whts[20]=0.91223129526075630302273111738031730e+00;
	  whts[21]=0.92820632453781548587556926577235572e+00;
	  whts[22]=0.94348428787583971111274649956612848e+00;
	  whts[23]=0.95819826511204764241114162359735928e+00;
	  whts[24]=0.97247741465258175086461278624483384e+00;
	  whts[25]=0.98645081746307394787720568274380639e+00;
	  whts[26]=0.10002522385876044808128426666371524e+01;
	  whts[27]=0.10140264169597386079146872361889109e+01;
	  whts[28]=0.10279380258279977589097597956424579e+01;
	  whts[29]=0.10421853786456081181199806451331824e+01;
	  whts[30]=0.10570235062817368021370612041209824e+01;
	  whts[31]=0.10728043592728884192410987452603877e+01;
	  whts[32]=0.10900568310850353714869243049179204e+01;
	  whts[33]=0.11096639729316220002175441550207324e+01;
	  whts[34]=0.11336107695019650432044500121264718e+01;
	  whts[35]=0.11734853243038561032562938635237515e+01;
  
  }
}

void numthetahalf(const int nlambs, int *numtets)
/*
  numthetahalf returns the number of Fourier modes needed in the phi
integral for each of the discrete lambda values given by Norman's quadratures. 
 */
{
  if ( nlambs == 9 ) {
    numtets[ 0] = 2;
    numtets[ 1] = 4;
    numtets[ 2] = 4;
    numtets[ 3] = 6;
    numtets[ 4] = 6;
    numtets[ 5] = 4;
    numtets[ 6] = 6;
    numtets[ 7] = 4;
    numtets[ 8] = 2;
  } else if ( nlambs == 18 ) {
    numtets[ 0] = 4;
    numtets[ 1] = 6;
    numtets[ 2] = 6;
    numtets[ 3] = 8;
    numtets[ 4] = 8;
    numtets[ 5] = 8;
    numtets[ 6] =10;
    numtets[ 7] =10;
    numtets[ 8] =10;
    numtets[9] =10;
    numtets[10] =12;
    numtets[11] =12;
    numtets[12] =12;
    numtets[13] =12;
    numtets[14] =12;
    numtets[15] =12;
    numtets[16] = 8;
    numtets[17] = 2;
  }else if(nlambs == 27){
    numtets[ 0] = 4;
	numtets[ 1] = 6;
	numtets[ 2] = 8;
	numtets[ 3] = 8;
	numtets[ 4] =10;
	numtets[ 5] =10;
    numtets[ 6] =12;
	numtets[ 7] =12;
	numtets[ 8] =14;
	numtets[9] =14;
	numtets[10] =14;												            
	numtets[11] =16;
	numtets[12] =16;
	numtets[13] =16;
	numtets[14] =16;
	numtets[15] =16;
	numtets[16] =18;	
	numtets[17] =18;
	numtets[18] =18;
	numtets[19] =18;
	numtets[20] =18;
	numtets[21] =18;
	numtets[22] =18;
	numtets[23] =16;
	numtets[24] =16;
	numtets[25] =12;
	numtets[26] = 2;
  
  }else if(nlambs == 36){
    numtets[ 0] = 6;
	numtets[ 1] = 8;
    numtets[ 2] = 8;
    numtets[ 3] =10;
    numtets[ 4] =12;
    numtets[ 5] =12;
    numtets[ 6] =14;
    numtets[ 7] =14;
    numtets[ 8] =16;
    numtets[ 9] =16;
    numtets[10] =16;
    numtets[11] =18;
    numtets[12] =18;
    numtets[13] =20;
    numtets[14] =20;
    numtets[15] =20;
    numtets[16] =20;
    numtets[17] =22;
    numtets[18] =22;
    numtets[19] =22;
    numtets[20] =22;
    numtets[21] =22;
    numtets[22] =24;
	numtets[23] =24;
    numtets[24] =24;
    numtets[25] =24;
    numtets[26] =24;
    numtets[27] =24;
    numtets[28] =24;
    numtets[29] =24;
    numtets[30] =24;
    numtets[31] =22;
    numtets[32] =22;
    numtets[33] =18;
    numtets[34] =16;
    numtets[35] = 2;
  }
}

void numthetafour(const int nlambs, int *numtets)
/*
  numthetafour returns the number of fourier modes needed in 
the phi-integral for each of the discrete lambda values given by Norman's
quadrature
*/
{
  if ( nlambs == 9 ) {
    numtets[ 0] =  4;
    numtets[ 1] =  8;
    numtets[ 2] = 12;
    numtets[ 3] = 16;
    numtets[ 4] = 20;
    numtets[ 5] = 20;
    numtets[ 6] = 24;
    numtets[ 7] =  8;
    numtets[ 8] =  2;
  } else if ( nlambs == 18 ) {
    numtets[ 0] =  6;
    numtets[ 1] =  8;
    numtets[ 2] = 12;
    numtets[ 3] = 16;
    numtets[ 4] = 20;
    numtets[ 5] = 26;
    numtets[ 6] = 30;
    numtets[ 7] = 34;
    numtets[ 8] = 38;
    numtets[9] = 44;
    numtets[10] = 48;
    numtets[11] = 52;
    numtets[12] = 56;
    numtets[13] = 60;
    numtets[14] = 60;
    numtets[15] = 52;
    numtets[16] =  4;
    numtets[17] =  2;
  }else if(nlambs == 27){
		 numtets[ 0] =  6;
         numtets[ 1] = 10;
         numtets[ 2] = 14;
         numtets[ 3] = 18;
         numtets[ 4] = 22;
         numtets[ 5] = 26;
         numtets[ 6] = 30;
         numtets[ 7] = 34;
         numtets[ 8] = 40;
         numtets[9] = 44;
         numtets[10] = 48;
         numtets[11] = 54;
         numtets[12] = 58;
         numtets[13] = 62;
         numtets[14] = 66;
         numtets[15] = 72;
         numtets[16] = 76;
         numtets[17] = 80;
         numtets[18] = 84;
         numtets[19] = 90;
         numtets[20] = 92;
         numtets[21] = 96;
         numtets[22] = 98;
         numtets[23] = 96;
         numtets[24] = 76;
         numtets[25] =  4;
         numtets[26] =  2;
  
  
  } else if(nlambs ==36){
		 numtets[ 0] =  8;
         numtets[ 1] = 12;
         numtets[ 2] = 16;
         numtets[ 3] = 20;
         numtets[ 4] = 24;
         numtets[ 5] = 28;
         numtets[ 6] = 32;
         numtets[ 7] = 36;
         numtets[ 8] = 40;
         numtets[9] = 44;
         numtets[10] = 48;
         numtets[11] = 54;
         numtets[12] = 58;
         numtets[13] = 62;
         numtets[14] = 68;
         numtets[15] = 72;
         numtets[16] = 76;
         numtets[17] = 82;
         numtets[18] = 86;
         numtets[19] = 90;
         numtets[20] = 96;
         numtets[21] =100;
         numtets[22] =104;
         numtets[23] =108;
         numtets[24] =114;
         numtets[25] =118;
         numtets[26] =122;
         numtets[27] =126;
         numtets[28] =130;
         numtets[29] =134;
         numtets[30] =138;
         numtets[31] =138;
         numtets[32] =128;
         numtets[33] =104;
         numtets[34] =  4;
         numtets[35] =  2;
  
  
  }
}

void ROTZ2X(const dcomplex *mpole, dcomplex *mrotate, const int direction)
{
  double *RD = (direction > 0 ? RDPLUS : RDMINUS);
  int ell, m, mp, offset1;

  offset1 = PTERMS*PGSZ;

  for ( m = 0; m <= PTERMS; m++ ) {
    int offset2 = m*(PTERMS+1);
    int offset3 = m*PGSZ+offset1;
    int offset4 = -m*PGSZ+offset1;
    for ( ell = m; ell <= PTERMS; ell++ ) {
      mrotate[ell + offset2] = mpole[ell]*RD[ell+offset3]; 
      for ( mp = 1; mp<= ell; mp++ ) {
	int offset5 = mp*(PTERMS+1);
	mrotate[ell+offset2]  += mpole[ell+offset5]*RD[ell+offset5+offset3]
	  + conj( mpole[ell+offset5] )*RD[ell + offset5 + offset4];
      }
    }
  }
}

void ROTZ2Y(const dcomplex *mpole, dcomplex *mrotate)
{
  dcomplex *mwork = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(1+PTERMS, sizeof(dcomplex));
  if ( mwork==0 || ephi==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  int m, ell, mp;
  ephi[0] = 1.0;
  for ( m =1; m <= PTERMS; m++) {
    ephi[m] = -ephi[m-1]*_Complex_I;
  }

  for ( m = 0; m <= PTERMS; m++ ) {
    int offset = m*(PTERMS+1);
    for ( ell = m; ell <= PTERMS; ell++) {
      int index = offset+ell;
      mwork[index] = ephi[m]*mpole[index];
    }
  }

  for ( m = 0; m <= PTERMS; m++ ) {
    int offset = m*(PTERMS+1);
    for ( ell = m; ell <= PTERMS; ell++ ) {
      int index = ell+offset;
      mrotate[index] = mwork[ell]*RDMINUS[ell+(m+PTERMS)*PGSZ];
      for ( mp = 1; mp <= ell; mp++ ) {
	int index1 = ell+mp*(PTERMS+1);
	mrotate[index] += mwork[index1]*RDMINUS[ell+mp*(PTERMS+1)+(m+PTERMS)*PGSZ]+
	  conj( mwork[index1] )*RDMINUS[ell+mp*(PTERMS+1)+(-m+PTERMS)*PGSZ];
      }
    }
  }

  free (ephi);
  free (mwork);
}

void ROTY2Z(const dcomplex *mpole, dcomplex *mrotate)
{
  int ell, m, mp;

  dcomplex *mwork = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(1+PTERMS, sizeof(dcomplex));
  if ( mwork==0 || ephi==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }
  ephi[0] = 1.0;
  for ( m = 1; m<= PTERMS; m++) {
    ephi[m] = ephi[m-1]*_Complex_I;
  }

  for ( m = 0; m <= PTERMS; m++ ) {
    int offset = m*(PTERMS+1);
    for ( ell = m; ell <= PTERMS; ell++ ) {
      int index = ell+offset;
      mwork[index] = mpole[ell]*RDPLUS[ell+(m+PTERMS)*PGSZ];
      for ( mp = 1; mp <= ell; mp++ ) {
	int index1 = ell+mp*(PTERMS+1);
	mwork[index] += mpole[index1]*RDPLUS[ell+mp*(PTERMS+1)+(m+PTERMS)*PGSZ]+
	  conj( mpole[index1] )*RDPLUS[ell+mp*(PTERMS+1)+(-m+PTERMS)*PGSZ];
      }
    }
  }

  
  for ( m = 0; m <= PTERMS; m++ ) {
    int offset = m*(PTERMS+1);
    for ( ell = m; ell <= PTERMS; ell++) {
      int index = ell+offset;
      mrotate[index] = ephi[m]*mwork[index];
    }
  }

  free (ephi);
  free (mwork);
}
