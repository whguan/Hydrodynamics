/*
  laplace_op.c: implementation of translation operators for fmm-laplace
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
#include <complex.h>
#include <math.h>
#include "fmm_ds.h"
#include "fmm_utils.h"
#include "fmm_global.h"
#include "fmm_laplace.h"

void kernel(const double *T, const double *S, const double charge, 
	    double *pot, double *field1, double *field2, double *field3, 
	    double *DXX, double *DYY, double *DZZ, double *DXY, double *DXZ,double *DYZ) 
{
  double rx = T[0] - S[0], ry = T[1] - S[1], rz = T[2] - S[2];
  double rr = rx*rx+ry*ry+rz*rz;
  double rdis = sqrt(rr);

  *pot = charge/rdis;
  double rmul = *pot/rr;

  *field1 = rmul*rx;
  *field2 = rmul*ry;
  *field3 = rmul*rz;
  
  *DXX = rmul * (2*rx*rx - ry*ry - rz*rz)/rr;
  *DYY = rmul * (2*ry*ry - rx*rx - rz*rz)/rr;
  *DZZ = rmul * (2*rz*rz - rx*rx - ry*ry)/rr;
  *DXY = rmul * (3*rx*ry)/rr;
  *DXZ = rmul * (3*rx*rz)/rr;
  *DYZ = rmul * (3*ry*rz)/rr;
}

void TSM(const int ibox)
{
  // TSM generates the multipole expansion for ibox using its particle information
  int ptr = BOXES[ibox].addr-1;
  int level = BOXES[ibox].level;
  dcomplex *mpole = &MPOLE[PGSZ*ibox];
  double *x0y0z0 = BOXES[ibox].center;
  double *ploc = &FMMLOC[3*ptr];
  double *pcharge = &FMMCHARGE[ptr];
  int nparts = BOXES[ibox].npts;
  double scale = SCALE[level];

  const double precis = 1.0e-14;
  
  double *powers = (double *)calloc(PTERMS+1, sizeof(double));
  double *p = (double *)calloc(PGSZ, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(PTERMS+1, sizeof(dcomplex));
  if ( powers==0 || p==0 || ephi == 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  int i;
  for ( i = 0; i < nparts; i++ ) {
    int ell, m;

    double rx = ploc[3*i] - x0y0z0[0];
    double ry = ploc[3*i+1] - x0y0z0[1];
    double rz = ploc[3*i+2] - x0y0z0[2];
    double proj = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = ( d <= precis ? 1.0 : rz/d );    
    ephi[0] = ( proj <= precis*d ? 1.0: rx/proj + _Complex_I*ry/proj );

    d = d*scale;
    powers[0] = 1.0;
    for ( ell = 1; ell <= PTERMS; ell++) {
      powers[ell] = powers[ell-1]*d;
      ephi[ell] = ephi[ell-1]*ephi[0];
    }

    mpole[0] += pcharge[i];

    lgndr (PTERMS, ctheta, p);

    for ( ell = 1; ell <= PTERMS; ell++) {
      double cp = pcharge[i]*powers[ell]*p[ell];
      mpole[ell] += cp;
    }

    for ( m = 1; m <= PTERMS; m++) {
      int offset = m*(PTERMS+1);
      for ( ell = m; ell <= PTERMS; ell++ ) {
	double cp = pcharge[i]*powers[ell]*YTOPC[ell+m*61]*p[ell+offset];
	mpole[ell+offset] += cp*conj(ephi[m-1]);
      }
    }
  }
  free (powers);
  free (ephi);
  free (p);
}

void TMM(const int pbox)
{
  // TMM computes the multipole expansion of a nonleaf box by translating
  // and accumulating each of its child's multipole expansion

  static dcomplex var[5] = {1, -1+_Complex_I, 1+_Complex_I, 1-_Complex_I, -1-_Complex_I};
  const double arg = sqrt(2)/2.0;

  double *powers = (double *)calloc(PTERMS+3, sizeof(double));
  dcomplex *mpolen = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(PTERMS+3, sizeof(dcomplex));
  if ( powers==0 || mpolen==0 || marray==0 || ephi==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  dcomplex *pMpole = &MPOLE[PGSZ*pbox];
  const int lev = BOXES[pbox].level;
  const double sc1 = SCALE[lev+1];
  const double sc2 = SCALE[lev];

  int i;
  for ( i = 0; i < 8; i++ ) {
    int cbox = BOXES[pbox].child[i];
    if ( cbox > 0 ) {      
      int ell, m, j, k, mp;

      int ifl = IFL_UP[i];
      double *RD = ( i < 4 ? RDSQ3: RDMSQ3);
      dcomplex *cMpole = &MPOLE[PGSZ*cbox];
      double dd = -sqrt(3)/2.0;

      ephi[0] = 1.0;
      ephi[1] = arg*var[ifl];
      powers[0] = 1.0;
      for ( ell = 1; ell <= PTERMS+1; ell++ ) {
		powers[ell] = powers[ell-1]*dd;
		ephi[ell+1] = ephi[ell]*ephi[1];
      }

      // a rotation of phi radians about the z-axis in the original coordinate system.
      for ( m = 0; m <= PTERMS; m++ ) {
	  int offset = m*(PTERMS+1);
	  for ( ell = m; ell <= PTERMS; ell++) {
	      int index = ell+offset; 
	      mpolen[index] = conj( ephi[m] )*cMpole[index];
	  }
      }

      // a rotation about the y'-axis in the rotated system.
      for ( m = 0; m <= PTERMS; m++) {
		int offset = m*(PTERMS+1);
		int offset1 = (m+PTERMS)*PGSZ;
		int offset2= (-m+PTERMS)*PGSZ;
		for ( ell = m; ell <= PTERMS; ell++) {
			int index = offset+ell;
			marray[index] = mpolen[ell]*RD[ell+offset1];
			for ( mp = 1; mp <= ell; mp++ ) {
				int index1 = ell+mp*(PTERMS+1);
				marray[index] += mpolen[index1]*RD[index1+offset1]+
				conj( mpolen[index1] )*RD[index1+offset2];
			}	
		}
      }
      
      // shift along z-axis
      for ( k = 0; k <= PTERMS; k++) {
		int offset = k*(PTERMS+1);
		for ( j = k; j<=PTERMS; j++) {
			int index = j+offset;
			mpolen[index] = marray[index];
			for ( ell = 1; ell <= j-k; ell++) {
				int index2 = j-k+ell*(2*PTERMS+1);
				int index3 = j+k+ell*(2*PTERMS+1);
				mpolen[index] += marray[index-ell]*powers[ell]*DC[index2]*DC[index3];
			}
		}
      }
      
      // reverse rotation about the y'-axis
      for ( m= 0; m <= PTERMS; m=m+2) {
		int offset = m*(PTERMS+1);
		int offset1 = (m+PTERMS)*PGSZ;
		int offset2 = (-m+PTERMS)*PGSZ;
		for ( ell = m; ell<=PTERMS; ell++) {
			int index = ell+offset;
			marray[index] = mpolen[ell]*RD[ell+offset1];
			for ( mp = 1; mp <= ell; mp = mp+2) {
				int index1 = ell+mp*(PTERMS+1);
				marray[index] = marray[index]- mpolen[index1]*RD[index1+offset1]-
					conj( mpolen[index1] )*RD[index1+offset2];
			}
			for ( mp = 2; mp <= ell; mp=mp+2) {
				int index1 = ell+mp*(PTERMS+1);
				marray[index] = marray[index] + mpolen[index1]*RD[index1+offset1]+
					conj( mpolen[index1] )*RD[index1+offset2];
			}
		}
      }
      
      for ( m= 1; m <= PTERMS; m=m+2) {
		int offset = m*(PTERMS+1);
		int offset1 = (m+PTERMS)*PGSZ;
		int offset2 = (-m+PTERMS)*PGSZ;
		for ( ell = m; ell <= PTERMS; ell++ ) {
			int index = ell+offset;
			marray[index] = -mpolen[ell]*RD[ell+offset1];
			for ( mp = 1; mp <= ell; mp=mp+2 ) {
				int index1 = ell+mp*(PTERMS+1);
				marray[index] = marray[index]+ mpolen[index1]*RD[index1+offset1] + 
				conj( mpolen[index1] )*RD[index1+offset2];
			}
			for ( mp = 2; mp <= ell; mp=mp+2) {
				int index1 = ell+mp*(PTERMS+1);
				marray[index] = marray[index] - mpolen[index1]*RD[index1+offset1]-
					conj( mpolen[index1] )*RD[index1+offset2];
			}
		}	
      }

      // rotate back phi radians about the z-axis in the above system and rescale
      powers[0] = 1.0;
      dd = sc2/sc1;
      for ( ell = 1; ell <= PTERMS+1; ell++ ) {
	     powers[ell] = powers[ell-1]*dd;
      }

      for ( m= 0; m <= PTERMS; m++ ) {
	int offset = m*(PTERMS+1);
	for ( ell = m; ell <= PTERMS; ell++ ) {
	  int index = ell+offset;
	  mpolen[index] = ephi[m] * marray[index]*powers[ell];
	}
      }

      for ( m = 0; m < PGSZ; m++) 
	pMpole[m] += mpolen[m];

    }
  }

  free (ephi);
  free (powers);
  free (mpolen);
  free (marray);
}

void TME(const int ibox)
{
  // TME converts the multipole expansion of ibox into six exponential expansions
  dcomplex *mw = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *mexpf1 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf2 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  if ( mw==0 || mexpf1==0 || mexpf2== 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }
  
  TME_PHASE1(&MPOLE[PGSZ*ibox], mexpf1, mexpf2, ibox);
  TME_PHASE2(mexpf1, &LEXPU[NEXPTOTP*ibox], ibox);
  TME_PHASE2(mexpf2, &LEXPD[NEXPTOTP*ibox], ibox);
  
  ROTZ2Y(&MPOLE[PGSZ*ibox], mw);
  TME_PHASE1(mw, mexpf1, mexpf2, ibox);
  TME_PHASE2(mexpf1, &LEXPN[NEXPTOTP*ibox], ibox);
  TME_PHASE2(mexpf2, &LEXPS[NEXPTOTP*ibox], ibox);
  
  ROTZ2X(&MPOLE[PGSZ*ibox], mw, 1);
  TME_PHASE1(mw, mexpf1, mexpf2, ibox);
  TME_PHASE2(mexpf1, &LEXPE[NEXPTOTP*ibox], ibox);
  TME_PHASE2(mexpf2, &LEXPW[NEXPTOTP*ibox], ibox);
  
  free(mw);
  free(mexpf1);
  free(mexpf2);
}

void TME_PHASE1(const dcomplex *mpole, dcomplex *mexpup, dcomplex *mexpdown, const int ibox)
{
  int nell, mth, ncurrent, nm;
  dcomplex zeyep, ztmp1, ztmp2;
  int ntot = 0;
  for ( nell = 0; nell < NLAMBS; nell++ ) {
    double sgn = -1.0;
    zeyep = 1.0;
    for ( mth = 0; mth <= NUMFOUR[nell]-1; mth++ ) {
      ncurrent = ntot + mth ;
      ztmp1 = 0.0;
      ztmp2 = 0.0;
      sgn = -sgn;
      int offset = mth*(PTERMS+1);
      int offset1 = offset+nell*PGSZ;
      for ( nm = mth; nm <= PTERMS; nm=nm+2 ) {
	ztmp1 += RLSC[nm+offset1]*mpole[nm+offset];
      }
      for ( nm = mth+1; nm<=PTERMS; nm=nm+2) {
	ztmp2 += RLSC[nm+offset1]*mpole[nm+offset];
      }
      mexpup[ncurrent] = (ztmp1+ztmp2)*zeyep;
      mexpdown[ncurrent] = sgn*(ztmp1 - ztmp2)*zeyep;
      zeyep = zeyep*_Complex_I;
    }
    ntot = ntot + NUMFOUR[nell];
  }
}

void TME_PHASE2(const dcomplex *mexpf, dcomplex *mexpphys, const int ibox)
{
  int i;
  int nftot = 0;
  int nptot = 0;
  int nexte = 0;
  int nexto = 0;

  for ( i = 0; i < NLAMBS; i++) {
    int ival;
    for ( ival = 0; ival < NUMPHYS[i]/2; ival++) {
      int nm;
      mexpphys[nptot+ival] = mexpf[nftot];
      for ( nm = 1; nm < NUMFOUR[i]; nm=nm+2) {
	double rt1 = cimag( FEXPE[nexte] )* creal( mexpf[nftot+nm] );
	double rt2 = creal( FEXPE[nexte] )*cimag( mexpf[nftot+nm] );
	double rtmp = 2*(rt1+rt2);
	nexte += 1;
	mexpphys[nptot+ival] += rtmp*_Complex_I;
      }

      for ( nm=2; nm < NUMFOUR[i]; nm=nm+2) {
	double rt1 = creal( FEXPO[nexto] )*creal( mexpf[nftot+nm] );
	double rt2 = cimag( FEXPO[nexto] )*cimag( mexpf[nftot+nm] );
	double rtmp = 2*(rt1-rt2);
	nexto += 1;
	mexpphys[nptot+ival] += rtmp;
      }
    }

    nftot += NUMFOUR[i];
    nptot += NUMPHYS[i]/2;
  }
}

void TSL(const int ibox, const int member)
{
  // TSL computes a local expansion about ibox's center using particle information contained
  // in the box member.

  double *x0y0z0 = BOXES[ibox].center;
  int ptr = BOXES[member].addr-1;
  double *ploc = &FMMLOC[3*ptr];
  double *pcharge = &FMMCHARGE[ptr];
  int nparts = BOXES[member].npts;
  dcomplex *local = &LOCAL[PGSZ*ibox];
  double scale = SCALE[BOXES[ibox].level];

  double *powers = (double *)calloc(PTERMS+3, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(PTERMS+2, sizeof(dcomplex));
  double *p = (double *)calloc(PGSZ, sizeof(double));
  if ( powers==0 || ephi==0 || p==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n",__FILE__,__LINE__);
    exit(-1);
  }

  int i;
  const double precis = 1.0e-14;

  for ( i = 0; i < nparts; i++ ) {
    int ell, m;

    double rx = ploc[3*i] - x0y0z0[0];
    double ry = ploc[3*i+1] - x0y0z0[1];
    double rz = ploc[3*i+2] - x0y0z0[2];
    double proj = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
    double d = sqrt(rr);
    double ctheta = ( d <= precis ? 1.0 : rz/d );
    ephi[0] = ( proj <= precis*d ? 1.0: rx/proj - _Complex_I*ry/proj );
    
    d = 1.0/d;
    powers[0] = 1.0;
    powers[1] = d;
    d = d/scale;

    for ( ell = 2; ell <= PTERMS+2; ell++ ) {
      powers[ell] = powers[ell-1]*d;
    }

    for ( ell = 1; ell <= PTERMS+1; ell++) {
      ephi[ell] = ephi[ell-1]*ephi[0];
    }

    local[0] += pcharge[i]*powers[1];

    lgndr (PTERMS, ctheta, p);

    for ( ell = 1; ell <= PTERMS; ell++ ) {
      double cp = pcharge[i]*p[ell]*powers[ell+1];
      local[ell] += cp;
    }

    for ( m = 1; m <= PTERMS; m++ ) {
      int offset = m*(PTERMS+1); 
      for ( ell = m; ell <= PTERMS; ell++ ) {
	double cp = pcharge[i]*powers[ell+1]*YTOPC[ell+m*61]*p[ell+offset];
	local[ell+offset] += cp*ephi[m-1];
      }
    }
  }

  free (powers);
  free (ephi);
  free (p);
}

void TLT(const int ibox, const int ifpot, const int ifgrad, const int ifhess)
{
  // TLT evaluates the local expansion at each particle location contained in a leaf box
  dcomplex *local  = &LOCAL[PGSZ*ibox];
  double   *x0y0z0 = BOXES[ibox].center;
  double    scale  = SCALE[BOXES[ibox].level];

  double *p = (double *)calloc(PGSZ, sizeof(double));
  double *powers = (double *)calloc(1+PTERMS, sizeof(double));
  dcomplex *ephi = (dcomplex *)calloc(1+PTERMS, sizeof(dcomplex));
  dcomplex *ylm = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  
  if ( p==0 || powers==0 || ephi==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  const double precis = 1.0e-14;

  int i;
  for ( i = 0; i < BOXES[ibox].npts; i++ ) {    
    int ell, m;
    
    int     ptr = BOXES[ibox].addr+i-1;
    double *point = &FMMLOC[3*ptr];
    double *pot = &FMMPOT[ptr];
    double *field = &FMMFIELD[3*ptr];
	
    // For second order derivative 
    double *dxx = &FMMDXX[ptr];
    double *dyy = &FMMDYY[ptr];
    double *dzz = &FMMDZZ[ptr]; 
    double *dxy = &FMMDXY[ptr];
    double *dxz = &FMMDXZ[ptr];
    double *dyz = &FMMDYZ[ptr];
    
    double rpotz = 0.0;
	
    dcomplex zs1 = 0.0;
    dcomplex zs2 = 0.0;
        
    double field0 = 0.0;
    double field1 = 0.0;
    double field2 = 0.0;

    double rx = point[0] - x0y0z0[0];
    double ry = point[1] - x0y0z0[1];
    double rz = point[2] - x0y0z0[2];
	
    double proj  = rx*rx + ry*ry;
    double rr = proj + rz*rz;
    proj = sqrt(proj);
	
    double d = sqrt(rr);
	
    double ctheta = ( d<=precis ? 0.0 : rz/d );
    ephi[0] = ( proj <= precis*d ? 1.0: rx/proj + _Complex_I*ry/proj );
    d = d*scale;
	
    double dd = d;
    powers[0] = 1.0;
	
    // Compute r^l = powers[ell] and exp(im\phi) = ephi[m-1]
    for ( ell = 1; ell <= PTERMS; ell++) {
      powers[ell] = dd;
      dd = dd*d;
      ephi[ell] = ephi[ell-1]*ephi[0];
    }
    
    //  Compute Legendre Polynomial p[ell+m*(pterms +1)]
    lgndr (PTERMS, ctheta, p);
    
    //**********************************************************************************
    // Spherical harmonic:  Compute Ylm
    // Note here, local coefficients contain the term sqrt((2*l+1)/(4pi))in 
    // orginal code. So here ylm below is computed without the 1/sqrt(4pi) term.
    //*********************************************************************************	
    ylm[0] = 1.0;
   	
    for ( ell = 1; ell <= PTERMS; ell++) {
	   	ylm[ell] = sqrt((2*ell+1))*p[ell];	
    }
    	
    for ( ell = 1; ell <= PTERMS; ell++) {
	double coeff = sqrt((2*ell+1));
	
	for (m =1; m<= ell; m++){
        	ylm[ell+m*(PTERMS+1)] = coeff * YTOPC[ell+m*61]*p[ell+m*(PTERMS+1)]* ephi[m-1];
	}
    }
    /***************************************************************************************	   
      Potential
      *************************************************************************************/
  if(ifpot!=0){
    //potential m = 0
    for (ell = 0; ell<= PTERMS; ell++){
	*pot += local[ell]*ylm[ell]*powers[ell]/sqrt(2*ell+1);
    }
	
    //potential m != 0
    for (ell = 1; ell <= PTERMS; ell++ ){
	for(m = 1; m <= ell; m++ ){
	    *pot += 2.0 *creal(local[ell+m*(PTERMS+1)]/sqrt(2*ell+1)* ylm[ell+m*(PTERMS+1)])* powers[ell];
	}
	
    }
  }
    /***************************************************************************************
      Field
    ****************************************************************************************/
    //Field D0 m = 0, zs1 =  \sum_{l=0}^{p-1} \sum{m=1}^{l} A(l,m) 
    //                zs2 =  \sum_{l=0}^{p-1} \sum{m=0}^{l} C(l,m) in notes 
    // bug fixed in March 16, 2015
  if(ifgrad!=0){ 
    for (ell = 0; ell<= PTERMS-1; ell++){
	double alpha0 = Alpha(ell+1,0); 
	field2 += alpha0* sqrt(2*ell+3)* local[ell+1]*ylm[ell]*powers[ell];
	double beta0 = Beta(ell+1, 1,-1);
	zs2 += beta0* sqrt(2*ell+3)*local[ell+1+PTERMS+1]*ylm[ell]*powers[ell];
    }
	
    // Field D0 m!=0
    for (ell = 1; ell <= PTERMS-1; ell++ ){
	for(m = 1; m <= ell; m++ ){
		dcomplex coeff = sqrt(2*ell+3)* ylm[ell+m*(PTERMS+1)]*powers[ell];
		double alpha = Alpha(ell+1,m);
		field2 += 2.0 *alpha*creal(local[ell+1+m*(PTERMS+1)]*coeff);
		double beta1 = Beta(ell+1,m-1, 1);
		double beta2 = Beta(ell+1,m+1,-1);
	        //Field D+, D-	
		zs1 += beta1*local[ell+1+(m-1)*(PTERMS+1)] *coeff;
		zs2 += beta2*local[ell+1+(m+1)*(PTERMS+1)] *coeff;
	}
    }
			
    field0 =  sqrt(2)* creal( zs1 - zs2 );
    field1 =  sqrt(2)* cimag( zs1 + zs2 );
    /* Here, notice that the 
      field =[-\partial{\phi}/\partial{x}, -\partial{\phi}/\partial{y}, -\partial{\phi}/\partial{z}] 
      
      field0 = \partial{\phi}/\partial{x}
      field1 = \partial{\phi}/\partial{y}
      field2 = \partial{\phi}/\partial{z}
    so we have to use -= below */

    field[0] -= field0*scale;
    field[1] -= field1*scale;
    field[2] -= field2*scale;
 }   
/***************************************************************************
// Below is to compute the second order derivative
// By Wenhua Guan
// Sep 12th, 2014 
****************************************************************************/
 if(ifhess!=0){
  // DPP, DMM, DPM 
	dcomplex DPP = 0.0;
	dcomplex DMM = 0.0;
	dcomplex DP0 = 0.0;
	dcomplex DM0 = 0.0;
	double D00 = 0.0;
	double DPM = 0.0;
	
	// Compute the second order derivative
        //1 . Use DPP, DPM, DMM to compute dxx, dyy. 
	for(ell = 0; ell <= PTERMS-2; ell++){
		double coeff = (2*ell+3)*(2*ell+5)*powers[ell];   // common coeff: (2l+3)(2l+5)r^l
		
		// First compute all the m = 0 case
		// -A(ell,1) is included in DPP,
		// L_{ell+2}^{-1}is not stored in system, obtained by L_{ell+2}^{1} 
		
		double beta0 = Beta(ell+2, -1, 1)* Beta(ell+1, 0, 1); 
				
		//DPP m =1
		dcomplex Local0 = creal(local[ell+2+(PTERMS+1)])-_Complex_I * cimag(local[ell+2+(PTERMS+1)]);
		DPP += -coeff * beta0 * Local0 * ylm[ell+(PTERMS+1)]/sqrt(2*(ell+2)+1);
					
		//DPM m =0
		double beta1 = Beta(ell+2,0,1)*Beta(ell+1,1,-1);
		DPM += coeff * beta1 * creal(local[ell+2]) * ylm[ell]/sqrt(2*(ell+2)+1);
		
		//D00 m =0
		double alpha0 = Alpha(ell+2, 0)*Alpha(ell+1,0);
		D00 += coeff * alpha0 * creal(local[ell+2]) * ylm[ell]/sqrt(2*(ell+2)+1);
		
		// DPP m!=1
		for(m = 2; m <= ell; m++){		
			double betaterm = Beta(ell+2, m-2, 1) * Beta(ell+1, m-1,1);
			dcomplex Alm = local[ell+2+(m-2)*(PTERMS+1)]*ylm[ell+m*(PTERMS+1)];
			DPP +=  coeff* betaterm * Alm/sqrt(2*(ell+2)+1);
			
		}
				
		// DMM
		for(m = 0; m <= ell; m++){
			double betaterm = Beta(ell+2, m+2, -1) * Beta(ell+1, m+1, -1);
			dcomplex Clm = local[ell+2+(m+2)*(PTERMS+1)]*ylm[ell+m*(PTERMS+1)];
			DMM +=  coeff * betaterm * Clm/sqrt(2*(ell+2)+1);			
		}
		
		
		//DPM  
		for(m = 1; m <= ell; m++){
			double betaterm = Beta(ell+2, m, 1) * Beta(ell+1, m+1, -1);
			dcomplex Elm = local[ell+2+m*(PTERMS+1)] * ylm[ell+m*(PTERMS+1)];
			DPM += 2.0 *coeff * betaterm * creal(Elm)/sqrt(2*(ell+2)+1);		
		}	
		
		//D00
		for(m = 1; m <= ell; m++){
			double alphaterm = Alpha(ell+2,m) * Alpha(ell+1,m);
			dcomplex Glm = local[ell+2+m*(PTERMS+1)] * ylm[ell+m*(PTERMS+1)];
			D00 += 2.0 * coeff * alphaterm * creal(Glm)/sqrt(2*(ell+2)+1); 
		}
		
		//DP0
		for(m = 1; m<=ell; m++){
			double alphabeta = Alpha(ell+2,m-1)*Beta(ell+1,m-1,1);
			dcomplex Plm1 = local[ell+2+(m-1)*(PTERMS +1)]* ylm[ell+m*(PTERMS+1)];
			DP0 += coeff * alphabeta * Plm1/sqrt(2*(ell+2)+1);
		}
				
		//DM0		
		for(m = 0; m<=ell; m++){
			double alphabeta = Alpha(ell+2,m+1)*Beta(ell+1,m+1,-1);
			dcomplex Plm = local[ell+2+(m+1)*(PTERMS +1)]* ylm[ell+m*(PTERMS+1)];
			DM0 += coeff * alphabeta * Plm/sqrt(2*(ell+2)+1);
		}			
	}

	*dxx +=  pow(scale,2) * (creal(DPP) + creal(DMM) - DPM);
	*dyy += -pow(scale,2) * (creal(DPP) + creal(DMM) + DPM);
	*dzz += pow(scale,2) *  D00;
	*dxy += pow(scale,2) * ( cimag(DPP) - cimag(DMM) );
	*dxz += pow(scale,2) * sqrt(2)*( creal(DP0) - creal(DM0));
	*dyz += pow(scale,2) * sqrt(2)*( cimag(DP0) + cimag(DM0));
    }
  
  }
 
  free (ylm);
  free (ephi);
  free (powers);
  free (p);

}

// Alpha returns the coefficient of sqrt((l^2-m^2)/(4*l^2-1))
double Alpha(int l, int m )
{
  double alm;
  alm = (l+m)*(l-m);
  alm /= ((2*l+1)*(2*l-1));
  alm = sqrt(alm);
  return alm;
 }
 
//Beta returns the coefficients
//Beta_+(l,m) = sqrt((l-m-1)(l-m)/(2(2l-1)(2l+1)))
//Beta_-(l,m) = sqrt((l+m-1)(l+m)/(2(2l-1)(2l+1)))
double Beta(int l, int m, int option)
{
	double blm = 0.0;
	double ans1 = 0.0;
	double ans2 = 0.0; 
	ans2 =(double)(2*(2*l-1)*(2*l+1));
	
	if (option == 1){
		ans1 = (double)((l-m)*(l-m-1));
		blm  = ans1/ans2;
		blm  = sqrt(blm);
	}else if (option == -1){
		ans1 = (double)((l+m)*(l+m-1));
		blm  = ans1/ans2;
		blm  = sqrt(blm);
	}else{
		printf("Error in the value of option: %d\n", option);
		exit(-1);
	}
	
	return blm;
}

void TMT(const int ibox, const int i, const int member)
{
  dcomplex *mpole = &MPOLE[PGSZ*member];
  double *x0y0z0 = BOXES[member].center;
  int ptr = BOXES[ibox].addr-1;
  double *point = &FMMLOC[3*(ptr+i)];
  double *pot = &FMMPOT[ptr+i];
  double *field = &FMMFIELD[3*(ptr+i)];
  double scale = SCALE[BOXES[member].level];
 
  double *p = (double *)calloc((PTERMS+2)*(PTERMS+2), sizeof(double));
  double *powers = (double *)calloc(PTERMS+4, sizeof(double));
  dcomplex *ephi= (dcomplex *)calloc(PTERMS+3, sizeof(dcomplex));
  if ( p==0 || powers==0 || ephi== 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  const double precis = 1.0e-14;
  int ell, m;

  double rpotz = 0.0;
  double field1 = 0;
  double field2 = 0; 
  double field3 = 0;
  dcomplex zs1 = 0.0;
  dcomplex zs2 = 0.0;
  dcomplex zs3 = 0.0;
  double rx = point[0] - x0y0z0[0];
  double ry = point[1] - x0y0z0[1];
  double rz = point[2] - x0y0z0[2];
  double proj = rx*rx + ry*ry;
  double rr = proj + rz*rz;
  proj = sqrt(proj);

  double d = sqrt(rr);
  double ctheta = ( d<=precis ? 0.0: rz/d );
  ephi[0] = ( proj<= precis*d ? 1.0: rx/proj+_Complex_I*ry/proj );

  d = 1.0/d;
  powers[0] = 1.0;
  powers[1] = d;
  d = d/scale;
  for ( ell = 1; ell <= PTERMS+2; ell++ ) {
    powers[ell+1] = powers[ell]*d;
    ephi[ell] = ephi[ell-1]*ephi[0];
  }

  lgndr (PTERMS+1, ctheta, p);

  double rtemp;

  // add contributions from monopole moment
  double rmp = creal( mpole[0] );
  *pot = *pot + rmp*powers[1];
  dcomplex cpz = ephi[0]*rmp*powers[2]*YTOPC[1+61]*p[1+PTERMS+2]*YTOPCSINV[1+61];
  zs1 = zs1 + cpz;
  double cp = rmp*powers[2]*p[1]*YTOPCS[0]*YTOPCSINV[1];
  field3 = cp;

  
  // add contributions from legendre polynomials
  for ( ell = 1; ell <= PTERMS; ell++ ) {
    rmp = creal ( mpole[ell] );
    cp = rmp*powers[ell+1]*p[ell];
    *pot = *pot + cp;
    zs1 += ephi[0]*rmp*powers[ell+2]*YTOPC[ell+1+61]*p[ell+1+PTERMS+2]*YTOPCS[ell]*YTOPCSINV[ell+1+61];
    cpz = mpole[ell+PTERMS+1];
    rtemp = powers[ell+2]*p[ell+1]*YTOPCSINV[ell+1];
    zs2 = zs2 + cpz*rtemp*YTOPCS[ell+61];
    cp = rmp*rtemp*YTOPCS[ell];
    field3 = field3 + cp;
  }

 
  // add contributions from associated legendre functions
  for ( ell = 1; ell <= PTERMS; ell++ ) {
    for ( m = 1; m <= ell; m++ ) {
      cpz = mpole[ell+m*(PTERMS+1)]*powers[ell+1]*YTOPC[ell+m*61]*p[ell+m*(PTERMS+2)];
      rpotz += creal(cpz*ephi[m-1]);
      cpz = mpole[ell+m*(PTERMS+1)]*YTOPCS[ell+m*61]*powers[ell+2];
      zs1 += cpz*ephi[m]*YTOPCSINV[ell+1+(m+1)*61]*YTOPC[ell+1+(m+1)*61]*p[ell+1+(m+1)*(PTERMS+2)];
      if ( m > 1 ) {
	zs2 += cpz*ephi[m-2]*YTOPCSINV[ell+1+(m-1)*61]*YTOPC[ell+1+(m-1)*61]*p[ell+1+(m-1)*(PTERMS+2)];
      }
      zs3 += cpz*ephi[m-1]*YTOPC[ell+1+m*61]*p[ell+1+m*(PTERMS+2)]*YTOPCSINV[ell+1+m*61];
    }
  }
  

  *pot += 2.0*rpotz;
  field1 = creal ( zs2 - zs1);
  field2 = -cimag ( zs2+zs1);
  field3 = field3 + 2.0*creal(zs3);

  field[0] += field1*scale;
  field[1] += field2*scale;
  field[2] += field3*scale;
  

  free (ephi);
  free (powers);
  free (p);
}

void TEL(const int ibox) 
{
  // TEL converts exponential expansion into local expansions. The implementation applies
  // the merge-and-shift technique, which requires ibox to be a nonleaf box.
  if ( BOXES[ibox].nchild ) {
    TEL_UD (ibox);
    TEL_NS (ibox);
    TEL_EW (ibox);
  }
}

void TEL_UD(const int ibox)
{
  int *child = BOXES[ibox].child;
  int ilev = BOXES[ibox].level;
  double scale = SCALE[ilev+1];

  int iuall[36], nuall, ixuall[36], iyuall[36], iu1234[16], nu1234, ix1234[16], iy1234[16],
    idall[36], ndall, ixdall[36], iydall[36], id5678[16], nd5678, ix5678[16], iy5678[16];
  create_list_up (ibox, iuall, &nuall, ixuall, iyuall, iu1234, &nu1234, ix1234, iy1234);
  create_list_down (ibox, idall, &ndall, ixdall, iydall, id5678, &nd5678, ix5678, iy5678);

  dcomplex *mexuall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexu1234 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexdall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexd5678 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mw1 = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *temp = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf1 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf2 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  if ( mexuall==0 || mexu1234==0 || mexdall== 0 || mexd5678==0 ||
       mw1==0 || temp==0 || mexpf1==0 || mexpf2==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  int i, jj, offset, iexp1, iexp2;
  dcomplex zmul;

  //if exists, process upall list  
  if ( nuall > 0 ) {
    for ( i = 0; i < nuall; i++ ) {
      offset = iuall[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixuall[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ixuall[i]-1];

	if ( ixuall[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ixuall[i]-1]); 

	if ( iyuall[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iyuall[i]-1]; 

	if ( iyuall[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iyuall[i]-1]); 

	mexuall[jj] += zmul*LEXPD[offset+jj];
      }
    }
  }

  //if exists, process up1234 list  
  if ( nu1234 > 0 ) {
    for ( i = 0; i < nu1234; i++ ) {
      offset = iu1234[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix1234[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix1234[i]-1];

	if ( ix1234[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix1234[i]-1]); 

	if ( iy1234[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy1234[i]-1]; 

	if ( iy1234[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy1234[i]-1]); 

	mexu1234[jj] += zmul*LEXPD[offset+jj];
      }
    }
  }

  // if exists, process dall list  
  if ( ndall > 0 ) {
    for ( i = 0; i < ndall; i++ ) {
      offset = idall[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixdall[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ixdall[i]-1]); 

	if ( ixdall[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ixdall[i]-1];

	if ( iydall[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iydall[i]-1]);

	if ( iydall[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iydall[i]-1];

	mexdall[jj] += zmul*LEXPU[offset+jj];
      }
    }
  }

  // if exists, process d5678 list  
  if ( nd5678 > 0 ) {
    for ( i = 0; i < nd5678; i++ ) {
      offset = id5678[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix5678[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix5678[i]-1]);

	if ( ix5678[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix5678[i]-1];

	if ( iy5678[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy5678[i]-1]);

	if ( iy5678[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy5678[i]-1];

	mexd5678[jj] += zmul*LEXPU[offset+jj];
      }
    }
  }

  // process child #1 
  if ( child[0] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    iexp1 = 0;
    iexp2 = 0;

    // receive upall & up1234 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+2]*scale;
      iexp1 = iexp1+1;
    }

    if ( nu1234 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj]+mexu1234[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);
    
    // receive dall
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0); 
    }

    // convert into local expansion of child #1
    if ( iexp1 +iexp2 ) {
      TEL_PHASE2(mw1, iexp2, mexpf2, iexp1, mexpf1, 0);
      vecAdd(mw1, &LOCAL[child[0]*PGSZ], PGSZ);
    }
  }

  // process child #2
  if ( child[1] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive upall & up1234 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+2]*conj(XS[3*jj])*scale;   
      iexp1 = iexp1+1;
    }

    if ( nu1234 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexu1234[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;   
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive dall
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      TEL_PHASE1(mexpf2,temp,0);
      iexp2 = iexp2+1;
    }

    // convert into local expansion of child #2
    if ( iexp1 +iexp2) {
      TEL_PHASE2(mw1, iexp2, mexpf2, iexp1, mexpf1, 0);
      vecAdd(mw1, &LOCAL[child[1]*PGSZ], PGSZ);
    }
  }

  // process child #3
  if ( child[2] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive upall & up1234 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+2]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nu1234 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexu1234[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1,temp,0);
    
    // receive dall
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #3
    if ( iexp1 > 0 || iexp2 > 0 ) {
      TEL_PHASE2(mw1, iexp2, mexpf2, iexp1, mexpf1, 0);
      vecAdd(mw1, &LOCAL[child[2]*PGSZ], PGSZ);
    }
  }

  // process child #4
  if ( child[3] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive upall & up1234 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+2]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nu1234 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexu1234[jj]*ZS[3*jj+1]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive dall 
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+1]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #4
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1,  iexp2, mexpf2, iexp1, mexpf1, 0);
      vecAdd(mw1, &LOCAL[child[3]*PGSZ], PGSZ);
    }
  }

  // process child #5
  if ( child[4] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive upall 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0); 
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive dall & d5678 
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+2]*scale;
      iexp2 = iexp2+1;
    }

    if ( nd5678 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexd5678[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #5
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      vecAdd(mw1, &LOCAL[child[4]*PGSZ], PGSZ);
    }
  }

  // process child #6
  if ( child[5] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive upall 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1,temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive dall & d5678
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+2]*XS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nd5678 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexd5678[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #6
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1,  iexp2, mexpf2, iexp1, mexpf1,0);
      vecAdd(mw1, &LOCAL[child[5]*PGSZ], PGSZ);
    }
  }

  // process child #7
  if ( child[6] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive upall 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive dall & d5678
    if ( ndall > 0 )  {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+2]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nd5678 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexd5678[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #7
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      vecAdd(mw1, &LOCAL[child[6]*PGSZ], PGSZ);
    }
  }

  // process child #8
  if ( child[7] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive upall 
    if ( nuall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexuall[jj]*ZS[3*jj+1]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive dall & d5678
    if ( ndall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexdall[jj]*ZS[3*jj+2]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nd5678 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexd5678[jj]*ZS[3*jj+1]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #8
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1,  iexp2, mexpf2, iexp1, mexpf1,0);
      vecAdd(mw1, &LOCAL[child[7]*PGSZ], PGSZ);
    }
  }   

  free (mexuall);
  free (mexu1234);
  free (mexdall);
  free (mexd5678);
  free (mw1);
  free (temp);
  free (mexpf1);
  free (mexpf2);

  return;
}

/*
    after rotation, the layout of the boxes are the following 
       ^ z' 
       |             ---------
       |            / 4    8 /
       |           / 3    7 /
       |          ----------
       |  y'  
       |   /      ---------
       |  /      / 2    6 /
       | /      / 1    5 /
       |/      ----------            
----------------------------------> x' 
*/

void TEL_NS(const int ibox)
{
  int inall[36], nnall, ixnall[36], iynall[36], in1256[16], nn1256, ix1256[16], iy1256[16],
    in12[4], nn12, ix12[4], iy12[4], in56[4], nn56, ix56[4], iy56[4],
    isall[36], nsall, ixsall[36], iysall[36], is3478[16], ns3478, ix3478[16], iy3478[16],
    is34[4], ns34, ix34[4], iy34[4], is78[4], ns78, ix78[4], iy78[4];

  create_list_north (ibox, inall, &nnall, ixnall, iynall, in1256, &nn1256, ix1256, iy1256,
		     in12, &nn12, ix12, iy12, in56, &nn56, ix56, iy56);
  create_list_south (ibox, isall, &nsall, ixsall, iysall, is3478, &ns3478, ix3478, iy3478,
		     is34, &ns34, ix34, iy34, is78, &ns78, ix78, iy78);

  dcomplex *mexnall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexn1256 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexn12 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexn56 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexsall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexs3478 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexs34 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexs78 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mw1 = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *mw2 = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *temp =(dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf1 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf2 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  if ( mexnall==0 || mexn1256==0 || mexn12==0 || mexn56==0 ||
       mexsall==0 || mexs3478==0 || mexs34==0 || mexs78==0 ||
       mw1==0 || mw2==0 || temp==0 || mexpf1==0 || mexpf2==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__,__LINE__);
    exit(-1);
  }

  int *child = BOXES[ibox].child;
  double scale = SCALE[BOXES[ibox].level+1];
  int i, jj, offset, iexp1, iexp2;
  dcomplex zmul;

  // if exists, process nall list  
  if ( nnall > 0 ) {
    for ( i = 0; i < nnall; i++ ) {
      offset = inall[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixnall[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ixnall[i]-1];

	if ( ixnall[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ixnall[i]-1]);

	if ( iynall[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iynall[i]-1];

	if ( iynall[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iynall[i]-1]);

	mexnall[jj] += zmul*LEXPS[offset+jj];
      }
    }
  }

  // if exists, process n1256 list  
  if ( nn1256 > 0 ) { 
    for ( i = 0; i < nn1256; i++ ) {
      offset = in1256[i] *NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix1256[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix1256[i]-1];

	if ( ix1256[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix1256[i]-1]);

	if ( iy1256[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy1256[i]-1];

	if ( iy1256[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy1256[i]-1]);

	mexn1256[jj] += zmul*LEXPS[offset+jj];
      }
    }
  }

  // if exists, process n12 list  
  if ( nn12 > 0 ) {
    for ( i = 0; i < nn12; i++ ) {
      offset = in12[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix12[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix12[i]-1];

	if ( ix12[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix12[i]-1]);

	if ( iy12[i] > 0 )
	  zmul = zmul*YS[3*jj+iy12[i]-1];

	if ( iy12[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy12[i]-1]);

	mexn12[jj] += zmul*LEXPS[offset+jj];
      }
    }
  }

  // if exists, process n56 list  
  if ( nn56 > 0 ) {
    for ( i = 0; i < nn56; i++ ) {
      offset = in56[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix56[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix56[i]-1];

	if ( ix56[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix56[i]-1]);

	if ( iy56[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy56[i]-1];

	if ( iy56[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy56[i]-1]);

	mexn56[jj] += zmul*LEXPS[offset+jj];
      }
    }
  }

  // if exists, process sall list  
  if ( nsall > 0 ) {
    for ( i = 0; i < nsall; i++ ) {
      offset = isall[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixsall[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ixsall[i]-1]);

	if ( ixsall[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ixsall[i]-1];

	if ( iysall[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iysall[i]-1]);

	if ( iysall[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iysall[i]-1];

	mexsall[jj] += zmul*LEXPN[offset+jj];
      }
    }
  }

  // if exists, process s3478 list  
  if ( ns3478 > 0 ) {
    for ( i = 0; i < ns3478; i++ ) {
      offset = is3478[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix3478[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix3478[i]-1]);

	if ( ix3478[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix3478[i]-1];

	if ( iy3478[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy3478[i]-1]);

	if ( iy3478[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy3478[i]-1];

	mexs3478[jj] += zmul*LEXPN[offset+jj];
      }
    }
  }

  // if exists, process s34 list  
  if ( ns34 > 0 ) { 
    for ( i = 0; i < ns34; i++ ) {
      offset = is34[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix34[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix34[i]-1]);

	if ( ix34[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix34[i]-1];

	if ( iy34[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy34[i]-1]);

	if ( iy34[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy34[i]-1];

	mexs34[jj] += zmul*LEXPN[offset+jj];
      }
    }
  }

  // if exists, process s78 list  
  if ( ns78 > 0 ) {
    for ( i = 0; i < ns78; i++ ) {
      offset = is78[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix78[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix78[i]-1]);

	if ( ix78[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix78[i]-1];

	if ( iy78[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy78[i]-1]);

	if ( iy78[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy78[i]-1];

	mexs78[jj] += zmul*LEXPN[offset+jj];
      }
    }
  }

  // process child #1 
  if ( child[0] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    iexp1 = 0;
    iexp2 = 0;

    // receive nall, n1256, & n12 
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+2]*scale;
      iexp1 = iexp1+1;
    }

    if ( nn1256 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn1256[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( nn12 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn12[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive sall 
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++) 
	temp[jj] = mexsall[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;  
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #1 
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[child[0]*PGSZ], PGSZ);
    }
  }

  // process child #2
  if ( child[1] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive nall, n1256, & n12 
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+2]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn1256 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn1256[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn12 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn12[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive sall
    if ( nsall > 0) {
      for ( jj = 0; jj < NEXPTOTP; jj++) 
	temp[jj] = mexsall[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #2
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[1]], PGSZ);
    }
  }

  // process child #3 
  if ( child[2] > 0 ) {
    iexp1 = 0; 
    iexp2 = 0;

    // receive nall 
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1,temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive sall, s3478, & s34 
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+2]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns3478 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ )
	temp[jj] = temp[jj] + mexs3478[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns34 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ )
	temp[jj] = temp[jj] + mexs34[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) {
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #3
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[child[2]*PGSZ], PGSZ);
    }
  }

  // process child #4
  if ( child[3] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive nall 
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive sall, s3478, & s34
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+2]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns3478 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs3478[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns34 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs34[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 )
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #4
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[3]], PGSZ);
    }
  }

  // process child #5
  if ( child[4] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0; 
    iexp2 = 0;

    // receive nall, n1256, & n56
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+2]*conj(XS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn1256 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn1256[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn56 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn56[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive sall
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      TEL_PHASE1(mexpf2, temp, 0);
      iexp2 = iexp2+1;
    }

    // convert into local expansion of child #5
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[4]], PGSZ);
    }
  }

  // process child #6
  if ( child[5] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive nall, n1256, & n56
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+2]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn1256 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn1256[jj]*ZS[3*jj+1]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( nn56 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexn56[jj]*ZS[3*jj+1]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive sall
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+1]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #6
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[5]], PGSZ);
    }
  }

  // process child #7
  if ( child[6] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive nall
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1,temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive sall, s3478, & s78 
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+2]*XS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns3478 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs3478[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns78 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs78[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #7
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[6]], PGSZ);
    }
  } 

  // process child #8
  if ( child[7] > 0 ) {
    iexp1 = 0; 
    iexp2 = 0;

    // receive nall
    if ( nnall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexnall[jj]*ZS[3*jj+1]*conj(XS[3*jj]*YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive sall, s3478, & s78
    if ( nsall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexsall[jj]*ZS[3*jj+2]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns3478 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs3478[jj]*ZS[3*jj+1]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( ns78 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexs78[jj]*ZS[3*jj+1]*XS[3*jj]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #8
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1,  iexp2, mexpf2, iexp1, mexpf1,0);
      ROTY2Z(mw1, mw2);
      vecAdd(mw2, &LOCAL[PGSZ*child[7]], PGSZ);
    }
  }
  
  free (mexnall);
  free (mexn1256);
  free (mexn12);
  free (mexn56);
  free (mexsall);
  free (mexs3478);
  free (mexs34);
  free (mexs78);
  free (mw1);
  free (mw2);
  free (temp);
  free (mexpf1);
  free (mexpf2);

  return;
}

//  receive east and west direction exponentials 
/*
    after rotation, the layout of the boxes are the following 
       ^ z' 
       |             ---------
       |            / 8    4 /
       |           / 6    2 /
       |          ----------
       |  y'  
       |   /      ---------
       |  /      / 7    3 /
       | /      / 5    1 /
       |/      ----------            
----------------------------------> x' 
*/

void TEL_EW(const int ibox)
{
  int ieall[36], neall, ixeall[36], iyeall[36], ie1357[16], ne1357, ix1357[16], iy1357[16],
    ie13[4], ne13, ix13[4], iy13[4], ie57[4], ne57, ix57[4], iy57[4],
    ie1[4], ne1, ix1[4], iy1[4], ie3[4], ne3, ix3[4], iy3[4], 
    ie5[4], ne5, ix5[4], iy5[4], ie7[4], ne7, ix7[4], iy7[4], 
    iwall[36], nwall, ixwall[36], iywall[36], iw2468[16], nw2468, ix2468[16], iy2468[16],
    iw24[4], nw24, ix24[4], iy24[4], iw68[4], nw68, ix68[4], iy68[4],
    iw2[4], nw2, ix2[4], iy2[4], iw4[4], nw4, ix4[4], iy4[4], 
    iw6[4], nw6, ix6[4], iy6[4], iw8[4], nw8, ix8[4], iy8[4];
  int *child = BOXES[ibox].child;
  double scale = SCALE[BOXES[ibox].level+1];
  
  create_list_east (ibox, ieall, &neall, ixeall, iyeall, ie1357, &ne1357, ix1357, iy1357,
		    ie13, &ne13, ix13, iy13, ie57, &ne57, ix57, iy57,
		    ie1, &ne1, ix1, iy1, ie3, &ne3, ix3, iy3, 
		    ie5, &ne5, ix5, iy5, ie7, &ne7, ix7, iy7);
  create_list_west (ibox, iwall, &nwall, ixwall, iywall, iw2468, &nw2468, ix2468, iy2468,
		    iw24, &nw24, ix24, iy24, iw68, &nw68, ix68, iy68,
		    iw2, &nw2, ix2, iy2, iw4, &nw4, ix4, iy4, 
		    iw6, &nw6, ix6, iy6, iw8, &nw8, ix8, iy8);

  dcomplex *mexeall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe1357 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe13 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe57 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe1 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe3 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe5 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexe7 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexwall = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw2468 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw24 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw68 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw2 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw4= (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw6 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexw8 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf1 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mexpf2 = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  dcomplex *mw1 = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *mw2 = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *temp = (dcomplex *)calloc(NEXPMAX, sizeof(dcomplex));
  if ( mexeall==0 || mexe1357==0 || mexe13==0 || mexe57==0 ||
       mexe1==0 || mexe3==0 || mexe5==0 || mexe7==0 ||
       mexwall==0 || mexw2468==0 || mexw24==0 || mexw68==0 ||
       mexw2==0 || mexw4==0 || mexw6==0 || mexw8==0 ||
       mexpf1==0 || mexpf2==0 || mw1==0 || mw2==0 || temp==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  dcomplex zmul;
  int i, jj, offset, iexp1, iexp2;

  // if exists, process eall list  
  if ( neall > 0 ) {
    for ( i = 0; i < neall; i++ ) {
      offset = ieall[i] *NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixeall[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ixeall[i]-1];

	if ( ixeall[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ixeall[i]-1]);

	if ( iyeall[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iyeall[i]-1];

	if ( iyeall[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iyeall[i]-1]);

	mexeall[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  //if exists, process e1357 list  
  if ( ne1357 > 0 ) {
    for ( i = 0; i < ne1357; i++ ) {
      offset = ie1357[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix1357[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix1357[i]-1];

	if ( ix1357[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix1357[i]-1]);

	if ( iy1357[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy1357[i]-1];

	if ( iy1357[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy1357[i]-1]);

	mexe1357[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  // if exists, process e13 list  
  if ( ne13 > 0 ) {
    for ( i = 0; i < ne13; i++ ) {
      offset = ie13[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix13[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix13[i]-1];

	if ( ix13[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix13[i]-1]);

	if ( iy13[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy13[i]-1];

	if ( iy13[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy13[i]-1]);

	mexe13[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  //if exists, process e57 list  
  if ( ne57 > 0 ) {   
    for ( i = 0; i < ne57; i++ ) {
      offset = ie57[i] *NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix57[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix57[i]-1];

	if ( ix57[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix57[i]-1]);

	if ( iy57[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy57[i]-1];

	if ( iy57[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy57[i]-1]);

	mexe57[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  // if exists, process e1 list  
  if ( ne1 > 0 ) {
    for ( i = 0; i < ne1; i++ ) {
      offset = ie1[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix1[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix1[i]-1];

	if ( ix1[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix1[i]-1]);

	if ( iy1[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy1[i]-1];

	if ( iy1[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy1[i]-1]);

	mexe1[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  // if exists, process e3 list  
  if ( ne3 > 0 ) {
    for ( i = 0; i < ne3; i++ ) {
      offset = ie3[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix3[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix3[i]-1];

	if ( ix3[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix3[i]-1]);

	if ( iy3[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy3[i]-1];

	if ( iy3[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy3[i]-1]);
	
	mexe3[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }
  
  // if exists, process e5 list  
  if ( ne5 > 0 ) {
    for ( i = 0; i < ne5; i++ ) {
      offset = ie5[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix5[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix5[i]-1];

	if ( ix5[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix5[i]-1]);

	if ( iy5[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy5[i]-1];

	if ( iy5[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy5[i]-1]);

	mexe5[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  // if exists, process e7 list  
  if ( ne7 > 0 ) {
    for ( i = 0; i < ne7; i++ ) {
      offset = ie7[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix7[i] > 0 ) 
	  zmul = zmul*XS[3*jj+ix7[i]-1];

	if ( ix7[i] < 0 ) 
	  zmul = zmul*conj(XS[3*jj-ix7[i]-1]);

	if ( iy7[i] > 0 ) 
	  zmul = zmul*YS[3*jj+iy7[i]-1];

	if ( iy7[i] < 0 ) 
	  zmul = zmul*conj(YS[3*jj-iy7[i]-1]);

	mexe7[jj] += zmul*LEXPW[offset+jj];
      }
    }
  }

  // if exists, process wall list  
  if ( nwall > 0 ) {
    for ( i = 0; i < nwall; i++ ) {
      offset = iwall[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ixwall[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ixwall[i]-1]);
	
	if ( ixwall[i] < 0 )
	  zmul = zmul*XS[3*jj-ixwall[i]-1];
	
	if ( iywall[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iywall[i]-1]);
	
	if ( iywall[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iywall[i]-1];
	
	mexwall[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // if exists, process w2468 list  
  if ( nw2468 > 0 ) {
    for ( i = 0; i < nw2468; i++ ) {
      offset = iw2468[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix2468[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix2468[i]-1]);
	
	if ( ix2468[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix2468[i]-1];
	
	if ( iy2468[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy2468[i]-1]);
	
	if ( iy2468[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy2468[i]-1];
	
	mexw2468[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }
  
  // if exists, process w24 list  
  if ( nw24 > 0 ) {
    for ( i = 0; i < nw24; i++ ) {
      offset = iw24[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix24[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix24[i]-1]);
	
	if ( ix24[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix24[i]-1];

	if ( iy24[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy24[i]-1]);
	
	if ( iy24[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy24[i]-1];

	mexw24[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // if exists, process w68 list  
  if ( nw68 > 0 ) {
    for ( i = 0; i < nw68; i++ ) {
      offset = iw68[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix68[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix68[i]-1]);

	if ( ix68[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix68[i]-1];

	if ( iy68[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy68[i]-1]);

	if ( iy68[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy68[i]-1];

	mexw68[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  //if exists, process w2 list  
  if ( nw2 > 0 ) {
    for ( i = 0; i < nw2; i++ ) {
      offset = iw2[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix2[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix2[i]-1]);

	if ( ix2[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix2[i]-1];

	if ( iy2[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy2[i]-1]);

	if ( iy2[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy2[i]-1];

	mexw2[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // if exists, process w4 list  
  if ( nw4 > 0 ) {
    for ( i = 0; i < nw4; i++ ) {
      offset = iw4[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix4[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix4[i]-1]);

	if ( ix4[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix4[i]-1];

	if ( iy4[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy4[i]-1]);

	if ( iy4[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy4[i]-1];

	mexw4[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // if exists, process w6 list  
  if ( nw6 > 0 ) {
    for ( i = 0; i < nw6; i++ ) {
      offset = iw6[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix6[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix6[i]-1]);

	if ( ix6[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix6[i]-1];

	if ( iy6[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy6[i]-1]);

	if ( iy6[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy6[i]-1];

	mexw6[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // if exists, process w8 list  
  if ( nw8 > 0 ) {
    for ( i = 0; i < nw8; i++ ) {
      offset = iw8[i]*NEXPTOTP;
      for ( jj = 0; jj < NEXPTOTP; jj++ ) {
	zmul = 1;
	if ( ix8[i] > 0 ) 
	  zmul = zmul*conj(XS[3*jj+ix8[i]-1]);

	if ( ix8[i] < 0 ) 
	  zmul = zmul*XS[3*jj-ix8[i]-1];

	if ( iy8[i] > 0 ) 
	  zmul = zmul*conj(YS[3*jj+iy8[i]-1]);

	if ( iy8[i] < 0 ) 
	  zmul = zmul*YS[3*jj-iy8[i]-1];

	mexw8[jj] += zmul*LEXPE[offset+jj];
      }
    }
  }

  // process child #1
  if ( child[0] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive eall, e1357, e13, e1
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+2]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne1357 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe1357[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne13 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe13[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne1 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe1[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive wall
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #1
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[0]], PGSZ);
    }
  }

  // process child #2
  if ( child[1] > 0 ) {  
    iexp1 = 0;
    iexp2 = 0;

    // receive eall
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+1]*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive wall, w2468, w24, & w2
    if ( nwall > 0 ) {      
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+2]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw2468 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw2468[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw24 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw24[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw2 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw2[jj]*ZS[3*jj+1]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #2
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[1]], PGSZ);
    }
  }

  // process child #3
  if ( child[2] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive eall, e1357, e13, & e3
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+2]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne1357 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe1357[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne13 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe13[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne3 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe3[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1,temp,0);

    // receive wall
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #3
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[2]], PGSZ);
    }
  }

  // process child #4 
  if ( child[3] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive eall
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+1]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    // receive wall, w2468, w24, & w4
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+2]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw2468 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw2468[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw24 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw24[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw4 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw4[jj]*ZS[3*jj+1]*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #4
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[3]], PGSZ);
    }
  }

  // process child #5 
  if ( child[4] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive eall, e1357, e57, & e5
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+2]*XS[3*jj]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne1357 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe1357[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne57 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe57[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp1 = iexp1+1;
    }

    if ( ne5 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe5[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp1 = iexp1 + 1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive wall
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #5
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[4]], PGSZ);
    }
  }

  // process child #6 
  if ( child[5] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive eall
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+1]*XS[3*jj]*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive wall, w2468, w68, & w6
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+2]*conj(XS[3*jj])*scale;
      iexp2 = iexp2+1;
    }

    if ( nw2468 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw2468[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp2 = iexp2+1;
    }

    if ( nw68 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw68[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp2 = iexp2+1;
    }

    if ( nw6 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw6[jj]*ZS[3*jj+1]*conj(XS[3*jj])*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #6
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[5]], PGSZ);
    }
  }

  // process child #7
  if ( child[6] > 0 ) {
    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;
    iexp1 = 0;
    iexp2 = 0;

    // receive eall, e1357, e57, & e7
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+2]*XS[3*jj]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne1357 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe1357[jj]*ZS[3*jj+1]*XS[3*jj]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne57 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe57[jj]*ZS[3*jj+1]* XS[3*jj]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( ne7 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexe7[jj]*ZS[3*jj+1]*XS[3*jj]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
    }

    if ( iexp1 > 0 ) 
      TEL_PHASE1(mexpf1, temp,0);

    // receive wall
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+1]*conj(XS[3*jj])*YS[3*jj]*scale;
      iexp2 = iexp2+1;
      TEL_PHASE1(mexpf2, temp, 0);
    }

    // convert into local expansion of child #7
    if ( iexp1+iexp2) {
      TEL_PHASE2( mw1, iexp2, mexpf2, iexp1, mexpf1,0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[6]], PGSZ);
    }
  }

  // process child #8
  if ( child[7] > 0 ) {
    iexp1 = 0;
    iexp2 = 0;

    // receive eall
    if ( neall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexeall[jj]*ZS[3*jj+1]*XS[3*jj]*conj(YS[3*jj])*scale;
      iexp1 = iexp1+1;
      TEL_PHASE1(mexpf1, temp,0);
    }

    for ( jj = 0; jj < NEXPTOTP; jj++ ) 
      temp[jj] = 0;

    // receive wall, w2468, w68, & w8
    if ( nwall > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = mexwall[jj]*ZS[3*jj+2]*conj(XS[3*jj])*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw2468 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw2468[jj]*ZS[3*jj+1]*conj(XS[3*jj])*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw68 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw68[jj]*ZS[3*jj+1]*conj(XS[3*jj])*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( nw8 > 0 ) {
      for ( jj = 0; jj < NEXPTOTP; jj++ ) 
	temp[jj] = temp[jj] + mexw8[jj]*ZS[3*jj+1]*conj(XS[3*jj])*YS[3*jj]*scale;
      iexp2 = iexp2+1;
    }

    if ( iexp2 > 0 ) 
      TEL_PHASE1(mexpf2, temp, 0);

    // convert into local expansion of child #8
    if ( iexp1+iexp2) {
      TEL_PHASE2(mw1, iexp2, mexpf2, iexp1, mexpf1, 0);
      ROTZ2X(mw1, mw2, 0);
      vecAdd(mw2, &LOCAL[PGSZ*child[7]], PGSZ);
    }
  }

  free (mexeall);
  free (mexe1357);
  free (mexe13);
  free (mexe57);
  free (mexe1);
  free (mexe3);
  free (mexe5);
  free (mexe7);
  free (mexwall);
  free (mexw2468);
  free (mexw24);
  free (mexw68);
  free (mexw2);
  free (mexw4);
  free (mexw6);
  free (mexw8);
  free (mexpf1);
  free (mexpf2);
  free (temp);
  free (mw1);
  free (mw2);

  return;
}

void TEL_PHASE1(dcomplex *mexpf, const dcomplex *mexpphys, const int level)
{
  int i;
  int nftot = 0;
  int nptot = 0;
  int next = 0;

  for ( i = 0; i < NLAMBS; i++ ) {
    int nalpha = NUMPHYS[i];
    int nalpha2 = nalpha/2;
    mexpf[nftot] = 0;
    int ival;
    for ( ival = 0; ival < nalpha2; ival++) {
      mexpf[nftot] += 2.0*creal( mexpphys[nptot+ival] );
    }
    mexpf[nftot] /= nalpha;

    int nm;
    for ( nm = 2; nm < NUMFOUR[i]; nm=nm+2) {
      mexpf[nftot+nm] = 0;
      for ( ival = 0; ival < nalpha2; ival++) {
	double rtmp = 2*creal( mexpphys[nptot+ival] );
	mexpf[nftot+nm] += FEXPBACK[next]*rtmp;
	next += 1;
      }
      mexpf[nftot+nm] /= nalpha;
    }

    for ( nm = 1; nm < NUMFOUR[i]; nm=nm+2 ) {
      mexpf[nftot+nm] = 0;
      for ( ival = 0; ival < nalpha2; ival++) {
	dcomplex ztmp = 2*cimag(mexpphys[nptot+ival])*_Complex_I;
	mexpf[nftot+nm] += FEXPBACK[next]*ztmp;
	next += 1;
      }
      mexpf[nftot+nm] /= nalpha;
    }
    nftot += NUMFOUR[i];
    nptot += NUMPHYS[i]/2;
  }
}

void TEL_PHASE2(dcomplex *local, const int iexpu, const dcomplex *mexpup, 
		const int iexpd, const dcomplex *mexpdown, const int level)
{
  double *rlampow = (double *)calloc(1+PTERMS, sizeof(double));
  dcomplex *zeye = (dcomplex *)calloc(1+PTERMS, sizeof(dcomplex));
  dcomplex *mexpplus = (dcomplex *)calloc(NEXPTOT, sizeof(dcomplex));
  dcomplex *mexpminus = (dcomplex *)calloc(NEXPTOT, sizeof(dcomplex));
  if ( rlampow==0 || zeye==0 || mexpplus==0 || mexpminus==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }
  
  int i;
  zeye[0] = 1.0;
  for ( i = 1; i <= PTERMS; i++ ) {
    zeye[i] = zeye[i-1]*_Complex_I;
  }
  for ( i = 0; i < PGSZ; i++ ) 
    local[i] = 0.0;

  // compute sum and difference of mexpup and mexpdown
  int nm;
  for ( nm = 0; nm < NEXPTOT; nm++ ) {
    if ( iexpu <= 0 ) {
      mexpplus[nm] = mexpdown[nm];
      mexpminus[nm] = mexpdown[nm];
    } else if ( iexpd <= 0 ) {
      mexpplus[nm] = mexpup[nm];
      mexpminus[nm] = -mexpup[nm];
    } else {
      mexpplus[nm] = mexpdown[nm] + mexpup[nm];
      mexpminus[nm] = mexpdown[nm] - mexpup[nm];
    }
  }

  // loop over multipole order to generate mexp values
  int ntot = 0;
  int nell;
  for ( nell = 0; nell < NLAMBS; nell++ ) {
    rlampow[0] = WHTS[nell];
    double rmul = RLAMS[nell];
    int j;
    for ( j = 1; j <= PTERMS; j++ ) {
      rlampow[j] = rlampow[j-1]*rmul;
    }

    // add contribution to local expansion
    int mmax = NUMFOUR[nell]-1;
    int mth;
    int ncurrent;

    for (mth = 0; mth<=mmax; mth=mth+2) {
      int offset = mth*(PTERMS+1);
      for ( nm = mth; nm <= PTERMS; nm = nm+2 ) {
	int index = offset+nm;
	ncurrent = ntot+mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpplus[ncurrent];
      }

      for ( nm = mth+1; nm <= PTERMS; nm=nm+2) {
	int index = offset+nm;
	ncurrent = ntot+mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpminus[ncurrent];
      }
    }

    for ( mth = 1; mth<=mmax; mth=mth+2) {
      int offset = mth*(PTERMS+1);
      for ( nm = mth+1; nm<=PTERMS; nm=nm+2) {
	int index = nm+offset;
	ncurrent = ntot+mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpplus[ncurrent];
      }

      for ( nm = mth; nm <= PTERMS; nm=nm+2) {
	int index = nm+offset;
	ncurrent = ntot+mth;
	rmul = rlampow[nm];
	local[index] += rmul*mexpminus[ncurrent];
      }
    }
    ntot = ntot + NUMFOUR[nell];
  }

  // scale the expansions according to formula
  int mth;
  for ( mth = 0; mth <= PTERMS; mth++ ) {
    int offset = mth*(PTERMS+1);
    for ( nm = mth; nm <= PTERMS; nm++ ) {
      int index = nm+offset;
      local[index] = local[index]*zeye[mth]*YTOPCS[nm+mth*61];
    }
  }

  free (rlampow);
  free (zeye);
  free (mexpplus);
  free (mexpminus);
}


void TLL(const int pbox)
{
  // TLL translate the local expansion of a nonleaf box to its children.

  static dcomplex var[5] = {1, 1-_Complex_I, -1-_Complex_I, -1+_Complex_I, 1+_Complex_I};
  const double arg = sqrt(2)/2.0;
  dcomplex *local = &LOCAL[PGSZ*pbox];

  dcomplex *localn = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *marray = (dcomplex *)calloc(PGSZ, sizeof(dcomplex));
  dcomplex *ephi = (dcomplex *)calloc(1+PTERMS, sizeof(dcomplex));
  double *powers = (double *)calloc(1+PTERMS, sizeof(double));
  if ( localn==0 || marray==0 || ephi==0 || powers==0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  int i;
  for ( i = 0; i < 8; i++ ) {
    int cbox = BOXES[pbox].child[i];
    int ell, m, mp, j, k;
    if ( cbox > 0 ) {     
      int ifl = IFL_DN[i]; 
      double *rd = ( i < 4 ? RDSQ3 : RDMSQ3);
      int level = BOXES[cbox].level;
      double sc1 = SCALE[level-1];
      double sc2 = SCALE[level];

      double dd;

      ephi[0] = 1.0;
      ephi[1] = arg*var[ifl];

      dd = -sqrt(3)/4.0;
      powers[0] = 1.0;
      for ( ell = 1; ell <= PTERMS; ell++ ) {
	powers[ell] = powers[ell-1]*dd;
      }
      for ( ell = 2; ell <= PTERMS; ell++ ) {
	ephi[ell] = ephi[ell-1]*ephi[1];
      }

      // rotation of phi radians about the z-axis in the original coordinate system
      for ( m = 0; m <= PTERMS; m++ ) {
	int offset = m*(PTERMS+1);
	for ( ell = m; ell <= PTERMS; ell++ ) {
	  int index = ell+offset;
	  localn[index] = conj( ephi[m] )*local[index];
	}
      }
      
      // a rotation about the y'-axis to align z' axis
      for ( m = 0; m <= PTERMS; m++ ) {
	int offset = m*(PTERMS+1);
	int offset1 = (m+PTERMS)*PGSZ;
	int offset2 = (-m+PTERMS)*PGSZ;
	for ( ell = m; ell <= PTERMS; ell++ ) {
	  int index = ell+offset;
	  marray[index] = localn[ell]*rd[ell+offset1];
	  for ( mp = 1; mp <= ell; mp++) {
	    int index1 = ell+mp*(PTERMS+1);
	    marray[index] = marray[index] + localn[index1]*rd[index1+offset1]+
	      conj( localn[index1] )*rd[index1+offset2];
	  }
	}
      }
 
      // shift along z' -axis
      for (  k = 0; k <= PTERMS; k++) {
	int offset = k*(PTERMS+1);
	for ( j = k; j<=PTERMS; j++) {
	  int index = j+offset;
	  localn[index] = marray[index];
	  for ( ell = 1; ell <= PTERMS-j; ell++) {
	    int index1 = ell+index;
	    int index2 = ell+j+k+ell*(2*PTERMS+1);
	    int index3 = ell+j-k+ell*(2*PTERMS+1);
	    localn[index] = localn[index] + marray[index1]*powers[ell]*DC[index2]*DC[index3];
	  }
	}
      }

      // rotate back about the y'-axis
      for ( m= 0; m <= PTERMS; m++) {
	int offset = m*(PTERMS+1);
	int offset1 = (m+PTERMS)*PGSZ;
	int offset2 = (-m+PTERMS)*PGSZ;
	for ( ell = m; ell <= PTERMS; ell++) {
	  int index = ell+offset;
	  marray[index] = localn[ell]*rd[ell+offset1];
	  for ( mp = 1; mp <= ell; mp=mp+2) {
	    int index1 = ell+mp*(PTERMS+1);
	    marray[index] = marray[index] - localn[index1]*rd[index1+offset1] - 
	      conj( localn[index1] )*rd[index1+offset2];
	  }
	  for ( mp = 2; mp <= ell; mp=mp+2) {
	    int index1 = ell+mp*(PTERMS+1);
	    marray[index] = marray[index] + localn[index1]*rd[index1+offset1] + 
	      conj( localn[index1] )*rd[index1 + offset2];
	  }
	}
      }

      for ( m = 1; m <= PTERMS; m=m+2 ) {
	int offset = m*(PTERMS+1);
	int offset1 = (m+PTERMS)*PGSZ;
	int offset2 = (-m+PTERMS)*PGSZ;
	for ( ell = m; ell <= PTERMS; ell++) {
	  int index = ell+offset;
	  marray[index] = -localn[ell]*rd[ell+offset1];
	  for ( mp = 1; mp <= ell; mp = mp+2) {
	    int index1 = ell+mp*(PTERMS+1);
	    marray[index] = marray[index] + localn[index1]*rd[index1+offset1]+
	      conj( localn[index1] )*rd[index1+offset2];
	  }
	  for ( mp = 2; mp <= ell; mp = mp+2) {
	    int index1 = ell+mp*(PTERMS+1);
	    marray[index] = marray[index] - localn[index1]*rd[index1+offset1] - 
	      conj( localn[index1] )*rd[index1+offset2];
	  }
	}
      }

      // rotate back about the z-axis and scale from sc1 to sc2
      powers[0] = 1.0;
      dd = sc1/sc2;
      for ( ell = 1; ell <= PTERMS; ell++) {
	powers[ell] = powers[ell-1]*dd;
      }
      
      for ( m= 0; m <= PTERMS; m++) {
	int offset = m*(PTERMS+1);
	for ( ell = m; ell <= PTERMS; ell++ ) {
	  int index = offset+ell;
	  localn[index] = ephi[m]*marray[index]*powers[ell];
	}
      }

      dcomplex *cLocal = &LOCAL[PGSZ*cbox];
      for ( m = 0; m < PGSZ; m++ ) 
	cLocal[m] += localn[m];
    }
  }

  free (ephi);
  free (powers);
  free (localn);
  free (marray);
}


