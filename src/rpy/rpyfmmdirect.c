#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "fmm_ds.h"
#include "fmm_utils.h"
#include "fmm_global.h"
#include "fmm_laplace.h"
#include "rpy.h"

/***********************************************************************
   Direct calculation of various (free space) RPY Green's function in R^3
          For D*v
************************************************************************/
void rpy3nearfield()
{
     int ibox;
     /*printf("BEFORE RPYPOT\n");
     int i;
      for(i=0;i<NPARTS; i++){
	     printf("%f %f %f \n", RPYPOT[3*i], RPYPOT[3*i+1], RPYPOT[3*i+2]);
      } */

     for( ibox =CONTENT[2]; ibox<=NBOXES; ibox++ ){
	 rpyLIST3(ibox); 
	 rpyLIST1(ibox);
     }

}
/***************************************************************************
 rpyLIST3 AND rpyLIST1
***************************************************************************/
void rpyLIST3(const int ibox)
{   
    if ( LIST[ibox].list3 != 0 ) {
		int nlist = LIST[ibox].list3[0];
		int j;
	
		for ( j = 1; j <= nlist; j++ ) {
			int member = LIST[ibox].list3[j];
			rpy_DIRECT_D(ibox, member, RPYPOT);
		}
    }
}

void rpyLIST1(const int ibox)
{
	if ( LIST[ibox].list1 != 0 ) {
	     int nlist = LIST[ibox].list1[0];
	     int j;
	     
	    for( j = 1; j <= nlist; j++ ) {
			int member = LIST[ibox].list1[j];
			if ( ibox == member ) {
				rpy_DIRECT_S(ibox, RPYPOT);
			} else {
				rpy_DIRECT_D(ibox, member, RPYPOT);
			}
	    }
	}
}

// Direct evaluation of rpy kernel for different iboxes
void rpy_DIRECT_D(const int ibox1, const int ibox2, double *rpypot)
{ 
	int start1 =BOXES[ibox1].addr - 1;
	int num1 = BOXES[ibox1].npts;
  	int end1 = start1+num1-1;

	int start2 = BOXES[ibox2].addr-1;
    	int num2 = BOXES[ibox2].npts;
        int end2 = start2+num2-1;
      	
	int i, j,i3;

	for ( i = start1; i <= end1; i++ ) {
	    double newpot[3], pot1[3];
	    i3 = i*3;
	    pot1[0]=0.0; pot1[1]=0.0; pot1[2] =0.0; 
	    for ( j = start2; j <= end2; j++ ) {
			rpy3sup_eval(&FMMLOC[i3], &FMMLOC[3*j], &RPYCHARGE[3*j], newpot);
			pot1[0] += newpot[0];	  
			pot1[1] += newpot[1];
			pot1[2] += newpot[2];
	    }

	    rpypot[i3] += pot1[0];
	    rpypot[i3+1] += pot1[1];
	    rpypot[i3+2] += pot1[2];
        
	}
}

/***********************************************************************
***********************************************************************/
void rpy_DIRECT_S(const int ibox, double *rpypot)
{
// rpy_DIRECT_S computes the interaction between particles within the same box.
    int start = BOXES[ibox].addr-1;
    int num = BOXES[ibox].npts;
    int end = start+num-1;
    int i, j, i3;

    for ( i = start; i <= end; i++ ) {
	    i3 = i*3;
	    double newpot[3], pot1[3];
	    pot1[0] = 0.0; pot1[1] = 0.0; pot1[2]=0.0;
	 //   for ( j = start; j <= i-1; j++ ) { Aug 25 Wenhua
	    for ( j = start; j <= end; j++ ) {
		  rpy3sup_eval(&FMMLOC[i3], &FMMLOC[3*j], &RPYCHARGE[3*j], newpot);
		  pot1[0] += newpot[0];
		  pot1[1] += newpot[1];
		  pot1[2] += newpot[2];

	    }

/*	    for ( j = i+1; j <= end; j++ ) {
		  rpy3sup_eval(&FMMLOC[i3], &FMMLOC[3*j], &RPYCHARGE[3*j], newpot);
		  pot1[0] += newpot[0];
		  pot1[1] += newpot[1];
		  pot1[2] += newpot[2];

	    }*/
			
	    rpypot[i3] += pot1[0];
	    rpypot[i3+1] += pot1[1];
	    rpypot[i3+2] += pot1[2];
    }
}


/**********************************************************************
  rpy3sup_eval--pairwise interaction for rpy kernel
** ********************************************************************/

/*
RPY function

This function evaluates  the potential vector -pot at the location -target due
to the static force -charge at the source.

Input: T, S, cearge, ha
Output: pot
*/
void rpy3sup_eval(const double *T, const double *S, const double *charge, double pot[3]) 
{
  /*
  printf("Interaction:\n");
  printf("T: %f %f %f\n", T[0], T[1], T[2]);
  printf("S: %f %f %f\n", S[0], S[1], S[2]);*/
  double rx = T[0] - S[0], ry = T[1] - S[1], rz = T[2] - S[2];
  double rr = rx*rx + ry*ry + rz*rz;
  double rdis = sqrt(rr);
  double rdisinv = 1.0/rdis;
  double rdisinv3 = rdisinv/rr;
  double mindif = 1.0e-12;

  double dd = charge[0] * rx + charge[1] * ry + charge[2] * rz;
  dd = dd * rdisinv;
 
  double ha = 1.0;
  double done = 1.0;
  double pi = 4.0*atan(done);
  double cmu =1.0;
  double c0 = cmu/(6.0*pi*ha);
  double c1 = cmu/(8.0*pi);
  double c2 = cmu*ha*ha/(12.0*pi);

  // Rotne Prager tensor : direct evaluation 
  if(rdis >= 2*ha) {
    dd = dd/rr;
    pot[0] = charge[0]*rdisinv + dd *rx;
    pot[1] = charge[1]*rdisinv + dd *ry;
    pot[2] = charge[2]*rdisinv + dd *rz;

    double c2rinv3 = c2 * rdisinv3;
    double c2dd = 3.0 * c2 * dd /rr;
   
    pot[0] = c1*pot[0];
    pot[0] += c2rinv3 * charge[0] - c2dd *rx;
    pot[1] = c1*pot[1];
    pot[1] += c2rinv3 * charge[1] - c2dd *ry;
    pot[2] = c1*pot[2];
    pot[2] +=  c2rinv3 * charge[2] - c2dd *rz;

  } else if((rdis<2*ha) &&( rdis> mindif)) {

	  double c3 = 9.0/(32.0*ha), c4 = c3/3.0;
	  double ans1 = c0*(1.0 - c3 *rdis);
	  double ans2 = c4*dd*c0;

	  pot[0] = ans1 * charge[0] + ans2 * rx;
	  pot[1] = ans1 * charge[1] + ans2 * ry;
	  pot[2] = ans1 * charge[2] + ans2 * rz;
         
  }else if (rdis <= mindif){
	
	pot[0] = c0 * charge[0];
	pot[1] = c0 * charge[1];
	pot[2] = c0 * charge[2];
  } 
}
