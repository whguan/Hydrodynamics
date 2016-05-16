/*
  traverse_sym_graph.c: traverse adaptive fmm graph in parallel.
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

void adap_fmm_compute(int ifdirect, int ifpot, int ifgrad, int ifhess)
{
  // The compute function traverse the adaptive fmm graph in parallel. 
  // (1) Near-field and far-field interactions are processed in parallel. 
  // (2) Near-field includes the process of lists 1 and 3. 
  // (3) Far-field includes (a) an upward traversal of the source tree, 
  // including TSM, TMM, and TME operations; and (b) an downward traversal
  // of the target tree, including TEL, TSL, TLL, and TLT operations.

  // Note: when near-field and far-field are scheduled concurrently, it 
  // increases the memory consumption. Separate copies are allocated to hold
  // the potential and field results. To reduce synchronization between 
  // far-field and near-field computation, list 3 is processed using direct
  // pairwise interaction. (Evaluating multipole expansion can be more 
  // efficient in certain case, but requires more synchronization). 
// wenhua
  //cilk_spawn ProcFarField();
 // cilk_spawn ProcNearField();
 	// march 3th changed by wenhua
	// Note that adap_fmm_compute just calculate either far field or near filed
	//
	if (ifdirect ==1) {
	//	printf("begin near field for rpy\n");
		rpy3nearfield();
	//	printf("end near filed for rpy\n");
	} else {
		ProcFarField(ifpot, ifgrad, ifhess);
	}
}

void ProcFarField( int ifpot, int ifgrad, int ifhess)
{
  Aggregate(1);
  Disaggregate(1, ifpot, ifgrad,ifhess);
}

void ProcNearField(void)
{

  int ibox;
  for(ibox =CONTENT[2];ibox<=NBOXES; ibox++){
	LIST3(ibox); // list1 and list3 are changed to rpy_direct instead of direct evaluation for laplace kernel
  	LIST1(ibox);
  }
}

void Aggregate(const int ibox) 
{
  
  if ( BOXES[ibox].nchild == 0 ) {
  	TSM(ibox); // ibox is a leaf box, generate its multipole expansion
  } else {
    	int i;
   	 // wenhua 2015
   	//cilk_for (i = 0; i < 8; i++ ) {
    	for(i = 0; i<8;i++){
      		int child = BOXES[ibox].child[i];
      		if ( child ) {
			Aggregate(child);
      		}
    	}
    	TMM(ibox); // ibox is a nonleaf box, translate its children's mutlipole expansions.
  }

  // Note: this implementation assumes free space boundary condition, which means 
  // no box at levels 0 and 1 are well-separated and their multipole expansions do not
  // need to be translated into exponential expansions. 
  if ( ibox >= CONTENT[4] ) 
    TME(ibox); // convert ibox's multipole expansion into exponential expansion.
}


void Disaggregate (int ibox,  int ifpot, int ifgrad, int ifhess)
{
  LIST4(ibox);
  if ( BOXES[ibox].nchild == 0 ) {
    TLT(ibox,ifpot, ifgrad, ifhess); // ibox is a leaf box, evaluate its local expansion
  } else {
    TEL(ibox); // ibox is a nonleaef box, completes exponential-to-local using merge-and-shift.
    TLL(ibox); // ibox is a nonleaf box, translate its local expansion to its children
    int i;
    //wenhua Aug 7
   // cilk_for ( i = 0; i < 8; i++ ) {
    for(i = 0; i< 8; i++){
      int child = BOXES[ibox].child[i];
      if ( child ) {
	Disaggregate(child, ifpot, ifgrad, ifhess);
      }
    }
  }
}

void LIST4(const int ibox)
{
  if ( LIST[ibox].list4 != 0 ) {
    int nlist = LIST[ibox].list4[0];
    int choice =  BOXES[ibox].npts > PGSZ;
    int j;
    if ( choice ) {
      // ibox has too many points, updating its local expansion is more efficient
      for ( j = 1; j <= nlist; j++ )
	    TSL(ibox, LIST[ibox].list4[j]);
    } else {
      // ibox has very few points, updating potential/field result directly is more efficient
      for ( j = 1; j <= nlist; j++ )
	DIRECT_D(ibox, LIST[ibox].list4[j], FMMPOT, FMMFIELD,FMMDXX, FMMDYY, FMMDZZ, FMMDXY, FMMDXZ, FMMDYZ);
    }
  }
}

void LIST3(const int ibox)
{
  if ( LIST[ibox].list3 != 0 ) {
    int nlist = LIST[ibox].list3[0];
    int j;
    for ( j = 1; j <= nlist; j++ ) {
      int member = LIST[ibox].list3[j];
      DIRECT_D(ibox, member, FMMPOTN, FMMFIELDN, FMMDXXN, FMMDYYN, FMMDZZN, FMMDXYN, FMMDXZN, FMMDYZN);
      
      /* If near-field and far-field are not scheduled concurrently, 
         the following code is more efficient way to process list 3. 

      int choice = BOXES[member].npts > PGSZ;
      if ( choice ) {
	// the member has too many points, evaluting its multipole is more efficient
	int k;
	for ( k = 0; k < BOXES[ibox].npts; k++ )
	  TMT(ibox, k, member); // evlaute member's multipole at kth particle location in ibox
      } else {
	// the member has very few points, pairwise interaction with ibox is more efficient
	DIRECT(ibox, member, 'n');
      }
      */
    }
  }
}

void LIST1(const int ibox)
{
  if ( LIST[ibox].list1 != 0 ) {
    int nlist = LIST[ibox].list1[0];
    int j;
    for ( j = 1; j <= nlist; j++ ) {
      int member = LIST[ibox].list1[j];
      if ( ibox == member ) {
	DIRECT_S(ibox, FMMPOTN, FMMFIELDN, FMMDXXN, FMMDYYN, FMMDZZN, FMMDXYN, FMMDXZN, FMMDYZN);
      } else {
	DIRECT_D(ibox, member, FMMPOTN, FMMFIELDN, FMMDXXN, FMMDYYN, FMMDZZN, FMMDXYN, FMMDXZN, FMMDYZN);
      }
    }
  }
}
  
inline void DIRECT_D( const int ibox1, const int ibox2, double *fmmpot, double *fmmfield,
		double *fmmdxx, double *fmmdyy, double *fmmdzz, double *fmmdxy, double *fmmdxz, double *fmmdyz)
{
  // DIRECT_D computes the interaction between a target box (ibox1) and a source
  // box (ibox2). ibox1 and ibox2 are different boxes. 

  int start1 = BOXES[ibox1].addr-1;
  int num1 = BOXES[ibox1].npts;
  int start2 = BOXES[ibox2].addr-1;
  int num2 = BOXES[ibox2].npts;

  int i, j, end1, end2, i3;
  double newpot, newf1, newf2, newf3;
  double newdxx, newdyy, newdzz, newdxy, newdxz, newdyz;
  double pot1, field1, field2, field3;
  double dxx, dyy, dzz, dxy, dxz, dyz;

  end1 = start1+num1-1;
  end2 = start2+num2-1;

  for ( i = start1; i <= end1; i++ ) {
    pot1 = 0;
    field1 = 0; field2 = 0; field3 = 0;
    dxx = 0.0; 	dxy = 0.0;  dxz = 0.0;
    dyy = 0.0;  dyz = 0.0;  dzz = 0.0;
	
    i3 = i*3;
    
    for ( j = start2; j <= end2; j++ ) {
      kernel(&FMMLOC[i3], &FMMLOC[3*j], FMMCHARGE[j], &newpot, &newf1, &newf2, &newf3, 
      	         &newdxx, &newdyy, &newdzz, &newdxy, &newdxz, &newdyz);
      //rpy3sup_eval(&FMMLOC[i3], &FMMLOC[3*j], FMMCHARGE[j], &newpot);
      pot1 += newpot;
      field1 += newf1;
      field2 += newf2;
      field3 += newf3;
      dxx += newdxx;
      dyy += newdyy;
      dzz += newdzz;
      dxy += newdxy;
      dxz += newdxz;
      dyz += newdyz;	  
    }

    fmmpot[i] += pot1;
    fmmfield[i3] += field1;
    fmmfield[i3+1] += field2;
    fmmfield[i3+2] += field3;
    
    fmmdxx[i] += dxx;
    fmmdyy[i] += dyy;
    fmmdzz[i] += dzz;
    fmmdxy[i] += dxy;
    fmmdxz[i] += dxz;
    fmmdyz[i] += dyz;
  }
}
inline void DIRECT_S(const int ibox, double *fmmpot, double *fmmfield,
		double *fmmdxx, double *fmmdyy, double *fmmdzz, 
		double *fmmdxy, double *fmmdxz, double *fmmdyz)
{
  // DIRECT_S computes the interaction between particles within the same box.
  int start = BOXES[ibox].addr-1;
  int num = BOXES[ibox].npts;
  int end = start+num-1;
  int i, j, i3;
  double newpot, newf1, newf2, newf3;
  double pot1, field1, field2, field3;
  double newdxx, newdyy, newdzz, newdxy, newdxz, newdyz;
  double dxx, dyy, dzz, dxy, dxz, dyz;

  for ( i = start; i <= end; i++ ) {
    pot1 = 0;
    field1 = 0;
    field2 = 0;
    field3 = 0;
    i3 = i*3;
    dxx = 0.0; 	dxy = 0.0;  dxz = 0.0;
    dyy = 0.0;  dyz = 0.0;  dzz = 0.0;

    for ( j = start; j <= i-1; j++ ) {
      kernel(&FMMLOC[i3], &FMMLOC[3*j], FMMCHARGE[j], &newpot, &newf1, &newf2, &newf3,
 	         &newdxx, &newdyy, &newdzz, &newdxy, &newdxz,&newdyz);
     
      pot1 += newpot;
      field1 += newf1;
      field2 += newf2;
      field3 += newf3;
      dxx += newdxx;
      dyy += newdyy;
      dzz += newdzz;
      dxy += newdxy;
      dxz += newdxz;
      dyz += newdyz;
    }

    for ( j = i+1; j <= end; j++ ) {
      kernel(&FMMLOC[i3], &FMMLOC[3*j], FMMCHARGE[j], &newpot, &newf1, &newf2, &newf3,
	        &newdxx, &newdyy, &newdzz, &newdxy, &newdxz,&newdyz);
      pot1 += newpot;
      field1 += newf1;
      field2 += newf2;
      field3 += newf3;
      dxx += newdxx;
      dyy += newdyy;
      dzz += newdzz;
      dxy += newdxy;
      dxz += newdxz;
      dyz += newdyz;
    }

    fmmpot[i] += pot1;
    fmmfield[i3] += field1;
    fmmfield[i3+1] += field2;
    fmmfield[i3+2] += field3;
    fmmdxx[i] += dxx;
    fmmdyy[i] += dyy;
    fmmdzz[i] += dzz;
    fmmdxy[i] += dxy;
    fmmdxz[i] += dxz;
    fmmdyz[i] += dyz;
  }
}
