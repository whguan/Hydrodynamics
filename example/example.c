/*
  example.c: demonstrates the use of parallel adaptive fmm-laplace package
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

#include <stdlib.h>
#include <math.h>
#include "adap_fmm.h"

int main (int argc, char **argv) {
  int nparts, s, accuracy; 

  double beta, *ploc, *pcharge, *pot, *field; 

  beta = 0;    // laplace kernel
  nparts = 10000; // number of particles
  accuracy = 6;   // number of digits of accuracy (3 or 6 supported) 
  s = 80;         // max number of particles per partition box

  // allocate memory to hold particle information and output results
  ploc =    (double *)calloc(nparts*3, sizeof(double));
  pcharge = (double *)calloc(nparts, sizeof(double));
  pot =     (double *)calloc(nparts, sizeof(double));
  field =   (double *)calloc(nparts*3, sizeof(double));

  int i;
  double sqrtn = sqrt((double) nparts); 
  double PI = 3.14159;
  for (i=0; i<nparts; i++) {
    double x,y,z,r;

    r = 1.0 / (double) nparts * (double) i;
    x = cos(2.0*PI*i/sqrtn) * r; 
    y = sin(2.0*PI*i/sqrtn) * r; 
    z = sin(2.0*PI*i/sqrtn) * cos(2.0*PI*i/sqrtn) * r;

    ploc[3*i    ] = x;
    ploc[3*i + 1] = y;
    ploc[3*i + 2] = z;
    pcharge[i] = 1.0;
  }

  adap_fmm_init(accuracy, nparts);
  
  // Use in a time-marching scheme
  int t0 = 0, tn = 1, t;
  for ( t = t0; t < tn; t++ ) {
    // Construct the graph
    adap_fmm_graph(nparts, s, beta, ploc, pcharge);

    // Compute potential and field
    adap_fmm_compute();

    // Post processing 
    adap_fmm_post(pot, field);

    // Update particle locations
    // update_position(); 
  }

  adap_fmm_clean();


  free(ploc); free(pcharge); free(pot); free(field);

  return 0;
}
