#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rpyfmm.h"
#include "rpytest.h"

int main(int argc, char **argv)
{
  int nparts, s, accuracy, distribution; 
  double *ploc, *pcharge, *RPY;
  double length, beta, elapsed; 
  struct timeval tic, toc;
  
  test_parser(argc, argv, &length, &beta, &nparts, &s, &accuracy, &distribution);
  
  //allocate for rpy
  ploc = (double *)calloc(3*nparts, sizeof(double));
  pcharge = (double *)calloc(3*nparts, sizeof(double));
  RPY    = (double *)calloc(3*nparts, sizeof(double));
  
  //Generate data for ploc and pcharge
  test_data(length, nparts, distribution, ploc, pcharge);

  //call RPYfmm to calculate RPY
  RPYfmm(length,beta, s, accuracy, nparts, ploc, pcharge, RPY);
  
  //verify RPYfmm
  int ncheck = nparts;
  int i;
  double *drpy;
  drpy   = (double*)calloc(3*ncheck, sizeof(double));
  for(i = 0; i< ncheck; i++){
    rpydirect(nparts, ploc, pcharge, i, &drpy[3*i]);
  }
  double err = 0.0, total =0.0;
  for(i=0; i<ncheck; i++){
  //   printf("%.2f %.2f %.2f \n", RPY[3*i], RPY[3*i+1], RPY[3*i+2]);
  //   printf("%.2f %.2f %.2f \n", drpy[3*i], drpy[3*i+1], drpy[3*i+2]);
  //   printf("\n");
     err += pow(RPY[3*i]-drpy[3*i],2);
     err += pow(RPY[3*i+1]-drpy[3*i+1],2);
     err += pow(RPY[3*i+2]-drpy[3*i+2],2);
     total += pow(drpy[3*i],2)+ pow(drpy[3*i+1],2)+pow(drpy[3*i+2],2);
  }

  err = sqrt(err);
  total = sqrt(total);
  printf("err = %f, total = %f, err/total = %20.4e\n", err, total, err/total);
  
  free(ploc);
  free(pcharge);
  free(RPY);
  free(drpy);
  
  return 0;
}
  
void rpydirect(const int nparts, double *ploc, double *pcharge, const int i, double dpoti[3])
{
	int j;
	dpoti[0]=0; dpoti[1]=0; dpoti[2]=0;
	for(j=0; j< i; j++){
	    double potj[3];
	    potj[0]=0.0; potj[1]=0.0; potj[2]=0.0;
		rpy3sup_eval(&ploc[3*i], &ploc[3*j], &pcharge[3*j], potj);
		dpoti[0] += potj[0];
		dpoti[1] += potj[1];
		dpoti[2] += potj[2];

	}

	 for(j=i+1; j< nparts; j++){
		double potj[3];
		potj[0]=0; potj[1]=0; potj[2]=0;
		rpy3sup_eval(&ploc[3*i], &ploc[3*j], &pcharge[3*j], potj);
		dpoti[0] += potj[0];
		dpoti[1] += potj[1];
		dpoti[2] += potj[2];
         }
	 //added self interaction Aug 25
	  double ha = 1;
	  double done = 1.0;
	  double pi = 4.0*atan(done);
	  double cmu =1.0;
	  double c0 = cmu/(6.0*pi*ha);
	 dpoti[0] += c0 * pcharge[3*i]; 
         dpoti[1] += c0 * pcharge[3*i+1];
	 dpoti[2] += c0 * pcharge[3*i+2];
	
}	
void test_parser(int argc, char ** argv, double *length, double *beta, int *nparts, int *s, 
		 int *accuracy, int *distribution)
{
  // Setup default demo parameters: The kernel function is expressed as 
  // k(r) = exp(-beta*r)/r. When beta is zero, the demo is for Laplace kernel.
  // When beta is a positive number, the demo is for Yukawa kernel. 

  // Three types of distributions are provided for the test: (1) particles uniformly
  // distributed inside a box, referred to as CUBE; (2) particles uniformly distributed 
  // over a spherical surface, referred to as SPHERE; and (3) the restriction of type
  // (2) in the first octant, referred to as OCTANT. 

  *length = 10000;
  *beta = 0; 
  *nparts =5000;
  *s = 80;
  *accuracy = 9;
  *distribution = 1;

  // Parse command line flags when necessary. 
  int i = 0;
  if ( argc > 1 ) {
    for ( i = 1; i < argc; i++ ) {
      if ( argv[i][0] == '-' ) {
	switch ( argv[i][1] ) {
	case 'l':
	  *length = atof(argv[++i]);
	  break;
	case 'b':
	  *beta = atof(argv[++i]);
	  break;
	case 'n':
	  *nparts = atoi(argv[++i]);
	  break;
	case 's':
	  *s = atoi(argv[++i]);
	  break;
	case 'a':
	  *accuracy = atoi(argv[++i]);
	  break;
	case 'd':
	  *distribution = atoi(argv[++i]);
	  break;
	default:
	  break;
	}
      }
    }
  }
  return;
}
void test_data(double length, int nparts, int distribution, double *ploc, double *pcharge)
{
  if ( distribution == 1 ) {
    // generate data randomly distributed in a unit cube
    int i; 
    // wenhua
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      // wenhua Aug 6
      pcharge[j] = 1.0;// (double) rand()/RAND_MAX - 0.5;
      pcharge[j+1] =1.0;
      pcharge[j+2] =1.0;

      ploc[j] = ((double) rand()/RAND_MAX - 0.5)*length;
      ploc[j+1] = ((double) rand()/RAND_MAX - 0.5)*length;
      ploc[j+2] = ((double) rand()/RAND_MAX - 0.5)*length;
    }
  } else if ( distribution == 2 ) {
    // generate data randomly distributed over a spherical surface
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      pcharge[j] = (double) rand()/RAND_MAX - 0.5;
	  pcharge[j+1] = (double) rand()/RAND_MAX - 0.5;
	  pcharge[j+2] = (double) rand()/RAND_MAX - 0.5;
      double theta = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      double phi = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      ploc[j] = sin(theta)*sin(phi);
      ploc[j+1] = sin(theta)*cos(phi);
      ploc[j+2] = cos(theta);
    }
  } else if ( distribution == 3 ) {
    // generate data randomly distributed over a spherical surface that is in the first octant
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      pcharge[j] = (double) rand()/RAND_MAX - 0.5;
	  pcharge[j+1] = (double) rand()/RAND_MAX - 0.5;
	  pcharge[j+2] = (double) rand()/RAND_MAX - 0.5;
      double theta = ((double) rand()/RAND_MAX)*M_PI_2;
      double phi = ((double) rand()/RAND_MAX)*M_PI_2;
      ploc[j] = sin(theta)*sin(phi);
      ploc[j+1] = sin(theta)*cos(phi);
      ploc[j+2] = cos(theta);
    }
  }

  return;
}
