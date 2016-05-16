#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include "mkl.h"
int BCG_FMM(const double bradius, const double length, const double beta, const int s, 
	      const int accuracy, const int  nparts, double *ShellSphs, 
	      const double *T1, const double *rhs, double *solution)
{ 
  double *Amat; // the approximation of the matrix T'D^{-1}*T
  Amat = (double *)calloc(144, sizeof(double)); 
  int i, n = 3*nparts +12, numiter =0;
  int *ipar;
  double *fpar, *uvec, *rwork;
  
  ipar = (int *)calloc(16, sizeof(int));
  fpar = (double *)calloc(16, sizeof(double));

  ipar[0] = 0;  // To initialize the iterative solver. 
  ipar[1] = 0;  // Right preconditioning.
  ipar[2] = 2;  // Stopping criteria.
  ipar[4] = 100; // Krylov subspace can be as large as 100. 
  
  fpar[0] = 1e-8; // Relative tolerance 
  fpar[1] = 1e-8; // Absolute tolerance 
  fpar[10] = 0.0; // Reset flops
  
  ipar[3]  = (n+3)*(ipar[4]+2) +(ipar[4] +1)*ipar[4]/2;
  
  rwork = (double*)calloc(ipar[3], sizeof(double));
  uvec = (double *)calloc(n, sizeof(double));
  if (rwork == 0 || uvec == 0) {
    printf("memory allocation failure");
    exit(-1);
  }
//  for(i=0;i<12;i++) printf("F[6N+%d]= %3f\n", i, rhs[3*nparts+i]);
  /*---------------------------------------------------------------------------*/
  /* Initialize the initial guess                                              */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    solution[i] = 1.E-4;
  
  while(1) {
      gmres_(&n, rhs, solution, ipar, fpar, rwork);

      if(ipar[0]==1){
      	 if(numiter){
	         printf("...Iteration %2d residual norm : %30.5e\n", numiter, fpar[4]);
	 }
	 numiter++; 
         double *ptr = &rwork[ipar[7] - 1];
	 for (i = 0; i < n; i++) uvec[i] = ptr[i]; 
         double *V = &rwork[ipar[8] - 1];
	 GetRightSide(bradius, length, beta, s, accuracy, nparts, ShellSphs, uvec, T1, V);
      }else if (ipar[0] >= 2) {
         printf("nonsupported matrix-vector multiply requested by the GMRES");
         exit(-1);
      }else {
         printf("... Iteration %2d residual norm: %38.5e\n"
	     "%-50s%20d\n"
	     "----------------------------------------------------------------------\n",
	     numiter, fpar[4], 	     
	     "... Iterative solver exits with status:", ipar[0]); 
      break;
     }
  }

  free(rwork);
  free(uvec);
  free(ipar);
  free(fpar);
  return 0;
}
// Below is to compute the matrix-vector multiplication
// [D  -T]  * [ v1 ] = [ D*v1 - T*v2 ]
// [T' 0]    [ v2 ]   [ T'*v1       ]
// cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, 
// const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);

void GetRightSide(const double bradius, const double length, const double beta, const int s, const int accuracy, 
		  const int nparts, double *ShellSphs, double *uvec,  const double *T1, double *V )
{
    int n = 3*nparts+12,i;
    double *u; //u=T*v2
    u = (double *)calloc(3*nparts, sizeof(double));

    //Compute V[1:3*nparts]= D*uvec[1:3*nparts]
    RPYfmm(length, beta, s, accuracy, nparts, ShellSphs, uvec, V);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3*nparts/2, 6, 1.0, T1, 6, &uvec[3*nparts], 1, 0, u,1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3*nparts/2, 6, 1.0, T1, 6, &uvec[3*nparts+6], 1, 0, &u[3*nparts/2],1);
   // for(i=0;i<3*nparts; i++) 
	//    printf("u[%d]= %f\n ", i, u[3i]);
    for(i=0;i<3*nparts; i++)  
	    V[i] = V[i]- u[i];
    cblas_dgemv(CblasRowMajor, CblasTrans, 3*nparts/2, 6, 1.0, T1, 6, uvec, 1, 0, &V[3*nparts],1);
    cblas_dgemv(CblasRowMajor, CblasTrans, 3*nparts/2, 6, 1.0, T1, 6, &uvec[3*nparts/2], 1, 0, &V[3*nparts+6],1);

    free(u);
}
