#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpyfmm.h>
#include <twobody.h>
#include <mkl.h>

int main(){
    int N, nc = 4,i,j;
    double bradius = 1.0, gap = 4.0, spacing=0.0; 
    // bead size and distance between the two sphere
    //gap: distance between two rigid body;
    //spacing: distance between two beads;
    double *ShellSphs1, *ShellSphs2, ShellRadius;
    
    /*===============================================================================
      Step 1 : Construct the Shell representation
    ==============================================================================*/
    //ShellSphs1 has size 3*nparts but allocated much larger since the number of beads 
    //on each sphere is not known beforehead
    ShellSphs1 = (double*)calloc(pow(nc,3)+3*nc*pow(nc-1,2), sizeof(double));

    N = shell_model(nc, bradius,spacing, &ShellRadius, ShellSphs1);
    
    ShellSphs2 = (double*)calloc(3*N, sizeof(double));
    
    double offset = (2*ShellRadius + gap); // Distance between two sphere
    
    for(i=0; i< 3*N; i++){
    	ShellSphs2[i] = ShellSphs1[i]+ offset;
    }
    double *T, *Q, *b, *y;
    T = (double*)calloc(3*N*6, sizeof(double)); //3N*6
    Q = (double*)calloc(6*3*N, sizeof(double));// 6*3N
    // b = T*F
    b = (double*)calloc(3*N*2, sizeof(double));
    
    //Compute T and Q matrix
    TandQ(T, Q,ShellSphs1,N);

  /*================================================================================
      Step 2 : Given the external Force, compute the force on beads 
   			D*S = T*F = b
                     
                        => x = T'*S = (T'* D^{-1} * T )*F = A*F 
   Note : Here, since the two sphere are with the same no. of beads,T1 = T2 = T.
   So we basically use the same T but with different right hand side F to caculate vb
   or body 1 and body 2.
   ================================================================================*/
    double F[12]={0.0, 0.0, 1.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0};
    // cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);
    double *xx;
    xx = (double *)calloc(3*N, sizeof(double));

    cblas_dgemv(CblasRowMajor, 'n', 3*N, 6, 1, T, 3*N, xx, 1, 0, b,1);
    //1. Solve D*x = b using Krylov subspace method
	S = (double*)calloc(3*N*12,sizeof(double));
	KrylovFMM(S, ShellSphs1, b, 3*N);
	KrylovFMM(&S[6], ShellSphs1,&b[6],3*N);
   /*=============================================================================
     Step 3 : Use iterative method to compute V = Z^(-1)*y
    =============================================================================*/
    ComputeBodyVelocity();
    
    free(ShellSphs1);
    free(ShellSphs2);
    free(T);
    free(Q);
    free(S);
    free(b);
    free(xx);
    return 0;
}
