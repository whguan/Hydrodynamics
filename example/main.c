/*=========================================================================
  This program is for modeling the hydrodynamic interaction of two shell-
  representated sphere.
  By Wenhua Guan Feb 10,2015
  ========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <RPY.h>
#include <twobody.h>
#include <mkl.h>

int main(){
    int N, nc = 20,i;
    double bradius = 1.0, gap = 4.0; // bead size and distance between the two sphere
    double *ShellSphs,*ShellSphs1, ShellSphs2, ShellRadius;
    
  /*===============================================================================
      Step 1 : Construnt the Shell representation
   ==============================================================================*/
    //ShellSphs1 has size 3*nparts but allocated much larger since the nparts 
    //is not known beforehead
    ShellSphs = (double*)calloc(pow(nc,3)+3*nc*pow(nc-1,2), sizeof(double));

    N = shell_model(nc, bradius, ShellRadius, gap, ShellSphs);
    ShellSphs1 = (double*)calloc(3*N, sizeof(double));
    ShellSphs2 = (double*)calloc(3*N, sizeof(double));
    
    double offset = (2*ShellRadius + gap); // Distance between two sphere
    for(i=0; i< 3*N; i++){
    	ShellSphs1[i] = ShellSphs[i];
    }
    free(ShellSphs);
    for(i=0; i< 3*N; i++){
    	ShellSphs2[i] = ShellSphs1[i]+ offset;
    }
    
    // Allocate memory for matrix T, Q 
    double *T, *Q;
    T = (double*)calloc(3*N*6, sizeof(double));
    Q = (double*)calloc(6*3*N, sizeof(double));

    Torque(T, ShellSphs1,N);
    Qmatrix(Q,ShellSphs1, N )
   
  /*================================================================================
      Step 2 : Given the external Force, compute the force on beads 
                      bf = R*T*(QRT)^(-1) * F   <==> S = R * T
		      				     M = Q * S
						     M * Z = F => Z = M^(-1)F
						     C = T*Z
						     D * bf = C  =>  bf = R*T*Z 
   Note : Here, since the two sphere are with the same no. of beads,T1 = T2 = T.
   So we basically use the same T but with different right hand side F to caculate vb
   or body 1 and body 2.
   ================================================================================*/
    double F1[6]={0.0, 0.0, 1.0,  0.0, 0.0, 0.0};
    double F2[6]={0.0, 0.0, 0.0,  0.0, 0.0, 0.0};
	
    	//1. Solve D1*S = T1 using Krylov subspace method
	double *S, *tmp, *tmprhs;
	S = (double*)calloc(3*N*6,sizeof(double));
	tmpsol = (double*)calloc(3*N,sizeof(double));
	tmprhs = (double*)calloc(3*N,sizeof(double));
	for(int i=0; i<6; i++){
		for(int j = 0; j<3*N, j++) tmprhs[j] = T[6*j+i];
		KrylovFMM(tmpsol, ShellSphs1, tmprhs, 3*N);
		for(int j = 0; j<3*N, j++) S[6*j+i] = tmpsol[j];
	}
	free(tmpsol);
	
	//2. Compute Q*S using matrix multiplication of blas
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
        //         m, n, k, alpha, A, k, B, n, beta, C, n);
	double *M;
	M =(double*)calloc(6*6,sizeof(double));
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                 6, 6, 3*N, 1, Q, 3*N, S, 6, 0, M, 6);
	
	//3. Compute the linear equation of 6*6 matrix Z = M^(-1)F 
	//Here, notice that Z contains Z1= M^(-1)F1 and Z2 = (M)^(-1)F2

	double Z[12] = {0.0, 0.0, 1.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0};
	int info,ipiv[6];
	//dgetrs('N',6,2,M,lda, ipiv, Z, ldb, info);
	// dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
			                // double* b, int* ldb, int* info );
	dgesv(6,2,M,6,ipiv, Z, 6,info);
	if( info > 0 ) {
	         printf( "The diagonal element of the triangular factor of A,\n" );
		 printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
		 printf( "the solution could not be computed.\n" );
		 exit( 1 );
	 }
	//4. Compute the matrix vector multiplication C = T*Z 
	 // Again, C = [C1, C2] C1 = T*Z1,C2 = T*Z2 
	 double *C;
	 C = (double*)calloc(3*N*2,sizeof(double));
	 cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		3*N, 2, 6, 1, T, 6, Z, 6, 0, C, 6);
	
	//5. Use CG method to caculate the solution of D1*fb = T * (Q*S)^(-1)*F = C
         double *fb1,*fb2;
	 *fb1 = (double *)calloc(3*N, sizeof(double)); 
	 *fb2 = (double *)calloc(3*N, sizeof(double));
	 //*tmprhs = (double *)calloc(3*N, sizeof(double));
	 
	 for(int j=0; j< 3*N; j++) tmprhs[j] = C[ j];
	 KrylovFMM(fb1, ShellSphs1, tmprhs, 3*N); 

	 for(int j= 0; j<3*N; j++) tmprhs[j] = C[3*N+j];
         KrylovFMM(fb2, ShellSphs1, tmprhs, 3*N);
	 

   /*===============================================================================
      Step 3 : Compute bead velocities : vb = D* fb through FMM
    ===============================================================================*/
        double *vb, *vold, *vnew;
	vb = (double*)calloc(3*N*2, sizeof(double));
	vnew = (double*)calloc(12, sizeof(double));
	vold = (double*)calloc(12, sizeof(double));
	
	for(int i =0; i<12; i++) vold[i]=0;
	int iter=0;
	int maxiter = 20;
     	while(iter< maxiter){
		iter++;
		RPY(N, ShellSphs1, fb1,  tmprhs );
		for(int j =0; j<3*N; j++) vb[j] = tmprhs[j];
		
		RPY(N, ShellSphs1, fb2,  tmprhs );
		for(int j =0; j<3*N; j++) vb[3*N+j]= tmprhs[j];
   /*===============================================================================
      Step 4 : Compute rigid body motions through solving  least square problem 
      		     v = P * vb   <==> min|| T*v - vb ||
     Here, T is with size (3*N) * 6. QR/LQ  factorization is applied to caculate the least
     square problem for full rank matrix. 
	dgels  : Overdetermined with full rank matrix A using QR factorization;
	dgelss : Minimun solution with SVD 
     vb is overwritten during the calculation.
     If m > n, vb[1:n] is the solution v.
   ===============================================================================*/
    // void dgels( char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
    //		                    double* b, int* ldb, double* work, int* lwork, int* info );
   
    //dgelss ( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info );
     		double rcond = -1.0, s,wkopt, *work;
		int rank, lwork = -1, m = 3*N, n =6,nrhs =2, lda = m, ldb = m;
     /* Query and allocate the optimal workspace */
     		dgels("No transpose" &m, &n, &nrhs, T, &lda, vb, &ldb, &wkopt, &lwork, &info );
    		lwork = (int)wkopt;
     		work = (double*)calloc(lwork, sizeof(double));
     	//	dgelss( 3*N, 6, 2, T, 3*N, vb, 3*N, s, rcond, rank, work, lwork, info );
                dgels("No transpose" &m, &n, &nrhs, T, &lda, vb, &ldb, &work, &lwork, &info );
     /* Check for convergence */
     		if( info > 0 ) {
			printf( "The algorithm computing SVD failed to converge;\n" );
	        	printf( "the least squares solution could not be computed.\n" );
	     		exit( 1 );
     		}

		for(int i = 0; i< nrhs; i++){
			for(int j = 0; j<n; j++){
			    vnew[j] = vb[i*n +j];
			}
		}

	        free(work);
   /*===============================================================================
      Step 5 : Compare v with the previous values; if close, then done; otherwise, 
      project the rigid body components of bead velocities and compute bead forces 
      needed to counter remaining motions.
                     fb <-- fb - R*(vb- T*v)  <==>  D * delta_fb = vb - T*v
		                                    fb = fb - delata_fb
    ==============================================================================*/
                double *deltav, eps = 1.0e-3;
		int incx =1;
		deltav = (double *)calloc(12, sizeof(double));
		for(int i =0; i< 3*N*2; i++) deltav[i] = vnew[i] -vold[i];
		double residual = cblas_dnrm2(12, deltav, incx);
	        // Check the difference between the new and the old velocities for the body
		if(residual<eps){
			break;
		}
                // Continue to the for loop  T*vnew - vb
		// Call dgemv to compute T*vnew
		CBLAS_LAYOUT layout =  CblasRowMajor;
		CBLAS_TRANSPOSE trans = CblasNoTrans; 
		double *deltfb, *deltvb;
		deltfb1 = (double*)calloc(3*N, sizeof(double));
		deltfb2 = (double*)calloc(3*N, sizeof(double));
		deltvb1 = (double*)calloc(3*N*2, sizeof(double));
		cblas_dgemv(layout, trans, m, n, 1, T, lda, vnew, 1, 0, deltvb, incy);
		
		// Solve D*deltfb1 = T*vnew - vb
	        for(int j = 0; j<3*N, j++) tmprhs[j] = deltvb[j];
		KrylovFMM(deltfb1, ShellSphs1, tmprhs, 3*N);
		for(int j = 0; j<3*N, j++) tmprhs[j] = deltvb[3*N+j];
		KrylovFMM(deltfb2, ShellSphs2, tmprhs, 3*N);

		for(int i = 0; i<3*N){
	             fb1[i] += deltfb1[i];
		     fb2[i] += deltfb2[i];
		}

	}
	free(ShellSphs1);
	free(ShellSphs2);
	free(T);
	free(Q);
	free(S);
	free(M);
	free(C);
	free(tmprhs);
	free(fb1);
	free(fb2);
	free(vb);
	free(vold);
	free(vnew);
	free(deltv);
	free(deltfb1);
	free(deltfb2);
	free(deltvb);
	return 0;
}
