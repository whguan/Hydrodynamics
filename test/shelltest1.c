#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpyfmm.h>
#include <mkl.h>

//bradius: bead size 
//gap: distance between two rigid body;
//spacing: distance between two beads;
//N : The number of the beads on each body
int main(){
    int N, nc = 4,i,j;
    double bradius = 1.0; //1.0/(2*sqrt(2)),
    double gap = 20.0, spacing=0.0; 
	
    double beta = 0.0;
    int accuracy = 12, s = 80;
    int size =pow(nc,3)+3*nc*pow(nc-1,2); 
    size = size *6;
    double *ShellSphs, ShellRadius; 
    ShellSphs = (double*)calloc(size, sizeof(double));
    /*===============================================================================
      Step 1 : Construct the Shell representation
    ==============================================================================*/
    //Construct sphere 1
    N = shell_model(nc, bradius,spacing, &ShellRadius, ShellSphs);
    double length = 2*ShellRadius;
    
    printf("Sphs1, N=%d,ShellRadius = %f\n", N, ShellRadius);    
    //for (i=0;i<N; i++)
    //	 printf("%.2f, %.2f, %.2f\n", ShellSphs[3*i], ShellSphs[3*i+1], ShellSphs[3*i+2]);

    // Construct Sphere 2 from sphere 1
    double offset = (2*ShellRadius + gap); // Distance between two sphere
    printf("offset = %f\n", offset); 
    for(i=0; i< N; i++){
	ShellSphs[3*(N+i)] = ShellSphs[3*i];
	ShellSphs[3*(N+i)+1] = ShellSphs[3*i+1];
    	ShellSphs[3*(N+i)+2] = ShellSphs[3*i+2]+ offset;
    }
	
    printf("Sphs2 is \n");
    for(i=N;i<2*N; i++)
	printf("%.2f, %.2f, %.2f\n", ShellSphs[3*i], ShellSphs[3*i+1], ShellSphs[3*i+2]);

    // Matrix computation
    double *T1;
    T1 = (double*)calloc(3*N*6, sizeof(double)); //3N*6
    
    //Compute T matrix
    Tmat(T1, ShellSphs, N);  
 
  /*================================================================================
      Step 2 : Given the external body Force, compute the body velocity 
   		Assume A = T'*D^{-1}*T, 
		       T = ( T1
		                T2 )  
		compute AV = F 
                     
   Note : Here, since the two sphere are with the same no. of beads,T1 = T2 = T.
   ================================================================================*/
    double *F0, *V0, *F,*V;
    F0 = (double *)calloc(6*N+12, sizeof(double));
    V0 = (double *)calloc(6*N+12, sizeof(double));
    F = (double *)calloc(12, sizeof(double));
    V = (double *)calloc(12, sizeof(double));
    
    for (i=0; i<6*N+12; i++){
	    F0[i]=0.0;
    }
    F[2] = 1.0;

    F0[6*N+2] = 1.0; 

    GMRES_FMM(bradius, length, beta, s, accuracy, 2*N, ShellSphs, T1, F0, V0);
    BCG_FMM(bradius, length, beta, s, accuracy, 2*N, ShellSphs, T1, F, V);

    for (i=0;i<12; i++){
       printf("%f\n",V0[6*N+i]);
    }

    
    free(ShellSphs);
    free(T1);
    free(F0);
    free(V0);
    free(F);
    free(V);
    return 0;
}
