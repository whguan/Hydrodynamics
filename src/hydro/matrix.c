#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//Torque  Matrix T = T1 = T2
// N is the number of beads in each large sphere,
// ie. the dimension of shellsphs 
int  Tmat(double *T, double *ShellSphs,int N){
	int i;
	
	double *lc_eps;
	lc_eps = (double*)calloc(27,sizeof(double));
	LCepsilon(lc_eps);
	
	for(i = 0; i < N; i++){
		int i3 = i*3;
		double r[3];
		r[0] = ShellSphs[i3];
		r[1] = ShellSphs[i3+1];
		r[2] = ShellSphs[i3+2];
		double *rmat;
		rmat = (double*)calloc(9,sizeof(double));
		rotmat(&r,rmat,lc_eps);
		int kx,ky;
		for(kx=0;kx<3;kx++){
		    for( ky =0; ky<3; ky++){
		    	T[ (i3+kx)*6 + ky  ] = (kx==ky ? 1:0);
			T[ (i3+kx)*6 + (3+ky)]= rmat[3*kx+ky];
		//	Q[ kx*3*N + (i3+ ky)] = (kx==ky? 1:0);
		//	Q[(3+kx)*3*N +(i3 +ky)] = -rmat[3*kx +ky];
		    }
		}
		free(rmat);
	}
	free(lc_eps);
	return 0;

}
// code for Matrix Q can be deleted.
//Matrix Q that converts forces on individual bead to total force/torque
int Qmatrix(double *Q, double *ShellSphs, int N ){
	int i;
	for(i=0; i< N; i++){
		int i3 = i*3;
		double r[3];
		r[0] = ShellSphs[i3];
		r[1] = ShellSphs[i3+1];
		r[2] = ShellSphs[i3+2];
		double *rmat;
		rmat = (double*)calloc(9,sizeof(double));
		rotmat(&r, rmat);
		int kx,ky;
		for(kx=0; kx < 3; kx++){
		    for(ky = 0; ky < 3; ky++){
		    	Q[ kx*3*N + (i3+ ky)] = (kx==ky? 1:0);
			Q[(3+kx)*3*N +(i3 +ky)] = -rmat[3*kx +ky];
		    }
		}

	}	
	return 0;
}

/*==========================================================================
 Levi-Civitia epsilon 
===========================================================================*/
void LCepsilon(double *lc_eps){
	int i,j,k;
	for(i=0;i<3;i++){
	    for(j=0;j<3;j++){
	    	for(k=0; k<3;k++){
			lc_eps[i*9+j*3+k]= (i-j)*(j-k)*(k-i)/2;
		}
	    }
	}
}
/*===========================================================================
   Measure of size of 
 [velocity 1, angular velocity 1, velocity 2, angular velocity 2] vector
 The angular velocity has been scaled by radius so that units are consistent.
============================================================================*/
double sph_norm(double *v, double radius ){
	
	double sum = 0.0;
	int i;
	
	for(i=0; i<12; i++){
		if((i<3)||(i>8)){
			sum += v[i]*v[i];
	        }else{
			sum += v[i]*v[i]*radius*radius;
		
		}

	}
	
	return sqrt(sum);
}

/* ============================================================================
   This matrix is used to take cross products with respect to the bead's
   displacements from the large sphere center.

==============================================================================*/
int  rotmat(double *r, double *M, double *lc_eps){
	int i, j,k;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			double sum = 0.0;
			for(k = 0;k < 3; k++){

			    sum += lc_eps[i*9+j*3+k] *r[k];
			}
		
		M[i*3+j] = sum;
		}	
	}
	return 0;
}

