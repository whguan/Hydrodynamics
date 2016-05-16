#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <math.h>
#include "adap_fmm.h"
#include "rpy.h"
/*********************************************************************
rpyfmm -- kernel for rotne-prager model D*v
          $ Far field is calculated through 4 calls of laplace FMM
          $ Near field is directly evaluated by rpydirect routines

Input : nparts  --- number of particles 
        ploc    --- location of particles with size 3*nparts
		pcharge --- charge of particles with size 3*nparts

Output : RPY    --- potential for D*v with size 3*nparts

***********************************************************************/
int RPYfmm(double length, double  beta, int s, int accuracy, int nparts, 
	   double *ploc, double *pcharge, double *RPY){
  /*********************************************************************
  //step 1:  Set all parameters 
  **********************************************************************/
  //Laplace FMM
  double *pot, *field, *charge;
  double *dxx, *dyy, *dzz, *dxy, *dyz, *dxz;
  
  //rpy
  double ha, pi, c0, c1, cmu, c2;
  double elapsed;
  struct timeval tic, toc;
  
  // ha - raids of the particle
  ha = 1.0;// 0.1;
  pi = 4*atan(1.0); 
  
  cmu = 1.0; //common factor 
  c0 = cmu/(6.0*pi*ha); // Coeff of D_ii
  c1 = cmu/(8.0*pi);
  c2 = cmu*ha*ha/(12.0*pi);
  
  //set whether displacement and/or gradient to be computed
  int ifpot, ifgrad, ifhess, ifdirect;

  /*************************************************************************
    step 2: Evaluate D*v via Laplace FMM
    	    1) Far field contribution through 4 laplace fmm calls
            2) Direct evaluation for local interaction of RPY model 
  *************************************************************************/
  test_init( length, ha, beta, s,  accuracy, nparts,  &charge, &pot, &field,
	     &dxx, &dyy, &dzz, &dxy, &dxz, &dyz);
  
//  printf("\n\tPROGRESS\n");
//  printf("======================================================\n");
  gettimeofday(&tic, 0);
     
/*************************************************
    step 2.1: far field - 4 laplace fmm calls
 *************************************************
    step 2.1.1 First three FMM calling
 *************************************************/
  // FMM Initialization
  int i, m;
  ifpot = 1; ifgrad =1; ifhess =1; ifdirect = 0; 
  
  // Direct evaluation for rpy is just caculated on the 4th fmm calling
  // The direct evalution is ignored in the first three fmm calls
  adap_fmm_init(accuracy, nparts);
  adap_fmm_graph(nparts, s, ploc); 
  for(i = 0; i < 3; i++){
        for(m = 0; m < nparts; m++){
		pot[m]=0;
		field[3*m] =0; field[3*m+1] =0; field[3*m+2] = 0;
		dxx[m] = 0; dyy[m] = 0; dzz[m] = 0;
		dxy[m] = 0; dxz[m] = 0; dyz[m] = 0;
		charge[m] = pcharge[3*m + i];//charge[i*nparts + m]; //July 9th
        }   
//	printf("The %dth fmm call\n ", i+1);
        
	// assign charge 
	adap_fmm_charge(ifdirect, charge, pcharge);
	
    	// FMMS compute 
	adap_fmm_compute(ifdirect, ifpot, ifgrad, ifhess);
  
    	// FMM POST
	adap_fmm_post(ifdirect, RPY, ifpot, pot, ifgrad, field, ifhess, dxx, dyy, dzz, dxy, dxz, dyz);
	
	// Update for RPY
	for(m = 0; m < nparts; m++)
	{   	  
	    RPY[3*m+i] += c1* pot[m]; // pot[i * nparts + m] +=  c1* ppot[m];
            //march 31 change "-" to "+" for field   
	    RPY[3*m ]  += c1 * ploc[3*m +i] * field[3*m ];    
	    RPY[3*m + 1] += c1 * ploc[3*m + i] * field[3*m + 1];   
	    RPY[3*m + 2] += c1 * ploc[3*m + i] * field[3*m + 2];    
	
		if(i==0){
	
			RPY[3*m ]    -= c2 * dxx[m];
			RPY[3*m + 1] -= c2 * dxy[m];
			RPY[3*m + 2] -= c2 * dxz[m];
		
		}else if(i==1){
	
			RPY[3*m ]    -= c2 * dxy[m];
			RPY[3*m + 1] -= c2 * dyy[m];
			RPY[3*m + 2] -= c2 * dyz[m];
	
		}else{
	
			RPY[3*m ]    -= c2 * dxz[m];
			RPY[3*m + 1] -= c2 * dyz[m];
			RPY[3*m + 2] -= c2 * dzz[m];
		} 
	}

 }
  
 //*************************************************
 //   step 2.1.2: The 4th FMM calling
//**************************************************
    ifpot =0; ifgrad =1; ifhess=0;
    //update the charge
  //  printf("The 4th Laplace FMM call\n");
    for(m = 0; m < nparts; m++) {
	
	    charge[m] = 0;
	    pot[m] = 0;
	    field[3*m] = 0;  field[3*m+1] = 0;  field[3*m+2] = 0;
	    dxx[m] = 0;      dyy[m] = 0;        dzz[m] = 0.0;
	    dxy[m] = 0.0;    dxz[m] = 0.0;      dyz[m] = 0.0;
		
	    for(i = 0; i < 3; i++){
		
               charge[m]  -= c1* ploc[3*m + i] *pcharge[3*m + i];  //charge[i*nparts + m];
	    }   
	   
	   
     }
     // assign charge 
     adap_fmm_charge(ifdirect, charge, pcharge);
	
    // FMMS compute 
     adap_fmm_compute(ifdirect,ifpot, ifgrad, ifhess);
  
    // FMM POST
      adap_fmm_post(ifdirect, RPY, ifpot, pot, ifgrad, field, ifhess, dxx, dyy, dzz, dxy, dxz, dyz);

     // Update RPY       
     for(m = 0; m < nparts; m++){
	   for(i=0; i<3; i++){
		// April 24th change "+" to "-"
		RPY[3*m+i] += field[3*m+i];
	   }
//	   printf("%f %f %f \n", RPY[3*m], RPY[3*m+1], RPY[3*m+2]);

     }

    /*************************************************
    step 3: Near field - RPY kernel
    *************************************************/
    ifdirect = 1; ifpot =1; ifgrad =0; ifhess=0; 
  //  printf("Near field for rpy\n");
    // assign charge 
    adap_fmm_charge(ifdirect, charge, pcharge);
	
    // FMMS compute 
    adap_fmm_compute(ifdirect, ifpot, ifgrad, ifhess);
    
    // FMM POST
    adap_fmm_post(ifdirect, RPY,ifpot, pot, ifgrad, field, ifhess, dxx, dyy, dzz, dxy, dxz, dyz);
    
    // FMM CLEAN
    adap_fmm_clean(); 
    

    test_clean( charge, pot, field, dxx, dyy, dzz, dxy, dxz, dyz);
    gettimeofday(&toc, 0);
    elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  //  printf("t(RPY COMPUTING)          | %20.4e\n", elapsed);

    return 0;	
}


void test_init( const double length, const double ha, const double beta, 
		const int s, const int accuracy,const int nparts,
		double **charge, double **pot, double **field, 
		double **dxx, double **dyy, double **dzz,
	        double **dxy, double **dxz, double **dyz)
{
  // The function completes two tasks: (1) Allocate memory to hold particle locations
  // and the charges carried by them, and to hold the computed potential and field 
  // result; (2) Generate the input data as specified by the distribution variable. 


  // Allocate memory to hold particle information and output results. 
  (*charge) = (double *)calloc(nparts, sizeof(double));
  (*pot)     = (double *)calloc(nparts, sizeof(double));
  (*field)   = (double *)calloc(3*nparts, sizeof(double));

  (*dxx) = (double *)calloc(nparts, sizeof(double));
  (*dyy) = (double *)calloc(nparts, sizeof(double));
  (*dzz) = (double *)calloc(nparts, sizeof(double));
  (*dxy) = (double *)calloc(nparts, sizeof(double));
  (*dxz) = (double *)calloc(nparts, sizeof(double));
  (*dyz) = (double *)calloc(nparts, sizeof(double)); 
 
   int allocFailure = (*pot==0) || (*field==0) || (*charge ==0) ||
  (*dxx == 0) || (*dyy == 0) ||(*dzz == 0) ||(*dxy == 0) || (*dxz == 0) || (*dyz == 0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n",  __FILE__, __LINE__);
    exit(-1);
  }

  // Generate input data set for the demo
//  test_data(nparts, distribution, *ploc, *pcharge);

  // Print summary of the demo
/*  if ( beta > 0 ) {
    printf("\n\tDEMO FMM-YUKAWA RUN\n"
	   "\n\tSETUP\n\n"
	   "======================================================\n"
	   "Length of Cube       | %20.2f\n"
	   "# OF PARTICLE        | %20d\n"
	   "Radius of PARTICLE   | %20.2f\n"
	   "S                    | %20d\n"
	   "ACCURACY             | %20d\n"
	   "======================================================\n", 
	   length, nparts, ha, s, accuracy);
  } else if  ( beta == 0 ) {
    printf("\n\tDEMO FMM-LAPLACE RUN\n"
	   "\n\tSETUP\n\n"
	   "======================================================\n"
	   "Length of Cube       | %20.2f\n"
	   "# OF PARTICLE        | %20d\n"
	   "Radius of PARTICLE   | %20.2f\n"
	   "S                    | %20d\n"
	   "ACCURACY             | %20d\n"
	   "======================================================\n", 
	   length, nparts, ha, s, accuracy);
  }*/
  return;
}
void test_clean( double *charge, double *pot, double *field, 
		double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz )
{
  free(charge);
  free(pot);
  free(field);
  free(dxx);
  free(dyy);
  free(dzz);
  free(dxy);
  free(dxz);
  free(dyz);
}
