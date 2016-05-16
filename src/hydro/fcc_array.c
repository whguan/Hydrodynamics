/*This program is to generate a face-centered cubic lattice with nc
  beads in each direction

  By Wenhua Guan, Jan 27.2015.
  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//The number of elements in sphs is (nc)^3+(nc-1)^2*nc*3
int  fcc_array(const int nc, const double bradius, const double spacing, double *sphs)
{	
	double d = 2 * bradius + spacing ;
	double ld = sqrt(2) * d ;
	double hd = ld/2;
	double c = ( nc - 1 ) * hd ;

	int ix, iy, iz;

	int  num1 =3* pow(nc,3);

	for (ix = 0; ix < nc; ix++){
	     for (iy = 0; iy < nc; iy++){
		  for(iz = 0; iz < nc; iz++){
			int num = ix*nc*nc + iy*nc + iz;
			sphs[3*num] = ix * ld - c;
			sphs[3*num +1] = iy * ld - c;
			sphs[3*num +2] = iz * ld - c;			
		  }
	     }
	}

	for(ix = 0; ix < nc-1; ix++){
	    for (iy = 0; iy < nc-1; iy++){
		for (iz = 0; iz < nc; iz++){
		    int num = num1 + 9*(ix*(nc-1)*nc + iy*nc +iz);
		    double sx = ix * ld ;
		    double sy = iy * ld ;
		    double sz = iz * ld ;

		    sphs[num] = sx + hd - c;
		    sphs[num + 1] = sy + hd - c;
		    sphs[num + 2] = sz - c;
                     
		    sphs[num + 3] = sz - c;
		    sphs[num + 4] = sx + hd - c;
		    sphs[num + 5] = sy + hd - c;

		    sphs[num + 6] = sy + hd - c;
		    sphs[num + 7] = sz - c;
		    sphs[num + 8] = sx + hd - c;

		}
	    }
	}

        return 0;

}

//Shell_model is to generate a large sphere with just the surface made of beads
int shell_model(const int nc, const double bradius, const double spacing, double *ShellRadius, double *ShellSphs)
{
   int num = pow(nc,3) + 3 *nc * (nc-1)*(nc-1);

   *ShellRadius = sqrt(2)/2* (2*bradius+spacing) *(nc-1); //? spacing? 
   double *sphs;
   sphs = (double*)calloc(3*num, sizeof(double));

   fcc_array(nc, bradius, spacing, sphs);
   
   int i;
   int count = 0;
   for(i = 0; i < num ; i++ ){
   	double pt[3]={0.0,0.0,0.0};
	pt[0] = sphs[3*i];
	pt[1] = sphs[3*i+1];
	pt[2] = sphs[3*i+2];
	
	double r = pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];
	r = sqrt(r);
	
	if((r > (*ShellRadius - bradius))&&( r < (*ShellRadius + bradius))){
		ShellSphs[3*count] = pt[0];
		ShellSphs[3*count + 1] = pt[1];
		ShellSphs[3*count + 2] = pt[2];
		count++;
	
	}
   }

   free(sphs);
   
   return count;
}
