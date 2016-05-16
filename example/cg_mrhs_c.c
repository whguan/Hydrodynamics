#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

/*---------------------------------------------------------------------------*/
/*  Example program for solving symmetric positive definite system of equations.*/
/*  Simplest case: no preconditioning and no the user-defined stopping tests.*/
/*---------------------------------------------------------------------------*/
int CGFMMmrhs(MKL_INT *solution, double *ShellSphs, double *rhs, double *nRhs, MKL_INT n)
{
  /*---------------------------------------------------------------------------*/
  /* Define arrays for the upper triangle of the coefficient matrix and rhs vector */
  /* Compressed sparse row storage is used for sparse representation           */
  /*---------------------------------------------------------------------------*/
  MKL_INT rci_request, expected_itercount = 20, i, j;
  MKL_INT itercount[nRhs];
  /* Fill all arrays containing matrix data. */
  /*---------------------------------------------------------------------------*/
  /* Allocate storage for the solver ?par and temporary storage tmp            */
  /*---------------------------------------------------------------------------*/
  MKL_INT length = 128, method = 1;
  /*---------------------------------------------------------------------------*/
  /* Some additional variables to use with the RCI (P)CG solver                */
  /*---------------------------------------------------------------------------*/
  MKL_INT ipar[128 + 2 * nRhs];
  double euclidean_norm, dpar[128 + 2 * nRhs]; 
  double *tmp;
  tmp = (double*)calloc(n * (3 + nRhs),sizeof(double));
  double eone = -1.E0;
  MKL_INT ione = 1;

  /*---------------------------------------------------------------------------*/
  /* Initialize the initial guess                                              */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n*nRhs; i++)
    solution[i] = 1.E0;
  /*---------------------------------------------------------------------------*/
  /* Initialize the solver                                                     */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < (length + 2 * nRhs); i++)
    ipar[i] = 0;
  for (i = 0; i < (length + 2 * nRhs); i++)
    dpar[i] = 0.E0;

  dcgmrhs_init (&n, solution, &nRhs, rhs, &method, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Set the desired parameters:                                               */
  /* LOGICAL parameters:                                                       */
  /* do residual stopping test                                                 */
  /* do not request for the user defined stopping test                         */
  /* DOUBLE parameters                                                         */
  /* set the relative tolerance to 1.0D-5 instead of default value 1.0D-6      */
  /*---------------------------------------------------------------------------*/
  ipar[8] = 1;
  ipar[9] = 0;
  dpar[0] = 1.E-5;
  /*---------------------------------------------------------------------------*/
  /* Compute the solution by RCI (P)CG solver without preconditioning          */
  /* Reverse Communications starts here                                        */
  /*---------------------------------------------------------------------------*/

rci:dcgmrhs (&n, solution, &nRhs, rhs, &rci_request, ipar, dpar, tmp);
  /*---------------------------------------------------------------------------*/
  /* If rci_request=0, then the solution was found with the required precision */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 0)
    goto getsln;
  /*---------------------------------------------------------------------------*/
  /* If rci_request=1, then compute the vector A*tmp[0]                        */
  /* and put the result in vector tmp[n]                                       */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 1)
    {
      //mkl_dcsrsymv (&tr, &n, a, ia, ja, tmp, &tmp[n]);
	    //debug
       RPY(n, ShellSphs,tmp); // SPMV by 4 calls of  FMM
       goto rci;
    }
  /*---------------------------------------------------------------------------*/
  /* If rci_request=anything else, then dcg subroutine failed                  */
  /* to compute the solution vector: solution[n]                               */
  /*---------------------------------------------------------------------------*/
  goto failure;
  /*---------------------------------------------------------------------------*/
  /* Reverse Communication ends here                                           */
  /* Get the current iteration number into itercount                           */
  /*---------------------------------------------------------------------------*/

getsln:dcgmrhs_get (&n, solution, &nRhs, rhs, &rci_request, ipar, dpar, tmp, itercount);

  /*---------------------------------------------------------------------------*/
  /* Print solution vector: solution[n] and number of iterations: itercount    */
  /*---------------------------------------------------------------------------*/
  printf ("The system has been solved\n");
  printf ("The following solution obtained\n");
  for (i = 0; i < n / 2; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\n");
  for (i = n / 2; i < n; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\n");
  for (i = 0; i < n / 2; i++)
    printf ("%6.3f  ", solution[n + i]);
  printf ("\n");
  for (i = n / 2; i < n; i++)
    printf ("%6.3f  ", solution[n + i]);

  printf ("\nExpected solution is\n");
  for (i = 0; i < n / 2; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\n");
  for (i = n / 2; i < n; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\n");
  for (i = 0; i < n / 2; i++)
    {
      printf ("%6.3f  ", expected_sol[n + i]);
      expected_sol[n + i] -= solution[n + i];
    }
  printf ("\n");
  for (i = n / 2; i < n; i++)
    {
      printf ("%6.3f  ", expected_sol[n + i]);
      expected_sol[n + i] -= solution[n + i];
    }
  printf ("\n");

  i = 1;
  j = n * nRhs;
  euclidean_norm = dnrm2 (&j, expected_sol, &i);

  /*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();

  if (euclidean_norm < 1.0e-12)
    {
      printf ("This example has successfully PASSED through all steps of computation!\n");
      return 0;
    }
  else
    {
      printf ("This example may have FAILED as the computed solution differs\n");
      printf ("much from the expected solution (Euclidean norm is %e).\n", euclidean_norm);
      return 1;
    }
  /*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
failure:printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
  MKL_Free_Buffers ();
  return 1;
}
