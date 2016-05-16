#include <sys/time.h>
#include <stdio.h>
#include "fmm_ds.h"
// global parameter for rpy march 3th wenhua
double ha,c0,c1,c2; 

int PTERMS, NLAMBS, PGSZ, NLEV, NBOXES, *NUMFOUR, *NUMPHYS, 
  *PERM, *CONTENT, *LEAFBOX, NPARTS, PTERMS2, NEXPMAX;

double *FMMLOC, *FMMCHARGE, *FMMPOT, *FMMFIELD, *FMMPOTN, *FMMFIELDN, 
  *FMMDXX,  *FMMDYY,  *FMMDZZ,  *FMMDXY,  *FMMDXZ,  *FMMDYZ, 
  *FMMDXXN, *FMMDYYN, *FMMDZZN, *FMMDXYN, *FMMDXZN, *FMMDYZN, 
  SIZE,  *WHTS, *RLAMS, *RDPLUS, *RDMINUS, *RDSQ3, *RDMSQ3,
  *RLSC, *CARRAY, *ZS;

double *RPYCHARGE, *RPYPOT;

dcomplex *MPOLE, *LOCAL, *LEXPU, *LEXPD, *LEXPN, *LEXPS, *LEXPE, 
  *LEXPW, *XS, *YS, *FEXPE, *FEXPO, *FEXPBACK;

fmmbox *BOXES;

fmmlist *LIST;



