#include <complex.h>

#ifndef _FMM_DATA_STRUCTURE_
#define _FMM_DATA_STRUCTURE_

typedef double complex dcomplex;

typedef struct fmmbox {
  double center[3];
  int level, boxid, parent, nchild, child[8], npts, addr;
} fmmbox;

typedef struct fmmlist {
  int *list1, *list3, *list4, *colleague; 
} fmmlist;

typedef struct listnode {
  struct listnode *next;
  int data;
} listnode;
#endif

