#include <complex.h>
#include "fmm_ds.h"

#define MAX(a,b) ( (a) < (b) ? (b) : (a) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define NINT(a) ( (a) >= 0 ? (int)((a)+0.5) : (int)((-a)+0.5)*(-1) )
#define TRUE 1
#define FALSE 0

// Part I: Traversal and operator functions

//void ProcFarField(void);
void ProcFarField(int ifpot, int ifgrad, int ifhess);
void ProcNearField(void);

void Aggregate(const int ibox);

//void Disaggregate(const int ibox);
void Disaggregate (int ibox,  int ifpot, int ifgrad, int ifhess);

void TSM(const int ibox);

void TMM(const int ibox);

void TME(const int ibox);

void LIST4(const int ibox);

void TSL(const int ibox, const int member);

//void TLT(const int ibox);
void TLT(const int ibox, int ifpot, int ifgrad, int ifhess);

void TEL(const int ibox);

void TLL(const int ibox);

inline void DIRECT_D(const int ibox1, const int ibox2, double *fmmpot, double *fmmfield,	
       double *fmmdxx, double *fmmdyy, double *fmmdzz, double *fmmdxy, double *fmmdxz, double *fmmdyz);

inline void DIRECT_S(const int ibox, double *fmmpot, double *fmmfield,
       double *fmmdxx, double *fmmdyy, double *fmmdzz, double *fmmdxy, double *fmmdxz, double *fmmdyz); 

void LIST3(const int ibox);

void LIST1(const int ibox);

void TMT(const int ibox, const int pid, const int member);

void kernel(const double *T, const double *S, const double charge, 
	    double *pot, double *field1, double *field2, double *field3, 
	    double *DXX, double *DYY, double *DZZ, double *DXY, double *DXZ,double *DYZ); 

// Part II: Construct graph functions
int create_all_box(const double *ploc, const int nparts, const int s, int *perm, 
		int *nboxes, int *nlev, double *size,
		fmmbox **boxes, int **content);

void divide_box(double *bcenter, const double *ploc, int *start, 
		const int npoints, int *swap, int *addrs, int *counts);

void assign_particle(const double *ploc, int *start, const int npoints, 
		     const int direction, const double cutoff, int *swap, int *undercutoff);

int create_colleague(const fmmbox *boxes, const int *content, const int nlev, 
		 const double size, fmmlist *list);

int if_touch(const double size, const int levela, const double *centera, 
	     const int levelb, const double *centerb);

int create_list134(const fmmbox *boxes, const int *content, const int nlev, 
		   const int nboxes, const double size, fmmlist *list);

void create_list_up(const int mybox, int *uall, int *nuall, int *xuall, int *yuall, 
		    int *u1234, int *nu1234, int *x1234, int *y1234);

void create_list_down(const int mybox, int *dall, int *ndall, int *xdall, int *ydall, 
		      int *d5678, int *nd5678, int *x5678, int *y5678);

void create_list_north(const int mybox, int *nall, int *nnall, int *xnall, int *ynall,
		       int *n1256, int *nn1256, int *x1256, int *y1256,
		       int *n12, int *nn12, int *x12, int *y12,
		       int *n56, int *nn56, int *x56, int *y56);

void create_list_south(const int mybox, int *sall, int *nsall, int *xsall, int *ysall,
		       int *s3478, int *ns3478, int *x3478, int *y3478,
		       int *s34, int *ns34, int *x34, int *y34,
		       int *s78, int *ns78, int *x78, int *y78);

void create_list_east(const int mybox, int *eall, int *neall, int *xeall, int *yeall,
		      int *e1357, int *ne1357, int *x1357, int *y1357,
		      int *e13, int *ne13, int *x13, int *y13,
		      int *e57, int *ne57, int *x57, int *y57,
		      int *e1, int *ne1, int *x1, int *y1,
		      int *e3, int *ne3, int *x3, int *y3,
		      int *e5, int *ne5, int *x5, int *y5,
		      int *e7, int *ne7, int *x7, int *y7);

void create_list_west(const int mybox, int *wall, int *nwall, int *xwall, int *ywall,
		      int *w2468, int *nw2468, int *x2468, int *y2468,
		      int *w24, int *nw24, int *x24, int *y24,
		      int *w68, int *nw68, int *x68, int *y68,
		      int *w2, int *nw2, int *x2, int *y2,
		      int *w4, int *nw4, int *x4, int *y4,
		      int *w6, int *nw6, int *x6, int *y6,
		      int *w8, int *nw8, int *x8, int *y8);

void free_list(fmmlist *list, const int nboxes);

void vecAdd(const dcomplex *b, dcomplex *a, const int nelements);

