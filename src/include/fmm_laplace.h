#include <complex.h>

int NEXPTOT, NEXPTOTP, NTHMAX, IFL_UP[8], IFL_DN[8];

double *DC, *YTOPC, *YTOPCS, *YTOPCSINV, *SCALE;

void TME_PHASE1(const dcomplex *mpole, dcomplex *mexpup, dcomplex *mexpdown, const int ibox);

void TME_PHASE2(const dcomplex *mexpf, dcomplex *mexpphys, const int ibox);

void TEL_UD(const int ibox);

void TEL_NS(const int ibox);

void TEL_EW(const int ibox);

void TEL_PHASE1(dcomplex *mexpf, const dcomplex *mexpphys, const int level);

void TEL_PHASE2(dcomplex *local, const int iexpu, const dcomplex *mexpup, 
		const int iexpd, const dcomplex *mexpdown, const int level);

void ROTZ2Y(const dcomplex *mpole, dcomplex *mrotate);

void ROTY2Z(const dcomplex *mpole, dcomplex *mrotate);

void ROTZ2X(const dcomplex *mpole, dcomplex *mrotate, const int direction);

void vwts(const int nlambs, double *rlams, double *whts);

void numthetahalf(const int nlambs, int *numtets);

void numthetafour(const int nlambs, int *numtets);

void rlscini(const int nlambs, const int pterms, const double *rlams, double *rlsc);

void mkfexp(const int nlambs, const int *numfour, const int *numphys,
	    dcomplex *fexpe, dcomplex *fexpo, dcomplex *fexpback);

void mkexps(const double *rlams, const int nlambs, const int *numphys,
	    const int nexptotp, dcomplex *xs, dcomplex *ys, double *zs);

void lgndr(const int nmax, const double x, double *y);

void frmini (double *c, double *cs, double *csinv);

void rotgen(const int pterms, const double *carray, double *rdpi2, double *rdmpi2,
	    double *rdsq3, double *rdmsq3, double *dc);

void bnlcft(double *c, const int pterms);

void fstrtn(const int pterms, double *d, const double *sqc, const double theta);
double Alpha(int l, int m);
double Beta(int l, int m, int option);
