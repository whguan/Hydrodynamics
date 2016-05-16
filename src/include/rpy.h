/*
  rpyfmm.h
*/
void rpy3sup_eval(const double *T, const double *S, const double *charge, double *pot);

void rpy_DIRECT_D(const int ibox1, const int ibox2, double *fmmpot);

void rpy_DIRECT_S(const int ibox, double *fmmpot);

void rpy3nearfield();

void rpyLIST3(const int ibox);
void rpyLIST1(const int ibox);

void test_init( const double length, const double ha,const double beta, const int s, const int accuracy, const int nparts, 
		double **charge, double **pot, double **field, 
		double **dxx, double **dyy, double **dzz,
	        double **dxy, double **dxz, double **dyz);

void test_data(double length, int nparts, int distribution, double *ploc, double *pcharge);

void test_clean( double *charge, double *pot, double *field,
	        double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz );
void rpydirect(const int nparts, double *ploc, double *pcharge, const int i, double *dpot);
