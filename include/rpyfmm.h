/*
  rpyfmm.h
*/
int RPYfmm( const double length, const double beta, const int s,  const int accuracy, 
	const int nparts, const double *ploc, const double *pcharge,  double *RPY );
/*void rpy3sup_eval( const double *T, const double *S, const double *charge, double *pot);

inline void rpy_DIRECT_D(const int ibox1, const int ibox2, double *fmmpot);

inline void DIRECT_S(const int ibox, double *fmmpot);

void rpy3nearfield();

void rpyLIST3(const int ibox);
void rpyLIST1(const int ibox);

void test_init( const double beta, const int s, const int accuracy, const int nparts, 
		double *charge, double *pot, double *field, 
		double *dxx, double *dyy, double *dzz,
	        double *dxy, double *dxz, double *dyz);

void test_data(int nparts, int distribution, double *ploc, double *pcharge);

void test_clean( double *charge, double *pot, double *field,
	        double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz );
void rpydirect(const int nparts, double *ploc, double *pcharge, const int i, double *dpot); */
