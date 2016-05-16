/*
  adap_fmm.h
*/
void adap_fmm_init(const int accuracy, const int nparts);
//2015
//void adap_fmm_graph(const int nparts, const int s, const double beta, 
//		    const double *ploc, const double *pcharge);
//void adap_fmm_compute(void);

//March 3
void adap_fmm_graph(const int nparts, const int s, const double *ploc);
void adap_fmm_charge(const int ifdirect, const double *charge, const double *rpycharge);
void adap_fmm_compute(int ifdirect, int ifpot, int ifgrad, int ifhess);

void adap_fmm_post(const int ifdirect, double *RPY, const int ifpot, double *pot, 
		const int ifgrad, double *field, const int ifhess, double *dxx,
		double *dyy, double *dzz, double *dxy, double *dxz, double *dyz);

void adap_fmm_clean(void);

