// Shell representation
int  fcc_array(int nc, double bradius, double spacing, double *sphs);
int shell_model(const int nc, const double bradius, const double spacing, const double *ShellRadius, double *ShellSphs);
int Tmat(double *T,double *ShellSphs,int N);
//int Qmatrix(double *Q,double *ShellSphs, int N);
void LCepsilon(double *lc_eps);
double sph_norm(double *v, double radius );
int  rotmat(double *r, double *M, double *lc_eps);
//int KrylovFMM(const double length, const double beta, const int s,
//	      const int accuracy, const int nparts, double *ShellSphs,
//              double *rhs, double *solution);
int GMRES_FMM(const double bradius, const double length, const double beta, const int s, 
	      const int accuracy, const int  nparts, double *ShellSphs, 
	      const double *rhs, double *solution);
void GetRightSide(const double bradius, const double length, const double beta, const int s, const int  accuracy,
	   const int nparts, double *ShellSphs, double *uvec,  const double *T, double *V );
