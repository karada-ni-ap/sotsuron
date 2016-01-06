#include "const.h"

extern const int Inf = 100000;

extern		 int t=0;
extern const int T=40;

extern int d=5;
extern int m=10;
extern int n=15;

extern VectorXd Ux0 = 2.0	*VectorXd::Constant(d, 1.0);
extern VectorXd Lx0 = -2.0	*VectorXd::Constant(d, 1.0);
extern VectorXd Uy0 = 2.0	*VectorXd::Constant(m, 1.0);
extern VectorXd Ly0 = -2.0	*VectorXd::Constant(m, 1.0);

extern const double alp0 = 1.0;
extern const double eps_alp = 1.0e-10;
extern const double rho = 0.9;
extern const double c1 = 0.6;
extern const double c2 = 0.8;

extern const double eps_y = 1.0e-3;

extern const int    ite_bfgs = 100000;
extern const double eps_bfgs = 1.0e-2;

extern const int    ite_sdm = 100000;
extern const double eps_sdm = 1.0e-2;

extern const int num_of_start = 40;

extern       double mean=0;
extern const double sigma_thre = 1.0e-10;
extern const double xi = 2.0;

extern double maxf_BO   = -Inf;
extern double maxf_lsBO = -Inf;

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern MatrixXd** A = new MatrixXd*[d + 1];
extern MatrixXd*  B = new MatrixXd [m + 1];
extern MatrixXd*  C = new MatrixXd [d + 1];

extern const double beta = 8.0;

extern const int    ite_sev = 100000;
extern const double eps_sev = 1.0e-2;

extern const int    ite_local = 100000;
extern const double eps_local = 1.0e-2;

extern const int    ite_relax = 100000;
extern const double eps_relax = 1.0e-2;

extern const int    ite_bc = 20;
extern const double eps_bc = 1.0e-2;

extern int t_find = 0;
extern int l_find = 0;
extern int k_find = 0;

extern clock_t finding_time_BO   = 0;
extern clock_t finding_time_BC   = 0;
extern clock_t finding_time_lsBO = 0;

extern const int  select = 0;			// select = 0のときsevの値を返す
extern       bool BO_or_lsBO = true;	// アルゴリズムの切り替え
extern       bool SDMorBFGS = true;		// argmaxの切り替え（SDM⇔BFGS）