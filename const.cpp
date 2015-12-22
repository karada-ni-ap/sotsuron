#include <Eigen/Dense>
using namespace Eigen;

extern const int Inf = 100000;

extern		 int t=0;
extern       int t_find=0;
extern const int T=10;

extern const int d=2;
extern const int m=2;
extern const int n=2;

extern const VectorXd Ux0 = 3.0		*VectorXd::Constant(d, 1.0);
extern const VectorXd Lx0 = -3.0	*VectorXd::Constant(d, 1.0);
extern const VectorXd Uy0 = 2.0		*VectorXd::Constant(m, 1.0);
extern const VectorXd Ly0 = -2.0	*VectorXd::Constant(m, 1.0);

extern const double Alp0 = 0.3;
extern const double rho = 0.8;
extern const double c1_bfgs = 0.5;

extern const double eps_bfgs = 1.0e-12;
extern const int    ite_bfgs = 50;
extern const int	num_bfgs = 36;

extern const int    mean=0;
extern const double sigma_thre = 1.0e-14;

extern double maxf = -Inf;

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern MatrixXd** A = new MatrixXd*[d + 1];
extern MatrixXd*  B = new MatrixXd [m + 1];
extern MatrixXd*  C = new MatrixXd [d + 1];

extern const double beta = 8.0;

extern const int    ite_sev = 10000;
extern const double eps_sev = 1.0e-4;

extern const int    ite_local = 10000;
extern const double eps_local = 1.0e-4;

extern const int    ite_relax = 100;
extern const double eps_relax = 1.0e-6;
extern const double alpW0 = 1.0;

extern const int    ite_bc = 10;
extern const double eps_bc = 1.0e-4;

extern const int  select = 0;			// select = 0ÇÃÇ∆Ç´sevÇÃílÇï‘Ç∑
extern const bool bfgs_or_rand = true;	// trueÇ»ÇÁbfgsÇÅCfalseÇ»ÇÁbound_randópÇ¢ÇÈ