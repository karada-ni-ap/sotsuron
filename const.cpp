#include <Eigen/Dense>
using namespace Eigen;

extern const int Inf = 10000;

extern		 int t=0;
extern const int T=50;

extern const int d=3;
extern const int m=5;
extern const int n=5;

extern const VectorXd Ux0 = 2.0	*VectorXd::Constant(d, 1.0);
extern const VectorXd Lx0 = -2.0*VectorXd::Constant(d, 1.0);
extern const VectorXd Uy0 = 3.0	*VectorXd::Constant(m, 1.0);
extern const VectorXd Ly0 = -3.0*VectorXd::Constant(m, 1.0);

extern const double Alp0 = 0.3;
extern const double rho = 0.8;
extern const double c1_bfgs = 0.5;

extern const double eps_bfgs = 1.0e-12;
extern const int    ite_bfgs = 100;
extern const int	num_bfgs = 50;

extern const int    mean=0;
extern const double sigma_thre = 1.0e-14;

extern double maxf = -Inf;

extern VectorXd x_next = VectorXd::Zero(d);
extern VectorXd x_opt  = VectorXd::Zero(d);

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern MatrixXd** A = new MatrixXd*[d + 1];
extern MatrixXd*  B = new MatrixXd [m + 1];
extern MatrixXd*  C = new MatrixXd [d + 1];

extern const int    ite_subgrad = 169;
extern const double Alp_subgrad = 1.0;
extern const double eps_subgrad = 1.0e-8;

extern const int    ite_local = 40;
extern const double Alp_local = 0.5;
extern const double eps_local = 1.0e-8;

extern const int    ite_relax = 200;
extern const double Alp_relax = 0.5;
extern const double eps_relax = 1.0e-8;

extern const int    ite_bc = 50;
extern const double eps_bc = 1.0e-3;

extern const int select = 0;			// select = 0ÇÃÇ∆Ç´sevÇÃílÇï‘Ç∑
extern const bool bfgs_or_rand = true;	// trueÇ»ÇÁbfgsÇÅCfalseÇ»ÇÁbound_randópÇ¢ÇÈ