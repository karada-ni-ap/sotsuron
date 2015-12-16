#include <Eigen/Dense>
using namespace Eigen;

extern		 int t=0;
extern const int T=50;

extern const int d=2;
extern const int m=2;
extern const int n=2;

extern const double Ub = 2.0;
extern const double Lb = -2.0;

extern const double Alp0 = 0.2;
extern const double rho = 0.9;
extern const double c1_bfgs = 0.9;

extern const double eps_bfgs = 1.0e-12;
extern const int    ite_bfgs = 50;
extern const int	num_bfgs = 30;

extern const int    mean=0;
extern const double sigma_thre = 1.0e-14;

extern double maxf = -100;

extern VectorXd x_next = VectorXd::Zero(d);
extern VectorXd x_opt  = VectorXd::Zero(d);

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern MatrixXd** A = new MatrixXd*[d + 1];
extern MatrixXd*  B = new MatrixXd [m + 1];

extern const int    ite_subgrad = 30;
extern const int    Alp_subgrad = 1.0;
extern const double rho_subgrad = 0.9;

extern const int select = 3;
extern const bool bfgs_or_rand = true; // trueÇ»ÇÁbfgsÇÅCfalseÇ»ÇÁbound_randópÇ¢ÇÈ