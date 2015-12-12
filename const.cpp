#include <Eigen/Dense>
using namespace Eigen;

extern		 int t=0;
extern const int T=30;

extern const int select = 3;

extern const int d=2;
extern const int m=2;
extern const int n=2;

extern const double Ub = 3.0;
extern const double Lb = -3.0;

extern const double Alp0 = 0.1;
extern const double rho = 0.8;

extern const double eps_bfgs = 1.0e-12;
extern const int    ite_bfgs = 100;
extern const int	num_bfgs = 20;

extern const int    mean=0;
extern const double sigma_thre = 1.0e-14;

extern double maxf = -100;

extern VectorXd x_next = VectorXd::Zero(d);
extern VectorXd x_opt  = VectorXd::Zero(d);

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern const bool bfgs_or_rand = true; // trueÇ»ÇÁbfgsÇÅCfalseÇ»ÇÁbound_randópÇ¢ÇÈ