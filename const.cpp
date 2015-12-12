#include <Eigen/Dense>
using namespace Eigen;

extern		 int t=0;
extern const int T=50;

extern const int d=2;
extern const int m=2;
extern const int n=2;

extern const double Ub = 1.0;
extern const double Lb = -1.0;

extern const double Alp0 = 0.5;
extern const double rho = 0.9;

//extern const double Alp = 0.01;
//extern const double conv_thre = 1.0e-12;
extern const int    ite_bfgs = 10;
extern const int	num_bfgs = 50;

extern const int    mean=0;
extern const double sigma_thre = 1.0e-14;

extern double maxf = -100;

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);

extern const bool bfgs_or_rand = false; // trueÇ»ÇÁbfgsÇÅCfalseÇ»ÇÁbound_randópÇ¢ÇÈ