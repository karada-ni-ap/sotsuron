#include <Eigen/Dense>
using namespace Eigen;

extern		 int t=0;
extern const int T=150;

extern const int d=6;
extern const int m=2;
extern const int n=2;

extern const double ub = 2.0;
extern const double lb = -2.0;

extern const int mean=0;

extern double maxf = -100;

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);