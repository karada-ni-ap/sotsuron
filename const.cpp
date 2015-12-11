#include <Eigen/Dense>
using namespace Eigen;

extern		 int t=0;
extern const int T=50;

extern const int d=2;
extern const int m=2;
extern const int n=2;

extern const double ub = 1.0;
extern const double lb = -1.0;

extern const int mean=0;

extern double maxf = -100;

extern MatrixXd D_q = MatrixXd::Zero(d,T);
extern VectorXd f   = VectorXd::Zero(T);

extern MatrixXd K    = MatrixXd::Zero(T,T);
extern MatrixXd Kinv = MatrixXd::Zero(T,T);