#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "BO.h"

using namespace std;
using namespace Eigen;

bool judge(VectorXd x);
MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y);
double back_track(VectorXd X, VectorXd Grad, VectorXd dir);
VectorXd bfgs(VectorXd x0);
VectorXd argmax_u();