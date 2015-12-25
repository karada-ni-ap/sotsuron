#include "const.h"
#include "myfunc.h"
#include "BO.h"

MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y);
double back_track(VectorXd X, VectorXd Grad, VectorXd dir);
VectorXd bfgs(VectorXd x0);
VectorXd sdm(VectorXd x0);
VectorXd argmax_u();