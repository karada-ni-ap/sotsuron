#include "const.h"

double pdf(double x);
double cdf(double x);
double kernel(Eigen::VectorXd x1, Eigen::VectorXd x2);
VectorXd bound_rand();
VectorXd projection(VectorXd z, VectorXd U, VectorXd L);