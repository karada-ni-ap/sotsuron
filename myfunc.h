#include "const.h"

double pdf(double x);
double cdf(double x);
double kernel(Eigen::VectorXd x1, Eigen::VectorXd x2);
double h_rule(int k, double h_conv);
VectorXd bound_rand();
VectorXd projection(VectorXd z, VectorXd U, VectorXd L);