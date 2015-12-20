#include "const.h"
#include <Eigen/Dense>

extern VectorXd m1;

Eigen::VectorXd k(Eigen::VectorXd x);
void update_K(Eigen::VectorXd x);

double mu(Eigen::VectorXd x);
double sigma(Eigen::VectorXd x);
double u(Eigen::VectorXd x);

VectorXd u_over_k(VectorXd x);
MatrixXd k_over_x(VectorXd x);
VectorXd u_over_x(VectorXd x);
double BO();