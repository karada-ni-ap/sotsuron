#include "const.h"
#include <Eigen/Dense>

Eigen::VectorXd k(Eigen::VectorXd x);
void updateK(Eigen::VectorXd x);
double mu(Eigen::VectorXd x);
double sigma(Eigen::VectorXd x);
double u(Eigen::VectorXd x);
Eigen::VectorXd argmax_u();