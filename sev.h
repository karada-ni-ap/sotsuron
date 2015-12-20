#include <iostream>
#include <Eigen/Dense>
#include "const.h"

using namespace std;
using namespace Eigen;

void initA();
double half_ip(MatrixXd X, VectorXd p);
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L);
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L);