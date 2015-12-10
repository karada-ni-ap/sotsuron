#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"

using namespace std;
using namespace Eigen;

double pdf(Eigen::VectorXd x);
double cdf(double x);
double kernel(Eigen::VectorXd x1, Eigen::VectorXd x2);