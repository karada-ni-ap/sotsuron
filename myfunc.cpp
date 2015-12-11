#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"

using namespace std;
using namespace Eigen;

double pdf(double x){
	return exp(-0.5 * x * x) / sqrt(2 * M_PI);
}

double cdf(double x){
	return 0.5 * ( 1 + erf(x/sqrt(2)) );
}

double kernel(Eigen::VectorXd x1, Eigen::VectorXd x2){
	return exp(-0.5 * (x1-x2).squaredNorm());
}

VectorXd bound_rand(int num){
	return VectorXd::Constant(num, (ub + lb) / 2) + ((ub - lb)/2)* VectorXd::Random(num);
}