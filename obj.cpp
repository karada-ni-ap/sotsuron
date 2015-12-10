#include <iostream>
#include <Eigen/Dense>
#include "const.h"

using namespace std;
using namespace Eigen;

double obj(VectorXd x){
	VectorXd center = VectorXd::Constant(d,0.5);
	return (x - center).norm() * (-1);
}