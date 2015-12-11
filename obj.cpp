#include <iostream>
#include <Eigen/Dense>
#include "const.h"

using namespace std;
using namespace Eigen;

double obj(VectorXd x){
	return -sqrt(x.squaredNorm());
}