#include <iostream>
#include <Eigen/Dense>
#include "obj.h"
#include "const.h"

using namespace std;
using namespace Eigen;

double sev(VectorXd x){
	return obj(x,select);
}