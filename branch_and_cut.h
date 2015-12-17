#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"

using namespace std;
using namespace Eigen;

double local_opt(VectorXd Ux0, VectorXd Lx0, VectorXd Uy0, VectorXd Ly0);
double relaxation(VectorXd Ux0, VectorXd Lx0, VectorXd Uy0, VectorXd Ly0);