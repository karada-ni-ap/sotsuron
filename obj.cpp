#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "sev.h"

using namespace std;
using namespace Eigen;

double obj(VectorXd x, int select){
	double val = 0;

	switch (select)
	{
	case 0:
		return sev_x(x, VectorXd::Zero(m), Uy0, Ly0, Alp_subgrad).first;

	case 1:
		return -sqrt(x.squaredNorm());

	case 2:
		return x(0) * abs(sin(M_PI*x.norm()));

	case 3:
		for (int i = 0; i < d; i++){
			val += x(i)*x(i) - 10 * cos(0.5*M_PI*x(i));
		}
		return val + 10 * d;

	case 4:
		return 0.3 * x.norm() - sqrt(x.squaredNorm());

	default:
		return 0;
	}
}
