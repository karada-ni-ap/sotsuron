#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"
#include "local_and_relax.h"

using namespace std;
using namespace Eigen;

class classQ{
public:
	VectorXd Ux;
	VectorXd Lx;
	VectorXd Uy;
	VectorXd Ly;

	double Q_U;
	double Q_L;

	classQ* next;
	classQ* prev;

	void makeQ(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly);
	pair<classQ, classQ> devide();
	void calculate_lo();
};

double branch_and_cut();