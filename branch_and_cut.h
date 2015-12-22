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

	classQ(){
		next = NULL;
		prev = NULL;
	}

	void makeQ(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly);
	pair<classQ, classQ> devide();
	void calculate_lo();
};

class Qlist{
public:
	classQ root;

	Qlist(){
		root.next = &root;
		root.prev = &root;
		root.Q_L  = -Inf;
	}

	void add(classQ Q);
	classQ extract();
	void delete_tail(double L);
	double maxL();
};

double branch_and_cut();