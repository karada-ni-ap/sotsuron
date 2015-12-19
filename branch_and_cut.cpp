#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"
#include "local_and_relax.h"
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

void classQ::makeQ(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	Ux = ux;
	Lx = lx;
	Uy = uy;
	Ly = ly;

	Q_U = relaxation(ux, lx, uy, ly);
}

pair<classQ, classQ> classQ::devide(){
	classQ Q1;
	classQ Q2;

	int longest_x = 0;
	int longest_y = 0;
	double max_x = -1;
	double max_y = -1;

	for (int i = 0; i < d; i++){
		if (Ux(i) - Lx(i) > max_x){
			max_x = Ux(i) - Lx(i);
			longest_x = i;
		}
	}

	for (int j = 0; j < m; j++){
		if (Uy(j) - Ly(j) > max_y){
			max_y = Uy(j) - Ly(j);
			longest_y = j;
		}
	}

	VectorXd ex = VectorXd::Zero(d);
	VectorXd ey = VectorXd::Zero(m);

	ex(longest_x) = 1;
	ey(longest_y) = 1;

	max_x /= 2.0;
	max_y /= 2.0;

	Q1.makeQ(Ux, Lx + max_x*ex, Uy, Ly + max_y*ey);
	Q2.makeQ(Ux - max_x*ex, Lx, Uy - max_y*ey, Ly);

	return make_pair(Q1, Q2);
}

void classQ::calculate_lo(){
	Q_L = local_opt(Ux, Lx, Uy, Ly);
}

double branch_and_cut(){
	return 0;
}