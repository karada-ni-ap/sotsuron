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
	Q_Ux = ux;
	Q_Lx = lx;
	Q_Uy = uy;
	Q_Ly = ly;

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
		if (Q_Ux(i) - Q_Lx(i) > max_x){
			max_x = Q_Ux(i) - Q_Lx(i);
			longest_x = i;
		}
	}

	for (int j = 0; j < m; j++){
		if (Q_Uy(j) - Q_Ly(j) > max_y){
			max_y = Q_Uy(j) - Q_Ly(j);
			longest_y = j;
		}
	}

	VectorXd ex = VectorXd::Zero(d);
	VectorXd ey = VectorXd::Zero(m);

	ex(longest_x) = 1;
	ey(longest_y) = 1;

	max_x /= 2.0;
	max_y /= 2.0;

	Q1.makeQ(Q_Ux, Q_Lx + max_x*ex, Q_Uy, Q_Ly + max_y*ey);
	Q2.makeQ(Q_Ux - max_x*ex, Q_Lx, Q_Uy - max_y*ey, Q_Ly);

	return make_pair(Q1, Q2);
}

void classQ::calculate_lo(){
	Q_L = local_opt(Q_Ux, Q_Lx, Q_Uy, Q_Ly);
}

double branch_and_cut(){
	return 0;
}