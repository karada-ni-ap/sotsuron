#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

void debug_inside(VectorXd x_next, VectorXd x_opt){ //tはデータセットのサイズ
		cout << "t is " << t << endl;
		cout << x_next.transpose() << endl;
		
		cout << "-----------------------------------------------" << endl;
}

void debug_last(){
	VectorXd x0 = VectorXd::Zero(d);
	VectorXd y0 = VectorXd::Zero(m);

	MatrixXd Wmax = MatrixXd::Zero(d, m);
	MatrixXd Wmin = MatrixXd::Zero(d, m);

	pair<MatrixXd, MatrixXd> Pair = Wbound(Ux0, Lx0, Uy0, Ly0);

	Wmax = Pair.first;
	Wmin = Pair.second;

	//cout << branch_and_cut() << endl;

	clock_t start = clock();

	BO();
	cout << "t_find : " << t_find << endl;
	cout << "maxf : " << maxf << endl;

	clock_t end = clock();

	cout << "time : " << end - start << endl;

	//cout << Wbound(Ux0, Lx0, Uy0, Ly0).second << endl;
	//cout << relaxation(Ux0, Lx0, Uy0, Ly0) << endl;
}