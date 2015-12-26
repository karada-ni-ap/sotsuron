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
	//clock_t start_BC = clock();
	//double max_of_BC = branch_and_cut();
	//clock_t end_BC = clock();

	clock_t start_BO = clock();
	BO();
	clock_t end_BO = clock();

	//初期化
	t = 0;
	D_q = MatrixXd::Zero(d, T);
	f = VectorXd::Zero(T);
	K = MatrixXd::Zero(T, T);
	Kinv = MatrixXd::Zero(T, T);

	clock_t start_lsBO = clock();
	lsBO();
	clock_t end_lsBO = clock();

	cout << "t_find     : " << t_find << endl;
	cout << "max of BO  : " << maxf_BO << endl;
	cout << "time of BO : " << end_BO - start_BO << endl;
	cout << "find timeBO: " << finding_time_BO - start_BO << endl;

	cout << "l_find        : " << l_find << endl;
	cout << "max of lsBO   : " << maxf_lsBO << endl;
	cout << "time of lsBO  : " << end_lsBO - start_lsBO << endl;
	cout << "find time lsBO: " << finding_time_lsBO - start_lsBO << endl;

	//cout << "k_find     : " << k_find << endl;
	//cout << "max of BC  : " << max_of_BC << endl;
	//cout << "time of BC : " << end_BC - start_BC << endl;
	//cout << "find timeBC: " << finding_time_BC - start_BC << endl;
}