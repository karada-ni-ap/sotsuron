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
		//cout << "t is " << t << endl;
		cout << x_next.transpose() << endl;
		
		//cout << "-----------------------------------------------" << endl;
}

void debug_last(){
	;
}