#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	VectorXd x_next = VectorXd::Zero(d);

	for (t=0; t < T; t++){
		x_next = argmax_u();
		
		//データセットの更新
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		maxf = max(f(t), maxf);
		updateK(x_next);
	}

	MatrixXd A = K*Kinv;
	cout << A(T-1, T-2) << endl;
	return 0;
}