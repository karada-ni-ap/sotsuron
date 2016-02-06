#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "argmax.h"

using namespace std;
using namespace Eigen;

void initA_byGoh(){
	d = 1;
	m = 1;
	n = 3;

	Ux0 = 7.0	*VectorXd::Constant(d, 1.0);
	Lx0 = -1.0	*VectorXd::Constant(d, 1.0);

	Uy0 = 7.0	*VectorXd::Constant(m, 1.0);
	Ly0 = -3.0	*VectorXd::Constant(m, 1.0);

	A = new MatrixXd*[d + 1];
	B = new MatrixXd[m + 1];
	C = new MatrixXd[d + 1];

	for (int i = 0; i < (d + 1); i++)
		A[i] = new MatrixXd[m + 1];

	for (int i = 0; i < (d + 1); i++){
		for (int j = 0; j < (m + 1); j++){
			A[i][j] = MatrixXd::Zero(3,3);
		}
	}

	A[0][0] << -10, -0.5, -2, -0.5, 4.5, 0, -2, 0, 0;
	A[1][0] << -1.8, -0.1, -0.4, -0.1, 1.2, -1, -0.4, -1, 0;
	A[0][1] << 9, 0.5, 0, 0.5, 0, -3, 0, -3, -1;
	A[1][1] << 0, 0, 2, 0, -5.5, 3, 2, 3, 0;

	for (int j = 0; j < (m + 1); j++)
		B[j] = MatrixXd::Zero(n, n);

	for (int i = 0; i < (d + 1); i++)
		C[i] = MatrixXd::Zero(n, n);
}