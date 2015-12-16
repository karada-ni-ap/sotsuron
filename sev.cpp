#include <iostream>
#include <Eigen/Dense>
#include "obj.h"
#include "const.h"
#include "myfunc.h"

using namespace std;
using namespace Eigen;

void initA(){
	for (int i = 0; i < (d+1); i++)
		A[i] = new MatrixXd[m+1];

	for (int i = 0; i < (d + 1); i++){
		for (int j = 0; j < (m + 1); j++){
			MatrixXd tmp = MatrixXd::Random(n, n);
			A[i][j] = 0.5*(tmp + tmp.transpose());
		}
	}

	for (int j = 0; j< (m + 1); j++)
		B[j] = MatrixXd::Zero(n, n);
}

double half_ip(MatrixXd X, VectorXd p){
	double val = 0;
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			val += X(i, j)*p(i)*p(j);
		}
	}
	return val;
}

double sev(VectorXd x){
	return obj(x,select);

	for (int j = 0; j< (m + 1); j++){
		B[j] = A[0][j];
		for (int i = 1; i <= d; i++){
				B[j] += (x(i-1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);
	VectorXd y   = VectorXd::Zero(m);
	VectorXd dir = VectorXd::Zero(m);
	double   min = 0;
	double   Alp = Alp_subgrad;

	//最急降下法
	for (int k = 0; k < ite_subgrad; k++){
		Axy = B[0];
		for (int j = 1; j <= m; j++)
			Axy += y(j-1)*B[j];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		
		pn  = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		cout << min << endl;

		for (int j = 1; j <= m; j++){
			dir(j-1) = half_ip(B[j], pn);
		}
		
		y += Alp*dir;
		Alp *= rho_subgrad;	
	}

	cout << "------------------------------------" << endl;
	//return min;
}