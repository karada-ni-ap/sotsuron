#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "myfunc.h"
#include "argmax.h"

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

	for (int j = 0; j < (m + 1); j++)
		B[j] = MatrixXd::Zero(n, n);

	for (int i = 0; i < (d + 1); i++)
		C[i] = MatrixXd::Zero(n, n);
}

double half_ip(MatrixXd X, VectorXd p){
	int size = p.size();
	
	double val = 0;
	for (int i = 0; i < size; i++){
		for (int j = i; j < size; j++){
			val += X(i, j)*p(i)*p(j);
		}
	}
	return val;
}

pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L, double Alpha0){
	for (int j = 0; j< (m + 1); j++){
		B[j] = A[0][j];
		for (int i = 1; i <= d; i++){
				B[j] += (x(i-1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);
	VectorXd y   = y0; //劣勾配法の初期点
	VectorXd dir = VectorXd::Zero(m);
	double   min = 0;
	double   pre = Inf;
	double   Alp = Alpha0;

	//yについて劣勾配法
	for (int k = 0; k < ite_subgrad; k++){
		//yの更新
		y += Alp*dir;
		y = projection(y, U, L);

		//Axyの構成
		Axy = B[0];
		for (int j = 1; j <= m; j++)
			Axy += y(j-1)*B[j];


		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		
		pn  = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		for (int j = 1; j <= m; j++){
			dir(j-1) = half_ip(B[j], pn);
		}

		if (dir.norm() < eps_subgrad			//収束判定
			|| abs(pre - min) < eps_subgrad		//ほんとは良くない
		)
			break;

		Alp = Alpha0 / sqrt(k+1);	//ステップサイズの更新
		pre = min;
	}

	return make_pair(min, y);
}

pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L, double Alpha0){
	for (int i = 0; i< (d + 1); i++){
		C[i] = A[i][0];
		for (int j = 1; j <= m; j++){
			C[i] += (y(j - 1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);
	VectorXd x   = x0; //劣勾配法の初期点
	VectorXd dir = VectorXd::Zero(d);
	double   min = 0;
	double   pre = Inf;
	double   Alp = Alpha0;

	//xについて劣勾配法
	for (int k = 0; k < ite_subgrad; k++){
		//xの更新
		x += Alp*dir;
		x = projection(x, U, L);

		//Axyの構成
		Axy = C[0];
		for (int i = 1; i <= d; i++)
			Axy += x(i - 1)*C[i];


		SelfAdjointEigenSolver<MatrixXd> es(Axy);

		pn = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		for (int i = 1; i <= d; i++){
			dir(i - 1) = half_ip(C[i], pn);
		}

		if (dir.norm() < eps_subgrad			//収束判定
			|| abs(pre - min) < eps_subgrad		//ほんとはよくない
		)
			break;

		Alp = Alpha0 / sqrt(k + 1);	//ステップサイズの更新
	}

	return make_pair(min, x);
}