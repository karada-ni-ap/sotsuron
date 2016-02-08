#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "argmax.h"

using namespace std;
using namespace Eigen;

void initA(){
	A = new MatrixXd*[d + 1];
	B = new MatrixXd[m + 1];
	C = new MatrixXd[d + 1];

	for (int i = 0; i < (d+1); i++)
		A[i] = new MatrixXd[m+1];

	for (int i = 0; i < (d + 1); i++){
		for (int j = 0; j < (m + 1); j++){
			/*
			A[i][j] = MatrixXd::Zero(n, n);
			
			VectorXd tmp = VectorXd::Random(n);
			int evenodd = (i + j) % 2;

			for (int k = 0; k < n; k++){
				if (evenodd == 0)
					A[i][j](k, k) = abs(tmp(k));
				else
					A[i][j](k, k) = -abs(tmp(k));
			}
			*/
			
			
			MatrixXd tmp = MatrixXd::Random(n, n);
			A[i][j] = 0.5*(tmp + tmp.transpose()); // これだと一様分布にならない
			

			/*
			A[i][j] = MatrixXd::Random(n, n);
			for (int k = 1; k <= n - 1; k++)
				for (int l = 0; l <= n - 2; l++)
					A[i][j](k, l) = A[i][j](l, k);
			*/
		}
	}

	for (int j = 0; j < (m + 1); j++)
		B[j] = MatrixXd::Zero(n, n);

	for (int i = 0; i < (d + 1); i++)
		C[i] = MatrixXd::Zero(n, n);
}

pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L){
	for (int j = 0; j< (m + 1); j++){
		B[j] = A[0][j];
		for (int i = 1; i <= d; i++){
				B[j] += (x(i-1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	VectorXd y     = y0; //劣勾配法の初期点
	VectorXd y_opt = VectorXd::Zero(m);

	VectorXd dir = VectorXd::Zero(m);

	double min     = 0;
	double opt = -Inf;

	double subgrad_norm = 0;

	double h = 0.0;
	double h_conv = 0.01;

	double alp = 0.0;

	int k_update = 0;

	//yについて劣勾配法
	for (int k = 0; k < ite_sev; k++){
		//Axyの構成
		Axy = B[0];
		for (int j = 1; j <= m; j++)
			Axy += y(j-1)*B[j];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		
		pn  = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		//cout << k << " 's sev_x : " << min << endl;

		if (min > opt){
			//cout << "updating!" << endl;

			y_opt = y;
			k_update = k;
			opt = min;
		}

		for (int j = 1; j <= m; j++){
			dir(j - 1) = pn.transpose() * B[j] * pn;
		}

		//ステップサイズ
		h = h_rule(k, h_conv);
		subgrad_norm = dir.norm();
		alp = h / subgrad_norm;

		//yの更新
		y = projection(y + alp*dir, U, L);

		//停止条件の判定
		if (k - k_update > N_sc){
			//cout << "sev_x 's k : " << k << endl;
			break;
		}

	}

	//cout << "sev_x 's opt : " << opt << endl;

	return make_pair(opt, y_opt);
}

pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L){
	for (int i = 0; i < (d + 1); i++){
		C[i] = A[i][0];
		for (int j = 1; j <= m; j++){
			C[i] += (y(j - 1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	VectorXd x     = x0; //劣勾配法の初期点
	VectorXd x_pre = x0;
	VectorXd x_opt = VectorXd::Zero(d);

	VectorXd dir = VectorXd::Zero(d);

	double min     = 0;
	double min_pre = 0;

	double opt = -Inf;

	double h = 0.0;
	double h_conv = 0.01;

	double alp = 0.0;

	double subgrad_norm = 0.0;

	int k_update = 0;

	//xについて劣勾配法
	for (int k = 0; k < ite_sev; k++){
		//Axyの構成
		Axy = C[0];
		for (int i = 1; i <= d; i++)
			Axy += x(i - 1)*C[i];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);

		pn = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		//cout << k << " 's sev_y : " << min << endl;

		if (min > opt){
			//cout << "updating!" << endl;

			x_opt = x;
			k_update = k;
			opt = min;
		}

		for (int i = 1; i <= d; i++){
			dir(i - 1) = pn.transpose() * C[i] * pn;
		}

		//ステップサイズ
		h = h_rule(k, h_conv);
		subgrad_norm = dir.norm();
		alp = h / subgrad_norm;

		//xの更新
		x = projection(x + alp*dir, U, L);

		if (k - k_update > N_sc){
			//cout << "sev_y 's k : " << k << endl;
			break;
		}

	}

	return make_pair(opt, x_opt);
}

pair<double, VectorXd> local_search(VectorXd x0, VectorXd y0){
	VectorXd x = x0;
	VectorXd y = y0;

	pair<double, VectorXd> X;
	pair<double, VectorXd> Y;

	double tmp = 0;
	double pre = 0;
	//double max = -Inf;

	for (int k = 0; k < ite_local; k++){
		Y = sev_x(x, y, Uy0, Ly0);
		y = Y.second;
		
		X = sev_y(y, x, Ux0, Lx0);
		x = X.second;

		pre = tmp;
		tmp = X.first;

		//max = tmp>max ? tmp : max;

		//収束判定
		if (abs(tmp - pre) < eps_local)
			break;
	}

	return make_pair(tmp, x); //局所最適解(x*, y*)のx*を返す
}