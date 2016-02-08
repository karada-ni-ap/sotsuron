#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "sev.h"
#include "myfunc.h"

using namespace std;
using namespace Eigen;

double local_opt(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	VectorXd x = (ux + lx) / 2;
	VectorXd y = (uy + ly) / 2;

	pair<double, VectorXd> X;
	pair<double, VectorXd> Y;

	double tmp = 0;
	double pre = 0;
	double max = -Inf;
	
	//cout << "ls starts." << endl;

	for (int k = 0; k < ite_local; k++){
		X = sev_y(y, x, ux, lx);
		x = X.second;

		Y = sev_x(x, y, uy, ly);
		y = Y.second;

		pre = tmp;
		tmp = Y.first;

		max = tmp>max ? tmp : max;

		//cout << tmp << endl;

		//収束判定
		if (abs(tmp - pre) < eps_local)
			break;
	}

	//cout << "ls ends." << endl;

	return max;
}

pair<MatrixXd, MatrixXd> Wbound(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	MatrixXd W[4];

	MatrixXd Wmax = MatrixXd::Constant(d, m, -Inf);
	MatrixXd Wmin = MatrixXd::Constant(d, m,  Inf);

	W[0] = ux * uy.transpose();
	W[1] = ux * ly.transpose();
	W[2] = lx * uy.transpose();
	W[3] = lx * ly.transpose();

	for (int i = 0; i < d; i++){
		for (int j = 0; j < m; j++){
			for (int l = 0; l < 4; l++){
				     if (W[l](i, j) > Wmax(i, j)) Wmax(i, j) = W[l](i, j);
				else if (W[l](i, j) < Wmin(i, j)) Wmin(i, j) = W[l](i, j);
			}
		}
	}

	return make_pair(Wmax, Wmin);
}

double norm(VectorXd x, VectorXd y, MatrixXd W){
	double val = 0;

	for (int i = 0; i < d; i++){
		val += x(i)*x(i);
	}

	for (int j = 0; j < m; j++){
		val += y(j)*y(j);
		for (int i = 0; i < d; i++){
			val += W(i, j)*W(i, j);
		}
	}

	return sqrt(val);
}

MatrixXd projectionW(MatrixXd W, MatrixXd Wmax, MatrixXd Wmin){
	MatrixXd Wret = W;
	for (int i = 0; i < d; i++){
		for (int j = 0; j < m; j++){
			     if (W(i, j) > Wmax(i, j)) Wret(i, j) = Wmax(i, j);
			else if (W(i, j) < Wmin(i, j)) Wret(i, j) = Wmin(i, j);
		}
	}
	return Wret;
}

double relaxation(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	VectorXd dx = VectorXd::Zero(d);
	VectorXd dy = VectorXd::Zero(m);
	MatrixXd dW = MatrixXd::Zero(d, m);
	
	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	pair<MatrixXd, MatrixXd> Pair = Wbound(ux, lx, uy, ly);

	VectorXd x = (ux + lx) / 2;
	VectorXd y = (uy + ly) / 2;
	MatrixXd W = (Pair.first + Pair.second) / 2;

	double min = -Inf;
	double opt = -Inf;

	double subgrad_norm = 0;

	double alp = 0.0;

	double h = 0.0;
	double h_conv = 1.0;

	int k_update = 0;

	for (int k = 0; k < ite_relax; k++){

		//Axyの構成
		Axy = A[0][0];

		for (int i = 1; i <= d; i++){
			Axy += x(i - 1)*A[i][0];
		}
		for (int j = 1; j <= m; j++){
			Axy += y(j - 1)*A[0][j];
			for (int i = 1; i <= d; i++){
				Axy += W(i - 1, j - 1)*A[i][j];
			}
		}

		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		pn = es.eigenvectors().col(0);	//Axyの最小固有ベクトル
		min = es.eigenvalues()(0);		//Axyの最小固有値

		//cout << k << " 's relaxation : " << min << endl;

		if (min > opt){
			opt = min;
			k_update = k;
		}

		//上昇方向
		for (int i = 1; i <= d; i++){
			dx(i-1) = pn.transpose() * A[i][0] * pn;
		}
		for (int j = 1; j <= m; j++){
			dy(j - 1) = pn.transpose() * A[0][j] * pn;
			for (int i = 1; i <= d; i++){
				dW(i-1,j-1)= pn.transpose() * A[i][j] * pn;
			}
		}

		//ステップサイズ
		h = h_rule(k, h_conv);
		subgrad_norm = norm(dx, dy, dW);
		alp = h / subgrad_norm;


		//x,y,Wの更新
		x = projection(x + alp*dx, ux, lx);
		y = projection(y + alp*dy, uy, ly);
		W = projectionW(W + alp*dW, Pair.first, Pair.second);	

		//停止条件の判定
		if (k - k_update > N_sc){
			//cout << "relax 's k : " << k << endl;
			break;
		}

	}

	//cout << "relax opt : " << opt << endl;
	return opt;
}