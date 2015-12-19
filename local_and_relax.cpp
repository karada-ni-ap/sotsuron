#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"
#include "argmax.h"

using namespace std;
using namespace Eigen;

double local_opt(VectorXd Ux0, VectorXd Lx0, VectorXd Uy0, VectorXd Ly0){
	VectorXd x = (Ux0 + Lx0) / 2;
	VectorXd y = (Uy0 + Ly0) / 2;

	pair<double, VectorXd> X;
	pair<double, VectorXd> Y;

	double tmp = 0;
	double max = -Inf;
	double pre = -Inf;

	double Alp = Alp_local;
	
	for (int k = 0; k < ite_local; k++){
		X = sev_y(y, x, Ux0, Lx0, Alp);
		x = X.second;

		Y = sev_x(x, y, Uy0, Ly0, Alp);
		y = Y.second;
		
		tmp = Y.first;

		max = tmp>max ? tmp : max;
		if (abs(tmp - pre) < eps_local)	//ほんとはよくない
			break;

		Alp = Alp_local / sqrt(k + 1);
	}

	return max;
}

double relaxation(VectorXd Ux0, VectorXd Lx0, VectorXd Uy0, VectorXd Ly0){
	VectorXd x = (Ux0 + Lx0) / 2;
	VectorXd y = (Uy0 + Ly0) / 2;
	MatrixXd W = MatrixXd::Zero(d, m);

	VectorXd dx = VectorXd::Zero(d);
	VectorXd dy = VectorXd::Zero(m);
	MatrixXd dW = MatrixXd::Zero(d, m);
	
	MatrixXd Axy = MatrixXd::Zero(n, n);

	VectorXd pn = VectorXd::Zero(n);

	double min = 0;
	double pre = -Inf;
	double opt = -Inf;

	double Alp = Alp_relax;

	for (int k = 0; k < ite_relax; k++){
		//x,y,Wの更新
		x = projection(x + Alp*dx, Ux0, Lx0);
		y = projection(y + Alp*dy, Uy0, Ly0);
		W += Alp*dW;

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
		opt = min > opt ? min : opt;

		//上昇方向
		for (int i = 1; i <= d; i++){
			dx(i-1) = half_ip(A[i][0],pn);
		}
		for (int j = 1; j <= m; j++){
			dy(j - 1) = half_ip(A[0][j],pn);
			for (int i = 1; i <= d; i++){
				dW(i-1,j-1)= half_ip(A[i][j],pn);
			}
		}
		
		//収束判定
		if (dx.norm() < eps_relax
			|| dy.norm() < eps_relax
			|| dW.norm() < eps_relax
			|| abs(pre - min) < eps_relax		//ほんとは良くない
			)
			break;

		Alp = Alp_relax / sqrt(k + 1);	//ステップサイズの更新
		pre = min;
	}

	return opt;
}