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
	
	for (int k = 0; k < ite_local; k++){
		X = sev_y(y, x, ux, lx);
		x = X.second;

		Y = sev_x(x, y, uy, ly);
		y = Y.second;

		pre = tmp;
		tmp = Y.first;

		max = tmp>max ? tmp : max;

		//��������
		if (abs(tmp - pre) < eps_local)
			break;
	}

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

	VectorXd px = VectorXd::Zero(d);
	VectorXd py = VectorXd::Zero(m);
	MatrixXd pW = MatrixXd::Zero(d, m);
	
	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	pair<MatrixXd, MatrixXd> Pair = Wbound(ux, lx, uy, ly);

	VectorXd x = (ux + lx) / 2;
	VectorXd y = (uy + ly) / 2;
	MatrixXd W = (Pair.first + Pair.second) / 2;

	double min = -Inf;
	double pre = -Inf;
	double opt = -Inf;

	double alp0 = 1.0; //norm(ux - lx, uy - ly, Pair.first - Pair.second) / beta;
	double alp = alp0;

	for (int k = 0; k < ite_relax; k++){
		//�X�e�b�v�T�C�Y
		alp = alp0 / (k + 1);

		//prev
		pre = min;
		px  = x;
		py  = y;
		pW  = W;

		//Axy�̍\��
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
		pn = es.eigenvectors().col(0);	//Axy�̍ŏ��ŗL�x�N�g��
		min = es.eigenvalues()(0);		//Axy�̍ŏ��ŗL�l
		opt = min > opt ? min : opt;

		//�㏸����
		for (int i = 1; i <= d; i++){
			dx(i-1) = pn.transpose() * A[i][0] * pn;
		}
		for (int j = 1; j <= m; j++){
			dy(j - 1) = pn.transpose() * A[0][j] * pn;
			for (int i = 1; i <= d; i++){
				dW(i-1,j-1)= pn.transpose() * A[i][j] * pn;
			}
		}

		//x,y,W�̍X�V
		x = projection(x + alp*dx, ux, lx);
		y = projection(y + alp*dy, uy, ly);
		W = projectionW(W + alp*dW, Pair.first, Pair.second);	

		//��������
		if (norm(dx, dy, dW) < eps_relax){
			//cout << "�ɒl��break" << endl;
			break;
		}

		if (abs(pre - min) < eps_relax && norm(px-x, py-y, pW-W) < eps_relax){
			//cout << "���E��break" << endl;
			break;
		}

	}

	return opt;
}