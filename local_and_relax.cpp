#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"
#include "argmax.h"

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

double relaxation(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	VectorXd x = (ux + lx) / 2;
	VectorXd y = (uy + ly) / 2;
	MatrixXd W = x*y.transpose();

	VectorXd dx = VectorXd::Zero(d);
	VectorXd dy = VectorXd::Zero(m);
	MatrixXd dW = MatrixXd::Zero(d, m);

	VectorXd px = VectorXd::Zero(d);
	VectorXd py = VectorXd::Zero(m);
	MatrixXd pW = MatrixXd::Zero(d, m);
	
	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	double min = 0;
	double opt = -Inf;

	double alpx0 = (ux - lx).norm() / beta;
	double alpx  = alpx0;

	double alpy0 = (uy - ly).norm() / beta;
	double alpy  = alpy0;

	double alpW  = alpW0;

	for (int k = 0; k < ite_relax; k++){
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
			dx(i-1) = half_ip(A[i][0],pn);
		}
		for (int j = 1; j <= m; j++){
			dy(j - 1) = half_ip(A[0][j],pn);
			for (int i = 1; i <= d; i++){
				dW(i-1,j-1)= half_ip(A[i][j],pn);
			}
		}

		//x,y,W�̍X�V
		px = x;
		py = y;
		pW = W;

		x = projection(x + alpx*dx, ux, lx);
		y = projection(y + alpy*dy, uy, ly);
		W = W + alpW*dW;
		

		//��������
		if (dx.norm() < eps_relax && dy.norm() < eps_relax)
			break;

		if ((px - x).norm() < eps_relax && (py - y).norm() < eps_relax) //�����Ŗ���E�o���Ă���
			break;

		alpx = alpx0 / sqrt(k + 1);
		alpy = alpy0 / sqrt(k + 1);
		alpW = alpW0 / sqrt(k + 1);

	}

	return opt;
}