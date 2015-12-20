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

pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L){
	for (int j = 0; j< (m + 1); j++){
		B[j] = A[0][j];
		for (int i = 1; i <= d; i++){
				B[j] += (x(i-1))*A[i][j];
		}
	}

	MatrixXd Axy = MatrixXd::Zero(n, n);
	VectorXd pn  = VectorXd::Zero(n);

	VectorXd y     = y0; //����z�@�̏����_
	VectorXd y_pre = y0;
	VectorXd dir = VectorXd::Zero(m);

	double   min  = 0;
	double   alp0 = (U - L).norm() / beta;
	double   alp  = alp0;

	//y�ɂ��ė���z�@
	for (int k = 0; k < ite_sev; k++){
		//Axy�̍\��
		Axy = B[0];
		for (int j = 1; j <= m; j++)
			Axy += y(j-1)*B[j];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		
		pn  = es.eigenvectors().col(0);	//Axy�̍ŏ��ŗL�x�N�g��
		min = es.eigenvalues()(0);		//Axy�̍ŏ��ŗL�l

		for (int j = 1; j <= m; j++){
			dir(j-1) = half_ip(B[j], pn);
		}

		//y�̍X�V
		y_pre = y;
		y = projection(y + alp*dir, U, L);

		//��������
		if (dir.norm() < eps_sev || (y_pre - y).norm() < eps_sev)
			break;

		alp = alp0 / sqrt(k+1);	//�X�e�b�v�T�C�Y�̍X�V
	}

	return make_pair(min, y_pre);
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

	VectorXd x     = x0; //����z�@�̏����_
	VectorXd x_pre = x0;
	VectorXd dir = VectorXd::Zero(d);

	double min  = 0;
	double alp0 = (U - L).norm() / beta;
	double alp  = alp0;

	//x�ɂ��ė���z�@
	for (int k = 0; k < ite_sev; k++){
		//Axy�̍\��
		Axy = C[0];
		for (int i = 1; i <= d; i++)
			Axy += x(i - 1)*C[i];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);

		pn = es.eigenvectors().col(0);	//Axy�̍ŏ��ŗL�x�N�g��
		min = es.eigenvalues()(0);		//Axy�̍ŏ��ŗL�l

		for (int i = 1; i <= d; i++){
			dir(i - 1) = half_ip(C[i], pn);
		}

		//x�̍X�V
		x_pre = x;
		x = projection(x + alp*dir, U, L);

		//��������
		if (dir.norm() < eps_sev || (x_pre - x).norm() < eps_sev)
			break;

		alp = alp0 / sqrt(k + 1);	//�X�e�b�v�T�C�Y�̍X�V
	}

	return make_pair(min, x_pre);
}