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
	VectorXd y   = y0; //����z�@�̏����_
	VectorXd dir = VectorXd::Zero(m);
	double   min = 0;
	double   pre = Inf;
	double   Alp = Alpha0;

	//y�ɂ��ė���z�@
	for (int k = 0; k < ite_subgrad; k++){
		//y�̍X�V
		y += Alp*dir;
		y = projection(y, U, L);

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

		if (dir.norm() < eps_subgrad			//��������
			|| abs(pre - min) < eps_subgrad		//�ق�Ƃ͗ǂ��Ȃ�
		)
			break;

		Alp = Alpha0 / sqrt(k+1);	//�X�e�b�v�T�C�Y�̍X�V
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
	VectorXd x   = x0; //����z�@�̏����_
	VectorXd dir = VectorXd::Zero(d);
	double   min = 0;
	double   pre = Inf;
	double   Alp = Alpha0;

	//x�ɂ��ė���z�@
	for (int k = 0; k < ite_subgrad; k++){
		//x�̍X�V
		x += Alp*dir;
		x = projection(x, U, L);

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

		if (dir.norm() < eps_subgrad			//��������
			|| abs(pre - min) < eps_subgrad		//�ق�Ƃ͂悭�Ȃ�
		)
			break;

		Alp = Alpha0 / sqrt(k + 1);	//�X�e�b�v�T�C�Y�̍X�V
	}

	return make_pair(min, x);
}