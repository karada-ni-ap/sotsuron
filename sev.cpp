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
			// MatrixXd tmp = MatrixXd::Random(n, n);
			// A[i][j] = 0.5*(tmp + tmp.transpose()); // ���ꂾ�ƈ�l���z�ɂȂ�Ȃ�
			A[i][j] = MatrixXd::Random(n, n);
			for (int k = 1; k <= n - 1; k++)
				for (int l = 0; l <= n - 2; l++)
					A[i][j](k, l) = A[i][j](l, k);
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

	VectorXd y     = y0; //����z�@�̏����_
	VectorXd y_pre = y0;
	VectorXd dir = VectorXd::Zero(m);

	double min     = 0;
	double min_pre = 0;

	double   alp0 = 1.0;
	double   alp  = alp0;

	//y�ɂ��ė���z�@
	for (int k = 0; k < ite_sev; k++){
		alp     = alp0 / (k + 1);
		min_pre = min;
		y_pre   = y;

		//Axy�̍\��
		Axy = B[0];
		for (int j = 1; j <= m; j++)
			Axy += y(j-1)*B[j];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);
		
		pn  = es.eigenvectors().col(0);	//Axy�̍ŏ��ŗL�x�N�g��
		min = es.eigenvalues()(0);		//Axy�̍ŏ��ŗL�l

		for (int j = 1; j <= m; j++){
			dir(j - 1) = pn.transpose() * B[j] * pn;
		}

		//y�̍X�V
		y = projection(y + alp*dir, U, L);

		//��������
		if (dir.norm() < eps_sev){
			//cout << "�ɒl��break" << endl;
			break;
		}

		if ((y_pre - y).norm() < eps_sev && abs(min - min_pre) < eps_sev){
			//cout << "���E��break" << endl;
			break;
		}

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

	double min     = 0;
	double min_pre = 0;

	double alp0 = 1.0;
	double alp  = alp0;

	//x�ɂ��ė���z�@
	for (int k = 0; k < ite_sev; k++){
		alp = alp0 / (k + 1);
		min_pre = min;
		x_pre = x;

		//Axy�̍\��
		Axy = C[0];
		for (int i = 1; i <= d; i++)
			Axy += x(i - 1)*C[i];

		SelfAdjointEigenSolver<MatrixXd> es(Axy);

		pn = es.eigenvectors().col(0);	//Axy�̍ŏ��ŗL�x�N�g��
		min = es.eigenvalues()(0);		//Axy�̍ŏ��ŗL�l

		for (int i = 1; i <= d; i++){
			dir(i - 1) = pn.transpose() * C[i] * pn;
		}

		//x�̍X�V
		x = projection(x + alp*dir, U, L);

		//��������
		if (dir.norm() < eps_sev){
			//cout << "�ɒl��break" << endl;
			break;
		}
			
		if ((x_pre - x).norm() < eps_sev && abs(min - min_pre) < eps_sev){
			//cout << "���E��break" << endl;
			break;
		}

	}

	return make_pair(min, x_pre);
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

		//��������
		if (abs(tmp - pre) < eps_local)
			break;
	}

	return make_pair(tmp, x); //�Ǐ��œK��(x*, y*)��x*��Ԃ�
}