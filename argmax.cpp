#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "BO.h"

using namespace std;
using namespace Eigen;

MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y){
	if (y.norm() < eps_H){
		cout << "y << 1" << endl;
		//return MatrixXd::Zero(d, d);
		return MatrixXd::Identity(d, d);
	}
	
	else{
		VectorXd a = s.normalized();
		VectorXd b = y.normalized();
		double bta = b.dot(a);

		VectorXd Hb = H*b;

		return H
			- (a*Hb.transpose() + Hb*a.transpose()) / bta;
		+(s.norm() / (bta*y.norm()) + b.dot(Hb) / (bta*bta)) * a * a.transpose();
	}
}

double back_track(VectorXd X, VectorXd Grad, VectorXd dir){
	double Alp = Alp0;
	double uX = u(X);
	double Gd = Grad.dot(dir);

	//dir���~�������łȂ�
	//H�͐���l�ł���̂ŁC�{�����肦�Ȃ���...
	if (Gd < 0){
		cout << "Gd is negative" << endl;
		return 0;
	}

	//dir���~������
	else {
		while (true){
			if (Alp < sigma_thre && u(X + Alp*dir) < sigma_thre)
				return 0;
			else if ( u(X + Alp*dir) >= uX + c1*Alp*Gd ) //Armijo�̏����𖞂���
				break;
			else
				Alp *= rho;
		}
		return Alp;
	}
}

VectorXd bfgs(VectorXd x0){
	MatrixXd H = MatrixXd::Identity(d, d);

	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	VectorXd Gold = u_over_x(Xnew);
	VectorXd Gnew = Gold;

	if (Gnew.norm() < eps_bfgs){
		cout << "��break" << endl;
		return Xnew;
	}

	VectorXd dir  = VectorXd::Zero(d);
	double alp = 0;

	for (int k = 0; k < ite_bfgs; k++){
		//�㏸
		dir = H*Gnew;
		alp = back_track(Xnew, Gnew, dir);

		Xold = Xnew;
		Xnew = projection(Xold + alp*dir, Ux0, Lx0);

		if ((Xnew - Xold).norm() < eps_bfgs){
			cout << "���E��break" << endl;
			break;
		}

		Gold = Gnew;
		Gnew = u_over_x(Xnew);

		if (Gnew.norm() < eps_bfgs){
			cout << "�ɒl��break" << endl;
			break;
		}

		H = update_H(H, Xnew - Xold, Gnew - Gold);
	}
	return Xnew;
}

VectorXd sdm(VectorXd x0){
	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	VectorXd dir = VectorXd::Zero(d);

	double alp = 0;

	for (int k = 0; k < ite_sdm; k++){
		//�㏸
		dir = u_over_x(Xnew);

		if (dir.norm() < eps_sdm){
			cout << "�ɒl��break" << endl;
			break;
		}

		alp = back_track(Xnew, dir, dir);

		Xold = Xnew;
		Xnew = projection(Xold + alp*dir, Ux0, Lx0);

		if ((Xnew - Xold).norm() < eps_sdm){
			cout << "���E��break" << endl;
			break;
		}
	}

	return Xnew;
}

VectorXd argmax_u(){
	if (t == 0) //�����_�̓����_��
		return bound_rand();

	else if (!bfgs_or_rand) // Debug�p
		return bound_rand();

	else{
		//������u�̍ő剻���s��
		double   opt = -Inf;
		double   utmp = 0;
		VectorXd Xopt = VectorXd::Random(d);
		VectorXd Xtmp = VectorXd::Zero(d);

		for (int i = 0; i < num_of_start; i++){
			
			//Xtmp = sdm(bound_rand());
			Xtmp = bfgs(bound_rand());
			
			utmp = u(Xtmp);
			if (utmp > opt){
				opt = utmp;
				Xopt = Xtmp;
			}
		}
		return Xopt;
	}

}