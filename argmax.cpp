#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "BO.h"

using namespace std;
using namespace Eigen;

bool judge(VectorXd x){
	if (	(x.array() >= Lb).all() &&
			(x.array() <= Ub).all()		)
		return true;
	else
		return false;
}

MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y){
	double S = y.dot(s);
	double Sinv = 1 / S;

	VectorXd Hy = H*y;
	
	return H
		+ (Sinv + Sinv * Sinv * y.dot(Hy)) * (s * s.transpose())
		- Sinv * ( Hy*s.transpose() + s*Hy.transpose() );
}

//�����ɑ����o�O����I�I
double back_track(VectorXd X, VectorXd Grad, VectorXd dir){
	double Alp  = Alp0;
	double uX = u(X);
	double Gd = Grad.transpose()*dir;

	//cout << dir.transpose() << endl;
	//cout << Gd << endl;


	//dir���~�������łȂ�
	//H��O���C����͂��肦�Ȃ�(����l��)
	if (Gd <= 0){
		cout << "Gd is negative." << endl;
		return Alp0;
	}

	//dir���~������
	else {
		cout << "Gd is positive." << endl;
		while (true){
			if (u(X + Alp*dir) >= uX + c1_bfgs*Alp*Gd)
				break;
			else
				Alp *= rho;
		}
		return Alp;
	}
}

VectorXd projection(VectorXd x){
	VectorXd v = x;
	for (int i = 0; i < d; i++){
		if (x(i) > Ub)		v(i) = Ub;
		else if (x(i) < Lb)	v(i) = Lb;
	}
	return v;
}

VectorXd bfgs(VectorXd x0){
	MatrixXd H = MatrixXd::Identity(d, d);

	//VectorXd Xold = VectorXd::Zero(d);
	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	//VectorXd Gold = VectorXd::Zero(d);
	VectorXd Gold = u_over_x(Xnew);
	VectorXd Gnew = u_over_x(Xnew);

	VectorXd dir  = VectorXd::Zero(d);
	double Alp = 0.1;

	for (int k = 0; k < ite_bfgs; k++){
		dir = H*Gnew;
		//Alp = back_track(Xnew, Gnew, dir);

		Xold = Xnew;
		Xnew = Xold + Alp*dir;
		Xnew = projection(Xnew);

		Gold = Gnew;
		Gnew = u_over_x(Xnew);

		if (Gnew.norm() < eps_bfgs) //��������
			break;

		H = update_H(H, Xnew - Xold, Gnew - Gold);
	}
	return Xnew;
}

VectorXd argmax_u(){
	if (t == 0){
		//�����_�̓����_��
		return bound_rand(d);
	}

	else if (!bfgs_or_rand)
		return bound_rand(d);

	else{
		//������u�̍ő剻���s��
		double   opt = -1000;
		double   utmp = 0;
		VectorXd Xopt = VectorXd::Random(d);
		VectorXd Xtmp = VectorXd::Zero(d);

		for (int i = 0; i < num_bfgs; i++){
			Xtmp = bfgs(bound_rand(d));		//BFGS�@�̏����_�̓����_��
			if (judge(Xtmp)){				//Xtmp�����s�\��
				utmp = u(Xtmp);
				if (utmp > opt){			//u�̒l���X�V
					opt = utmp;
					Xopt = Xtmp;
				}
			}
		}
		return Xopt;
	}

}