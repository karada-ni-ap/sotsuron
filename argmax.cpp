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

//id
MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y){
	double S = y.dot(s);
	double Sinv = 1 / S;

	VectorXd Hy = H*y;
	
	return H
		+ (Sinv + Sinv * Sinv * y.dot(Hy)) * (s * s.transpose())
		- Sinv * ( Hy*s.transpose() + s*Hy.transpose() );
}

//id
double back_track(VectorXd X, VectorXd Grad, VectorXd dir){
	double Alp  = Alp0;
	double uX = u(X);
	double Gd = Grad.dot(dir);

	while(true){
		if (u(X + Alp*dir) >= uX + Alp*Gd)
			break;
		else
			Alp *= rho;
	}
	return Alp;
}

//id？
VectorXd bfgs(VectorXd x0){
	MatrixXd H = MatrixXd::Identity(d, d);

	VectorXd Xold = VectorXd::Zero(d);
	VectorXd Xnew = x0;

	VectorXd Gold = VectorXd::Zero(d);
	VectorXd Gnew = u_over_x(Xnew);

	VectorXd dir  = VectorXd::Zero(d);
	double Alp = 0.05;

	for (int k = 0; k < ite_bfgs; k++){
		dir = H*Gnew;

		//Alp = back_track(Xnew, Gnew, dir);

		Xold = Xnew;
		Xnew = Xold + Alp*dir;

		Gold = Gnew;
		Gnew = u_over_x(Xnew);

		if (Gnew.norm() < eps_bfgs)
			break;

		H = update_H(H, Xnew - Xold, Gnew - Gold);
	}
	return Xnew;
}

VectorXd argmax_u(){
	if (t == 0){
		//初期点はランダム
		return bound_rand(d);
	}

	else if (!bfgs_or_rand)
		return bound_rand(d);

	else{
		//ここでuの最大化を行う
		double   opt = -1000;
		double   utmp = 0;
		VectorXd Xopt = VectorXd::Random(d);
		VectorXd Xtmp = VectorXd::Zero(d);

		for (int i = 0; i < num_bfgs; i++){
			Xtmp = bfgs(bound_rand(d));		//BFGS法の初期点はランダム
			if (judge(Xtmp)){				//Xtmpが実行可能解
				utmp = u(Xtmp);
				if (utmp > opt){			//uの値を更新
					opt = utmp;
					Xopt = Xtmp;
				}
			}
		}
		return Xopt;
	}

}