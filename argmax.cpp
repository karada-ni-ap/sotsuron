#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "BO.h"

using namespace std;
using namespace Eigen;

MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y){
	if (s.norm() < sigma_thre || y.norm() < sigma_thre)
		return MatrixXd::Zero(d, d);
	
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

	//dirが降下方向でない
	//Hは正定値であるので，本来ありえないが...
	if (Gd < 0){
		//cout << "Gd is not positive..." << endl;
		return 0;
	}

	//dirが降下方向
	else {
		while (true){
			if (Alp < sigma_thre && u(X + Alp*dir) < sigma_thre)
				return 0;
			else if ( u(X + Alp*dir) >= uX + c1_bfgs*Alp*Gd ) //Armijoの条件を満たす
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
		if (x(i) > Ux(i))		v(i) = Ux(i);
		else if (x(i) < Lx(i))	v(i) = Lx(i);
	}
	return v;
}

VectorXd bfgs(VectorXd x0){
	MatrixXd H = MatrixXd::Identity(d, d);

	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	VectorXd Gold = u_over_x(Xnew);
	VectorXd Gnew = Gold;

	VectorXd dir  = VectorXd::Zero(d);
	double Alp = 0.1;

	for (int k = 0; k < ite_bfgs; k++){
		//収束判定
		if (Gnew.norm() < eps_bfgs		//勾配=0
			|| u(Xnew) < sigma_thre		//u(x)=0
			|| Alp < sigma_thre)		//α=0 
		{
			break;
		}	
		
		//上昇
		dir = H*Gnew;
		Alp = back_track(Xnew, Gnew, dir);

		Xold = Xnew;
		Xnew = Xold + Alp*dir;
		Xnew = projection(Xnew);

		Gold = Gnew;
		Gnew = u_over_x(Xnew);

		H = update_H(H, Xnew - Xold, Gnew - Gold);
	}
	return Xnew;
}

VectorXd argmax_u(){
	if (t == 0) //初期点はランダム
		return bound_rand(d);

	else if (!bfgs_or_rand) // Debug用
		return bound_rand(d);

	else{
		//ここでuの最大化を行う
		double   opt = -Inf;
		double   utmp = 0;
		VectorXd Xopt = VectorXd::Random(d);
		VectorXd Xtmp = VectorXd::Zero(d);

		for (int i = 0; i < num_bfgs; i++){
			Xtmp = bfgs(bound_rand(d));	//BFGS法の初期点はランダム
			utmp = u(Xtmp);
			if (utmp > opt){			//uの値を更新
				opt = utmp;
				Xopt = Xtmp;
			}
		}
		return Xopt;
	}

}