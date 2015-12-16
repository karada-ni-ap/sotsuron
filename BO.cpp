#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"
#include "obj.h"

using namespace std;
using namespace Eigen;

extern VectorXd m1 = VectorXd::Constant(T,mean);

Eigen::VectorXd k(Eigen::VectorXd x){
	if (t == 0){
		//【起こり得ない状況】//
		return VectorXd::Zero(1);
	}
	else{
		VectorXd kx = VectorXd::Zero(t);
		for (int i = 0; i < t; i++){
			kx(i) = kernel(x, D_q.col(i));
		}
		return kx;
	}
}

void update_K(Eigen::VectorXd x){
	if (t == 0){
		K(0,0)    = 1;
		Kinv(0,0) = 1;
	}
	else{
		//Kの更新
		VectorXd kx = k(x);
		K.block(0,t,t,1) = kx;
		K.block(t,0,1,t) = kx.transpose();
		K(t, t) = 1;

		//Ainv, Sinvの計算
		MatrixXd Ainv = Kinv.topLeftCorner(t,t);
		double S = 1 - kx.transpose()*Ainv*kx;
		double Sinv = 1 / S;
		
		//Kinvの更新
		Kinv.topLeftCorner(t, t) = Ainv + Sinv*(Ainv * kx * kx.transpose() * Ainv);
		VectorXd binv = -Sinv*Ainv*kx;
		Kinv.block(0, t, t, 1) = binv;
		Kinv.block(t, 0, 1, t) = binv.transpose();
		Kinv(t, t) = Sinv;
	}
}

double mu(Eigen::VectorXd x){
	if (t == 0)
		return mean;

	else{
		VectorXd kx = k(x);
		return mean + kx.transpose() * Kinv.topLeftCorner(t, t) * (f.head(t)-m1.head(t));
	}
}


double sigma(Eigen::VectorXd x){
	VectorXd kx = k(x);
	double sigma2 = 1 - kx.transpose() * Kinv.topLeftCorner(t, t) * kx;

	if (sigma2 < 0)
		return 0;
	else
		return sqrt(sigma2);
}


double u(Eigen::VectorXd x){
	double sigma_ = sigma(x);

	if (sigma_ < sigma_thre)
		return 0;

	else{
		double mu_ = mu(x);
		double gamma = (mu_ - maxf) / sigma_;
		return (mu_ - maxf)*cdf(gamma) + sigma_*pdf(gamma);
	}
}

VectorXd u_over_k(VectorXd x){
	double sigma_ = sigma(x);

	if (sigma_ < sigma_thre)
		return VectorXd::Zero(t);

	else{
		double mu_ = mu(x);
		double gamma = (mu_ - maxf) / sigma_;

		VectorXd v = VectorXd::Zero(t);
		v = cdf(gamma)*(f.head(t) - m1.head(t)) - (pdf(gamma) / sigma_)*k(x);

		return Kinv.topLeftCorner(t, t)*v;
	}
}

MatrixXd k_over_x(VectorXd x){
	MatrixXd A = MatrixXd::Zero(d, t);
	for (int j = 0; j < t; j++){
		A.col(j) = kernel(D_q.col(j), x) * (D_q.col(j) - x);
	}
	return A;
}

VectorXd u_over_x(VectorXd x){
	return k_over_x(x)*u_over_k(x);
}