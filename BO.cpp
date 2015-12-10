#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"

using namespace std;
using namespace Eigen;


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

void updateK(Eigen::VectorXd x){
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
	return 0;
}


double sigma2(Eigen::VectorXd x){
	return 0;
}


double u(Eigen::VectorXd x){
	return 0;
}

Eigen::VectorXd argmax_u(){
	return VectorXd::Random(d);
}