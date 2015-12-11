#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "myfunc.h"

using namespace std;
using namespace Eigen;

extern VectorXd m1 = VectorXd::Constant(T,mean);

Eigen::VectorXd k(Eigen::VectorXd x){
	if (t == 0){
		//�y�N���蓾�Ȃ��󋵁z//
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
		//K�̍X�V
		VectorXd kx = k(x);
		K.block(0,t,t,1) = kx;
		K.block(t,0,1,t) = kx.transpose();
		K(t, t) = 1;

		//Ainv, Sinv�̌v�Z
		MatrixXd Ainv = Kinv.topLeftCorner(t,t);
		double S = 1 - kx.transpose()*Ainv*kx;
		double Sinv = 1 / S;
		
		//Kinv�̍X�V
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
	return sqrt(1 - kx.transpose() * Kinv.topLeftCorner(t, t) * kx);
}


double u(Eigen::VectorXd x){
	double mu_ = mu(x);
	double sigma_ = sigma(x);
	double gamma = (mu_ - maxf) / sigma_;
	return (mu_ - maxf)*cdf(gamma) + sigma_*pdf(gamma);
}

VectorXd u_over_k(VectorXd x){
	double mu_ = mu(x);
	double sigma_ = sigma(x);
	double gamma = (mu_ - maxf) / sigma_;

	VectorXd v = VectorXd::Zero(t);
	v = cdf(gamma)*(f.head(t) - m1.head(t)) - (pdf(gamma) / sigma_)*k(x);

	return Kinv.topLeftCorner(t, t)*v;
}

Eigen::VectorXd argmax_u(){
	if (t == 0){
		//�����_�̓����_��
		return bound_rand(d);
	}
	else{
		//������u�̍ő剻���s��
		return bound_rand(d);
		//return VectorXd::Random(d);
	}
}