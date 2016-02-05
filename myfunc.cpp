#include "myfunc.h"

double pdf(double x){
	return exp(-0.5 * x * x) / sqrt(2 * M_PI);
}

double cdf(double x){
	return 0.5 * ( 1 + erf(x/sqrt(2)) );
}

double kernel(Eigen::VectorXd x1, Eigen::VectorXd x2){
	return exp( - (x1-x2).squaredNorm() / (2*theta*theta) );
}

VectorXd bound_rand(){
	VectorXd rand = VectorXd::Random(d);
	VectorXd v	  = VectorXd::Zero(d);

	for (int i = 0; i < d; i++){
		v(i) = 0.5 * (Ux0(i) - Lx0(i)) * rand(i) + 0.5 * (Ux0(i) + Lx0(i));
	}
	return v;
}

VectorXd projection(VectorXd z, VectorXd U, VectorXd L){
	VectorXd v = z;
	for (int i = 0; i < z.size(); i++){
		if (z(i) > U(i))		v(i) = U(i);
		else if (z(i) < L(i))	v(i) = L(i);
	}
	return v;
}