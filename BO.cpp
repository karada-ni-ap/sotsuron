#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "const.h"
#include "myfunc.h"
#include "obj.h"
#include "sev.h"
#include "argmax.h"
#include "debug.h"
#include "local_and_relax.h"

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

void update_m(){
	VectorXd ones = VectorXd::Constant(t, 1.0);
	VectorXd v = Kinv.topLeftCorner(t, t) * ones;

	double a   = v.transpose() * f.head(t);
	double b   = v.transpose() * ones;

	mean = a / b;
}

void update_K(Eigen::VectorXd x){ // x = x_nextを想定
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
		VectorXd m1 = VectorXd::Constant(t, mean);
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
		double gamma = (mu_ - maxf_BO - xi) / sigma_;
		return (mu_ - maxf_BO - xi)*cdf(gamma) + sigma_*pdf(gamma);
	}
}

VectorXd u_over_k(VectorXd x){
	double sigma_ = sigma(x);
	double mu_ = mu(x);
	double gamma = (mu_ - maxf_BO - xi) / sigma_;

	VectorXd m1 = VectorXd::Constant(t,mean);

	VectorXd v = VectorXd::Zero(t);
	v = cdf(gamma)*(f.head(t) - m1.head(t)) - (pdf(gamma) / sigma_)*k(x);

	return Kinv.topLeftCorner(t, t)*v;
}

MatrixXd k_over_x(VectorXd x){
	MatrixXd A = MatrixXd::Zero(d, t);
	for (int j = 0; j < t; j++){
		A.col(j) = kernel(D_q.col(j), x) * (D_q.col(j) - x);
	}
	return A;
}

VectorXd u_over_x(VectorXd x){
	if (sigma(x) < sigma_thre)
		return VectorXd::Zero(d);
	
	else
		return k_over_x(x)*u_over_k(x);
}

pair<double, VectorXd> BO(){
	VectorXd x_next = VectorXd::Zero(d);
	VectorXd x_opt  = VectorXd::Zero(d);

	for (t = 0; t < T; t++){
		//【tはこの時点におけるデータセットのサイズ】//

		update_m();
		x_next = argmax_u();
		debug_inside(x_next, x_opt);

		//データセットの更新
		D_q.col(t) = x_next;

		f(t) = sev_x(x_next, (Uy0 + Ly0) / 2, Uy0, Ly0).first;

		if (f(t)>maxf_BO){
			x_opt = x_next;

			maxf_BO = f(t);
			t_find = t;
			finding_time_BO = clock();
		}
		//【この時点でデータセットのサイズはt+1】//

		update_K(x_next);
		//【updateが行われた後，Kのサイズはt+1】//
	}

	return make_pair(maxf_BO, x_opt);
}

pair<double, VectorXd> lsBO(){
	VectorXd x_next = VectorXd::Zero(d);
	VectorXd x_lo   = VectorXd::Zero(d);

	pair<double, VectorXd> Pair_next;
	pair<double, VectorXd> Pair_lo;

	for (t = 0; t < T; t++){
		cout << "t : " << t << endl;
		//【tはこの時点におけるデータセットのサイズ】//

		update_m();

		cout << "x_next" << endl;

		x_next = argmax_u();
		Pair_next = sev_x(x_next, (Uy0 + Ly0) / 2, Uy0, Ly0); // Pair_next = <f(x_next), y0>

		f(t) = Pair_next.first;
		D_q.col(t) = x_next;
		update_K(x_next);

		if (f(t)>maxf_lsBO){
			t++;
			if (t == T) break;

			cout << "local search" << endl;

			Pair_lo = local_search(x_next, Pair_next.second); // Pair_lo = <f(x_lo), x_lo>
			x_lo = Pair_lo.second;

			f(t) = Pair_lo.first;
			D_q.col(t) = x_lo;
			update_K(x_lo);

			maxf_lsBO = f(t);
			l_find = t;
			finding_time_lsBO = clock();

			t++;
			if (t == T) break;

			cout << "opposite sampling" << endl;

			x_next = 2 * x_lo - x_next; //x_loを挟んで反対側の点
			Pair_next = sev_x(x_next, (Uy0 + Ly0) / 2, Uy0, Ly0);
			
			f(t) = Pair_next.first;
			D_q.col(t) = x_next;
			update_K(x_next);
		}

		cout << "-----------------------------------" << endl;

	}

	return make_pair(maxf_lsBO, x_lo);
}

void initialize_for_BO(){
	t = 0;

	D_q = MatrixXd::Zero(d, T);
	f = VectorXd::Zero(T);
	K = MatrixXd::Zero(T, T);
	Kinv = MatrixXd::Zero(T, T);

	maxf_BO = -Inf;
	maxf_lsBO = -Inf;
}