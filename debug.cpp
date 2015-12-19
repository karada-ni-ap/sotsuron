#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "local_and_relax.h"
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

void debug_inside(){ //tはデータセットのサイズ
	int interval = 1;

	VectorXd gntn = VectorXd::Constant(d, 0.0);
	VectorXd mgue = VectorXd::Constant(d, 0.5);
	VectorXd hdrst = VectorXd::Constant(d, -0.5);


	if (t > 0 && t % interval == 0){
		//cout << "sigma is " << sigma(gntn) << endl;
		//cout << "mu is " << mu(gntn) << endl;
		//cout << "k(x) is " << k(mgue).transpose() << endl;

		//cout << u(hdrst) << endl;

		//cout << u_over_k(VectorXd::Constant(d, -0.5)).transpose() << endl;
		//cout << k_over_x(VectorXd::Constant(d, -0.5)) << endl;
		//cout << u_over_x(VectorXd::Constant(d, -0.5)).transpose() << endl;

		//cout << bfgs(mgue).transpose() << endl;
		//cout << argmax_u().transpose() << endl;

		cout << "t is " << t << endl;
		cout << x_next.transpose() << endl;
		//cout << u_over_x(x_next) << endl;
		//cout << back_track(x_next, u_over_x(x_next), u_over_x(x_next)) << endl;
		//cout << "f(t) is " << f(t-1) << endl;

		cout << "-----------------------------------------------" << endl;
	}

}

void debug_last(){
	cout << "maxf  is " << maxf << endl;
	//cout << "relax is " << relaxation(Ux0, Lx0, Uy0, Ly0) << endl;
	//cout << "local is " << local_opt(Ux0, Lx0, Uy0, Ly0) << endl;

	classQ Q0;

	Q0.makeQ(Ux0, Lx0, Uy0, Ly0);
	pair<classQ, classQ> devided = Q0.devide();

	classQ Q1 = devided.first;
	classQ Q2 = devided.second;

	Q1.next = &Q2;
	cout << Q1.next->Ux << endl;

	cout << devided.first.Lx.transpose() << endl;
	cout << devided.second.Ux.transpose() << endl;
	


	for (int z = 0; z < 100; z++){
			//VectorXd x = VectorXd::Random(d);
			//cout << x.transpose() << endl;
			//cout << projection(x).transpose() << endl;

			//VectorXd x = bound_rand(d);
			//VectorXd y = bound_rand(d);
			//VectorXd y = u_over_x(x);
			//cout << y.transpose() << endl;
			//cout << back_track(x, y, y) << endl;
			//cout << u_over_x(x).transpose() << endl;

			//cout << "-----------------------------------------------" << endl;
	}
	
}