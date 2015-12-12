#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"

using namespace std;
using namespace Eigen;

void debug(){
	int interval = 3;

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

		cout << bfgs(mgue).transpose() << endl;
		//cout << argmax_u().transpose() << endl;

		cout << "-----------------------------------------------" << endl;
	}

}