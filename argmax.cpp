#include "argmax.h"

MatrixXd update_H(MatrixXd H, VectorXd s, VectorXd y){
	if (y.norm() < eps_H){
		cout << "|y| << 1" << endl;
		//return MatrixXd::Zero(d, d);
		return MatrixXd::Identity(d, d);
	}

	if (y.dot(s) > 0){
		VectorXd a = s.normalized();
		VectorXd b = y.normalized();
		double bta = b.dot(a);
		double c = s.norm() / y.norm();

		MatrixXd P = MatrixXd::Identity(d, d) - (b*a.transpose()) / bta;
		return P*H*P + c * a*a.transpose() / bta;
	}

	/*else{
		double yts = y.transpose()*s;
		VectorXd Hy = H*y;
		MatrixXd Hyst = Hy*s.transpose();
		return
			H
			+ ( (yts + y.transpose()*Hy) / yts )*s*s.transpose()
			- ( Hyst+ Hyst.transpose() ) / yts;
	}*/

	/*else{
		VectorXd a = s.normalized();
		VectorXd b = y.normalized();
		double bta = b.dot(a);

		VectorXd Hb = H*b;

		return H
			- (a*Hb.transpose() + Hb*a.transpose()) / bta;
		+(s.norm() / (bta*y.norm()) + b.dot(Hb) / (bta*bta)) * a * a.transpose();
	}*/

	else{ //正定値が破れる可能性あり
		cout << "<y,s> <= 0" << endl;
		return MatrixXd::Identity(d, d);
	}

}

double back_track(VectorXd X, VectorXd Grad, VectorXd dir){
	double alp = alp0;
	double uX = u(X);
	double Gd = Grad.dot(dir);

	//dirが降下方向でない⇔H＜O
	if (Gd < 0){
		cout << "!!! Gd is negative !!!" << endl;
		return 0;
	}

	//dirが降下方向
	else {
		while (true){
			if (alp < sigma_thre && u(X + alp*dir) < sigma_thre){
				cout << "check check check check check check check check" << endl;
				return 0;
			}
			else if (u(X + alp*dir) >= uX + c1*alp*Gd){ //Armijoの条件を満たす
				if (dir.transpose() * (u_over_x(X + alp*dir) - c2*Grad) > 0){ //Wolfeの条件を満たさない
					cout << "but it doesn't satisfy Wolfe..." << endl;
				}
				break;
			}
			else
				alp *= rho;
		}
		return alp;
	}
}

VectorXd bfgs(VectorXd x0){
	MatrixXd H = MatrixXd::Identity(d, d);

	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	VectorXd Gold = u_over_x(Xnew);
	VectorXd Gnew = Gold;

	if (Gnew.norm() < eps_bfgs){
		cout << "即break" << endl;
		return Xnew;
	}

	VectorXd dir  = VectorXd::Zero(d);
	double alp = 0;

	for (int k = 0; k < ite_bfgs; k++){
		//上昇
		dir = H*Gnew;
		alp = back_track(Xnew, Gnew, dir);

		Xold = Xnew;
		Xnew = projection(Xold + alp*dir, Ux0, Lx0);

		if ((Xnew - Xold).norm() < eps_bfgs){
			//cout << "境界でbreak" << endl;
			break;
		}

		Gold = Gnew;
		Gnew = u_over_x(Xnew);

		if (Gnew.norm() < eps_bfgs){
			//cout << "極値でbreak" << endl;
			break;
		}

		H = update_H(H, Xnew - Xold, Gnew - Gold);
	}
	return Xnew;
}

VectorXd sdm(VectorXd x0){
	VectorXd Xold = x0;
	VectorXd Xnew = x0;

	VectorXd dir = VectorXd::Zero(d);

	double alp = 0;

	for (int k = 0; k < ite_sdm; k++){
		//上昇
		dir = u_over_x(Xnew);

		if (dir.norm() < eps_sdm){
			cout << "極値でbreak" << endl;
			break;
		}

		alp = back_track(Xnew, dir, dir);

		Xold = Xnew;
		Xnew = projection(Xold + alp*dir, Ux0, Lx0);

		if ((Xnew - Xold).norm() < eps_sdm){
			cout << "境界でbreak" << endl;
			break;
		}
	}

	return Xnew;
}

VectorXd argmax_u(){
	if (t == 0) //初期点はランダム
		return bound_rand();

	else if (!bfgs_or_rand) // Debug用
		return bound_rand();

	else{
		//ここでuの最大化を行う
		double   opt = -Inf;
		double   utmp = 0;
		VectorXd Xopt = VectorXd::Random(d);
		VectorXd Xtmp = VectorXd::Zero(d);

		for (int i = 0; i < num_of_start; i++){
			
			//Xtmp = sdm(bound_rand());
			Xtmp = bfgs(bound_rand());
			
			utmp = u(Xtmp);
			if (utmp > opt){
				opt = utmp;
				Xopt = Xtmp;
			}
		}
		return Xopt;
	}

}