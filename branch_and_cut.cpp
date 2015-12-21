#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include "const.h"
#include "sev.h"
#include "local_and_relax.h"
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

void classQ::makeQ(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly){
	Ux = ux;
	Lx = lx;
	Uy = uy;
	Ly = ly;

	Q_U = relaxation(ux, lx, uy, ly);
}

pair<classQ, classQ> classQ::devide(){
	pair<classQ, classQ> Pair;

	int longest_x = 0;
	int longest_y = 0;
	double max_x = -1;
	double max_y = -1;

	for (int i = 0; i < d; i++){
		if (Ux(i) - Lx(i) > max_x){
			max_x = Ux(i) - Lx(i);
			longest_x = i;
		}
	}

	for (int j = 0; j < m; j++){
		if (Uy(j) - Ly(j) > max_y){
			max_y = Uy(j) - Ly(j);
			longest_y = j;
		}
	}

	VectorXd ex = VectorXd::Zero(d);
	VectorXd ey = VectorXd::Zero(m);

	ex(longest_x) = 1;
	ey(longest_y) = 1;

	max_x /= 2.0;
	max_y /= 2.0;

	Pair.first .makeQ(Ux, Lx + max_x*ex, Uy, Ly + max_y*ey);
	Pair.second.makeQ(Ux - max_x*ex, Lx, Uy - max_y*ey, Ly);

	return Pair;
}

void classQ::calculate_lo(){
	Q_L = local_opt(Ux, Lx, Uy, Ly);
}

//–â‘è‚È‚µ
void Qlist::add(classQ Q){
	classQ *Qins;
	Qins = new classQ();

	Qins->Ux = Q.Ux;
	Qins->Lx = Q.Lx;
	Qins->Uy = Q.Uy;
	Qins->Ly = Q.Ly;
	Qins->Q_U = Q.Q_U;

	classQ* Qtmp = root.next;
	classQ* Qprev = &root;

	while (Qtmp != &root){
		if (Qins->Q_U > Qtmp->Q_U){
			break;
		}

		Qprev = Qtmp;
		Qtmp = Qtmp->next;
	}

	Qprev->next = Qins;
	Qins->prev = Qprev;

	Qtmp->prev = Qins;
	Qins->next = Qtmp;
}

classQ Qlist::extract(){
	classQ* Q;
	Q = new classQ();

	classQ* Qret = root.next;

	Qret->next->prev = &root;
	root.next = Qret->next;

	Q->Ux = Qret->Ux;
	Q->Lx = Qret->Lx;
	Q->Uy = Qret->Uy;
	Q->Ly = Qret->Ly;
	
	Q->Q_U = Qret->Q_U;
	Q->Q_L = Qret->Q_L;

	delete Qret;
	return *Q;
}

void Qlist::delete_tail(double L){
	classQ* Qtmp;

	while (true){
		Qtmp = root.prev;

		if (Qtmp == &root || Qtmp->Q_U > L)
			break;

		else{
			Qtmp->prev->next = Qtmp->next;
			Qtmp->next->prev = Qtmp->prev;
			delete Qtmp;
		}
	}
}

double branch_and_cut(){
	classQ Q0;
	classQ Q, Qopt;
	classQ Q12[2];

	pair<classQ, classQ> Pair;
	
	Q0.makeQ(Ux0, Lx0, Uy0, Ly0);
	Q0.calculate_lo();

	Qopt = Q0;

	double maxU = Q0.Q_U;
	double maxL = Q0.Q_L;

	Qlist List;
	List.add(Q0);

	for (int k = 0; k < ite_bc; k++){
		cout << "k is " << k << endl;

		Q = List.extract();
		Pair = Q.devide();

		Q12[0] = Pair.first;
		Q12[1] = Pair.second;

		for (int i = 0; i < 2; i++){
			if (Q12[i].Q_U > maxL){
				Q12[i].calculate_lo();
				List.add(Q12[i]);

				//cout << "Q" << i+1 << " : " << Q12[i].Q_L << endl;

				if (Q12[i].Q_L > maxL){
					maxL = Q12[i].Q_L;
					Qopt = Q12[i];
					cout << "maxL is updated." << endl;
				}
			}
		}

		List.delete_tail(maxL);

		maxU = List.root.next->Q_U;

		//Debug
		cout << "maxU " << maxU << endl;
		//cout << "maxL " << maxL << endl;
		cout << "--------------------------------" << endl;

		if (maxU - maxL < eps_bc)
			break;
	}

	return maxL;
}