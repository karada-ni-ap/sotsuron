#include <time.h>
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

	//èÄîı
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

	//ï™äÑ
	if (max_x >= max_y){
		VectorXd ex = VectorXd::Zero(d);

		ex(longest_x) = 1;

		max_x /= 2.0;

		Pair.first .makeQ(Ux, Lx + max_x*ex, Uy, Ly);
		Pair.second.makeQ(Ux - max_x*ex, Lx, Uy, Ly);
	}

	else{
		VectorXd ey = VectorXd::Zero(m);

		ey(longest_y) = 1;

		max_y /= 2.0;

		Pair.first .makeQ(Ux, Lx, Uy, Ly + max_y*ey);
		Pair.second.makeQ(Ux, Lx, Uy - max_y*ey, Ly);
	}

	return Pair;
}

void classQ::calculate_lo(){
	Q_L = local_opt(Ux, Lx, Uy, Ly);
}

void Qlist::add(classQ Q){
	size++;

	classQ *Qins;
	Qins = new classQ();

	Qins->Ux = Q.Ux;
	Qins->Lx = Q.Lx;
	Qins->Uy = Q.Uy;
	Qins->Ly = Q.Ly;

	Qins->Q_U = Q.Q_U;
	Qins->Q_L = Q.Q_L;

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
	size--;

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
			size--;
		}
	}
}

double Qlist::maxL(){
	double  max = -Inf;
	classQ* Qtmp = root.next;

	while (Qtmp != &root){
		//cout << Qtmp->Q_U << " : " << Qtmp->Q_L << endl; //Debug
		if (Qtmp->Q_L > max){
			max = Qtmp->Q_L;
		}
		Qtmp = Qtmp->next;
	}

	return max;
}

double branch_and_cut(clock_t* sample_time, double* sample_val){
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

	clock_t start = clock();

	for (int k = 0; k < ite_bc; k++){
		cout << "k is " << k << endl;

		//cout << "List size before extraxt : " << List.size << endl;
		Q = List.extract();

		maxL = List.maxL();

		Pair = Q.devide();
		Q12[0] = Pair.first;
		Q12[1] = Pair.second;

		for (int i = 0; i < 2; i++){
			if (Q12[i].Q_U > maxL){
				Q12[i].calculate_lo();
				List.add(Q12[i]);
				//cout << "Q" << i + 1 << " : " << Q12[i].Q_U << " > " << Q12[i].Q_L << endl;

				//cout << Q12[i].Ux.transpose() << endl;
				//cout << Q12[i].Lx.transpose() << endl;
				//cout << Q12[i].Uy.transpose() << endl;
				//cout << Q12[i].Ly.transpose() << endl;

				if (Q12[i].Q_L > Q12[i].Q_U)
					cout << "L > U (;_;)" << endl;

				if (Q12[i].Q_L > maxL)
					maxL = Q12[i].Q_L;

				if (Q12[i].Q_L > Qopt.Q_L){
					Qopt = Q12[i];
					cout << "Qopt is updated." << endl;
				}
			}
		}

		List.delete_tail(maxL);

		if (List.maxL() != maxL){
			cout << "Q having maxL has been deleted..." << endl;
			maxL = List.maxL();
		}

		sample_time[k] = clock() - start;
		sample_val[k] = Qopt.Q_L;

		maxU = List.root.next->Q_U;

		//Debug
		cout << "maxU " << maxU << endl;
		cout << "maxL " << maxL << endl;
		cout << "--------------------------------" << endl;

		if (maxU - maxL < eps_bc)
			break;
	}

	return Qopt.Q_L;
}