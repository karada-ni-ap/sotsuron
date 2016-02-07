#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "argmax.h"
#include "sev.h"
#include "experiment.h"

using namespace std;
using namespace Eigen;

int main(void){
	EIorUCB  = true;
	BOorlsBO = true;

	int dmin = 5;
	int dmax = 5;
	
	int mmin = 10;
	int mmax = 10;

	int nmin = 7;
	int nmax = 7;

	int num = 10;


	for (int d_ = dmin; d_ <= dmax; d_ += 1){
		for (int m_ = mmin; m_ <= mmax; m_ += 1){
			for (int n_ = nmin; n_ <= nmax; n_ +=5){
				cout << "d : " << d_ << ", m : " << m_ << ", n : " << n_ << endl;

				if (BOorlsBO)
					vs_BO(d_, m_, n_, num);

				else
					vs_lsBO(d_, m_, n_, num);

				cout << "----------------------------------------------" << endl;
			}
		}
	}
	

	return 0;
}