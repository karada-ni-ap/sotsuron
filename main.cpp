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

	int dmin = 3;
	int dmax = 5;
	
	int mmin = 3;
	int mmax = 5;

	int nmin = 5;
	int nmax = 5;

	int num = 10;


	for (int d_ = dmin; d_ <= dmax; d_ += 2){
		for (int m_ = mmin; m_ <= mmax; m_ += 2){
			for (int n_ = nmin; n_ <= nmax; n_ +=1){
				
				if (d_ == 5 && m_ == 3) break;

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