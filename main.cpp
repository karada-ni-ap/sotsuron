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
	int dmin = 5;
	int dmax = 6;
	
	int mmin = 8;
	int mmax = 8;

	int nmin = 15;
	int nmax = 15;

	int num = 2;


	for (int d_ = dmin; d_ <= dmax; d_++){
		for (int m_ = mmin; m_ <= mmax; m_++){
			for (int n_ = nmin; n_ <= nmax; n_++){
				cout << "d : " << d_ << ", m : " << m_ << ", n : " << n_ << endl;

				vs_BO(d_, m_, n_, num);
				//vs_lsBO(d_, m_, n_);

				cout << "----------------------------------------------" << endl;
			}
		}
	}
	

	return 0;
}