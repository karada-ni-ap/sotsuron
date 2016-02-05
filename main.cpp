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
	int dmin = 10;
	int dmax = 10;
	
	int mmin = 15;
	int mmax = 15;

	int nmin = 10;
	int nmax = 10;

	int num = 3;


	for (int d_ = dmin; d_ <= dmax; d_++){
		for (int m_ = mmin; m_ <= mmax; m_++){
			for (int n_ = nmin; n_ <= nmax; n_++){
				cout << "d : " << d_ << ", m : " << m_ << ", n : " << n_ << endl;

				//vs_BO(d_, m_, n_, num);
				vs_lsBO(d_, m_, n_, num);

				cout << "----------------------------------------------" << endl;
			}
		}
	}
	

	return 0;
}