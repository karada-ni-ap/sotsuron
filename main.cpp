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

	for (int d_ = 4; d_ <= 5; d_++){
		for (int m_ = 8; m_ <= 8; m_++){
			for (int n_ = 12; n_ <= 12; n_++){
				cout << "d : " << d_ << ", m : " << m_ << ", n : " << n_ << endl;
				//vs_BO(d_, m_, n_);
				vs_lsBO(d_, m_, n_);
				cout << "----------------------------------------------" << endl;
			}
		}
	}
	

	return 0;
}