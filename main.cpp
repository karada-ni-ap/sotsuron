#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "debug.h"
#include "sev.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	double maxf_of_each_d[15];

	BO_or_lsBO = false;		// trueÅ®BO  / falseÅ®lsBO
	SDMorBFGS = true;		// trueÅ®SDM / falseÅ®BFGS

	for (d = 7; d <= 15; d++){
		initA();
		initialize_for_BO();

		if (BO_or_lsBO){
			BO();
			maxf_of_each_d[d - 1] = maxf_BO;
		}

		else{
			lsBO();
			maxf_of_each_d[d - 1] = maxf_lsBO;
		}

		deleting_for_BO();
	}

	for (int i = 1; i <= 15; i++){
		cout << "d=" << i << " : " << maxf_of_each_d[i-1] << endl;
	}

	return 0;
}