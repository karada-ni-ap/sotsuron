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
	double maxf_of_each_d[5];

	for (d = 1; d <= 5; d++){
		initA();
		initialize_for_BO();
		BO();
		maxf_of_each_d[d - 1] = maxf_BO;
	}

	for (int i = 0; i < 5; i++){
		cout << "d=" << i + 1 << " : " << maxf_of_each_d[i] << endl;
	}

	return 0;
}