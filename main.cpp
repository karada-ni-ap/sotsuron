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
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	double maxf_of_each_d[15];
	int    find_of_each_d[15];

	clock_t finding_time_of_each_d[15];
	clock_t time_of_each_d[15];

	clock_t start;
	clock_t end;

	SDMorBFGS = false;		// trueÅ®SDM / falseÅ®BFGS

	for (d = 5; d <= 15; d++){
		initA();
		initialize_for_BO();

		start = clock();
		maxf_of_each_d[d - 1] = branch_and_cut();
		end = clock();

		finding_time_of_each_d[d-1]		= end - finding_time_BC;
		time_of_each_d[d-1]				= end - start;
		find_of_each_d[d-1]				= k_find;

		deleting_for_BO();
	}

	for (int i = 5; i <= 15; i++){
		cout << "d=" << i << " : " << maxf_of_each_d[i-1] << endl;
		cout << "time : " << time_of_each_d[i - 1] << endl;
		cout << "find : " << finding_time_of_each_d[i - 1] << endl;
		cout << "find : " << find_of_each_d[i - 1] << endl;
		
	}

	return 0;
}