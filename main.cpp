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
	initA();

	double max_with_BO = 0;
	double max_with_BC = 0;

	max_with_BO = BO();
	cout << "max with BO is " << max_with_BO << endl;

	max_with_BC = branch_and_cut();
	cout << "max with BC is " << max_with_BC << endl;

	debug_last();

	return 0;
}