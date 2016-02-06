#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "argmax.h"
#include "sev.h"
#include "branch_and_cut.h"
#include "problem.h"

using namespace std;
using namespace Eigen;

void vs_BO(int d_, int m_, int n_, int num);
void vs_lsBO(int d_, int m_, int n_, int num);