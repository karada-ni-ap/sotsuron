#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "branch_and_cut.h"

using namespace std;
using namespace Eigen;

void debug_inside();
void debug_last();