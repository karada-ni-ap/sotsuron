#include <iostream>
#include <Eigen/Dense>
#include "const.h"

using namespace std;
using namespace Eigen;

void initA();
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L); //U,Y‚Íy0‚ÌBOX
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L); //U,Y‚Íx0‚ÌBOX