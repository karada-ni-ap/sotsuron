#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double local_opt(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly);
pair<MatrixXd, MatrixXd> Wbound(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly);
MatrixXd projectionW(MatrixXd W, MatrixXd Wmax, MatrixXd Wmin);
double norm(VectorXd x, VectorXd y, MatrixXd W);
double relaxation(VectorXd ux, VectorXd lx, VectorXd uy, VectorXd ly);