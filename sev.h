#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void initA();
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L); //U,Y��y0��BOX
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L); //U,Y��x0��BOX

pair<double, VectorXd> local_search(VectorXd x0, VectorXd y0); // lsBO�Ŏg�p