#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void initA();
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L); //U,Yはy0のBOX
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L); //U,Yはx0のBOX

pair<double, VectorXd> local_search(VectorXd x0, VectorXd y0); // lsBOで使用