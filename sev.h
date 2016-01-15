#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void initA();
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L); // U,Y��y0��BOX, <��(x, y*), y*>��Ԃ�
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L); // U,Y��x0��BOX, <��(x*, y), x*>��Ԃ�

pair<double, VectorXd> local_search(VectorXd x0, VectorXd y0); // <��(x*, y*), x*>��Ԃ�