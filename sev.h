#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void initA();
pair<double, VectorXd> sev_x(VectorXd x, VectorXd y0, VectorXd U, VectorXd L); // U,Y‚Íy0‚ÌBOX, <ƒÕ(x, y*), y*>‚ð•Ô‚·
pair<double, VectorXd> sev_y(VectorXd y, VectorXd x0, VectorXd U, VectorXd L); // U,Y‚Íx0‚ÌBOX, <ƒÕ(x*, y), x*>‚ð•Ô‚·

pair<double, VectorXd> local_search(VectorXd x0, VectorXd y0); // <ƒÕ(x*, y*), x*>‚ð•Ô‚·