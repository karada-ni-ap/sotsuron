#include <Eigen/Dense>

Eigen::VectorXd k(Eigen::VectorXd x);
void update_m();
void update_K(Eigen::VectorXd x);

double mu(Eigen::VectorXd x);
double sigma(Eigen::VectorXd x);
double u(Eigen::VectorXd x);

VectorXd u_over_k(VectorXd x);
MatrixXd k_over_x(VectorXd x);
VectorXd u_over_x(VectorXd x);

pair<double, VectorXd> BO();
pair<double, VectorXd> lsBO();

void initialize_for_BO();
void deleting_for_BO();