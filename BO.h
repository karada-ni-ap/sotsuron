#include <Eigen/Dense>

Eigen::VectorXd k(Eigen::VectorXd x);
void update_m();
void update_K(Eigen::VectorXd x);

double mu(Eigen::VectorXd x, Eigen::VectorXd Kinvk);
double sigma(Eigen::VectorXd x, Eigen::VectorXd Kinvk);
double u(Eigen::VectorXd x);

VectorXd u_over_k(VectorXd x, double sigma_, VectorXd Kinvk);
MatrixXd k_over_x(VectorXd x);
VectorXd u_over_x(VectorXd x);

pair<double, VectorXd> BO(clock_t* sample_time, double* sample_val);
pair<double, VectorXd> lsBO(clock_t* sample_time, double* sample_val);

void initialize_for_BO();
void deleting_for_BO();