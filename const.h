#include <Eigen/Dense>
using namespace Eigen;

extern		 int t; // BO�̔�����
extern const int T; // BO�̔����񐔂̏��

extern const int d; // x�̎���
extern const int m; // y�̎���
extern const int n; // �Ώ̍s��̃T�C�Y

extern const int mean; // GP�̎��O���z

extern double maxf; // �f�[�^�Z�b�g�̍ő�l

extern MatrixXd D_q; // �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;   // �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��