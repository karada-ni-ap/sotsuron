#include <Eigen/Dense>
using namespace Eigen;

extern		 int t; // BO�̔�����
extern const int T; // BO�̔����񐔂̏��

extern const int d; // x�̎���
extern const int m; // y�̎���
extern const int n; // �Ώ̍s��̃T�C�Y

extern const double Ub;  // �������̂̏�E
extern const double Lb;  // �������̂̉��E

extern const double Alp;		 // BFGS�̃X�e�b�v�T�C�Y
//extern const double conv_thre; // BFGS�̎��������臒l
extern const int    ite_bfgs;	 // BFGS�̔�����
extern const int	num_bfgs;	 // BFGS�̏����_�̐�

extern const int    mean;		// GP�̎��O���z
extern const double sigma_thre; // u=0�Ƃ���臒l

extern double maxf; // �f�[�^�Z�b�g�̍ő�l

extern MatrixXd D_q; // �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;   // �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��

extern const bool bfgs_or_rand; // Debug�p