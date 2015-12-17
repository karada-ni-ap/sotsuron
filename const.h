#include <Eigen/Dense>
using namespace Eigen;

extern const int Inf; //�傫���l

extern		 int t; // BO�̔�����
extern const int T; // BO�̔����񐔂̏��

extern const int d; // x�̎���
extern const int m; // y�̎���
extern const int n; // �Ώ̍s��̃T�C�Y

extern const VectorXd Ux;
extern const VectorXd Lx;
extern const VectorXd Uy;
extern const VectorXd Ly;

extern const double Ub;  // �������̂̏�E
extern const double Lb;  // �������̂̉��E

extern const double Alp0;		 // �o�b�N�g�b�v�@�̃��̏����l
extern const double rho;		 // �o�b�N�g���b�N�@�̃�
extern const double c1_bfgs;	 // �o�b�N�g���b�N�@��c1

extern const double eps_bfgs;	 // BFGS�̎��������臒l
extern const int    ite_bfgs;	 // BFGS�̔�����
extern const int	num_bfgs;	 // BFGS�̏����_�̐�

extern const int    mean;		// GP�̎��O���z
extern const double sigma_thre; // u=0�Ƃ���臒l

extern double maxf;			// �f�[�^�Z�b�g�̍ő�l

extern VectorXd x_next;		// ���̃T���v���_
extern VectorXd x_opt;		// �œK��

extern MatrixXd D_q;		// �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;			// �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��

extern MatrixXd** A; //���͍s��
extern MatrixXd*  B; //i�ɂ��Ă�sum
extern MatrixXd*  C; //j�ɂ��Ă�sum

extern const int    ite_subgrad; //����z�@�̔�����
extern const double Alp_subgrad; //����z�@�̃X�e�b�v�T�C�Y�̏����l
extern const double eps_subgrad; //����z�@�̃�

extern const int    ite_local; //�Ǐ��œK�l�֐��̔�����
extern const double Alp_local; //�Ǐ��œK�l�֐��̏����X�e�b�v�T�C�Y
extern const double eps_local; //�Ǐ��œK�l�֐��̃�

extern const int    ite_relax; //�ɘa���̔�����
extern const double Alp_relax; //�ɘa���̏����X�e�b�v�T�C�Y
extern const double eps_relax; //�ɘa���̃�

extern const int select;		// �ړI�֐��̐؂�ւ�
extern const bool bfgs_or_rand;	// argmax�̐؂�ւ��iBFGS�̃����_���j