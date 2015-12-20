#include <Eigen/Dense>
using namespace Eigen;

extern const int Inf; //�傫���l

extern		 int t;			// BO�̔�����
extern       int t_find;	// maxf��T�������Ƃ��̔�����
extern const int T;			// BO�̔����񐔂̏��

extern const int d; // x�̎���
extern const int m; // y�̎���
extern const int n; // �Ώ̍s��̃T�C�Y

extern const VectorXd Ux0;
extern const VectorXd Lx0;
extern const VectorXd Uy0;
extern const VectorXd Ly0;

extern const double Alp0;		 // �o�b�N�g���b�N�@�̃��̏����l
extern const double rho;		 // �o�b�N�g���b�N�@�̃�
extern const double c1_bfgs;	 // �o�b�N�g���b�N�@��c1

extern const double eps_bfgs;	 // BFGS�̎��������臒l
extern const int    ite_bfgs;	 // BFGS�̔�����
extern const int	num_bfgs;	 // BFGS�̏����_�̐�

extern const int    mean;		// GP�̎��O���z
extern const double sigma_thre; // u=0�Ƃ���臒l

extern double maxf;			// �f�[�^�Z�b�g�̍ő�l

extern MatrixXd D_q;		// �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;			// �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��

extern MatrixXd** A; //���͍s��
extern MatrixXd*  B; //i�ɂ��Ă�sum
extern MatrixXd*  C; //j�ɂ��Ă�sum

extern const double beta; //����z�@�̃X�e�b�v�T�C�Y�̔{��

extern const int    ite_sev; //����z�@�̔�����
extern const double eps_sev; //����z�@�̃�

extern const int    ite_local; //�Ǐ��œK�l�֐��̔�����
extern const double eps_local; //�Ǐ��œK�l�֐��̃�

extern const int    ite_relax; //�ɘa���̔�����
extern const double eps_relax; //�ɘa���̃�
extern const double alpW0;     //�ɘa����W�̏����X�e�b�v�T�C�Y

extern const int    ite_bc;	//���}����@�̔�����
extern const double eps_bc; //���}����@�̃�

extern const int  select;		// �ړI�֐��̐؂�ւ�
extern const bool bfgs_or_rand;	// argmax�̐؂�ւ��iBFGS�̃����_���j