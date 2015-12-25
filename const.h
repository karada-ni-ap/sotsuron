#include <iostream>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

extern const int Inf; //�傫���l

extern		 int t;			// BO�̔�����
extern const int T;			// BO�̔����񐔂̏��

extern const int d; // x�̎���
extern const int m; // y�̎���
extern const int n; // �Ώ̍s��̃T�C�Y

extern const VectorXd Ux0;
extern const VectorXd Lx0;
extern const VectorXd Uy0;
extern const VectorXd Ly0;

extern const double alp0;	// �o�b�N�g���b�N�@�̃��̏����l
extern const double rho;	// �o�b�N�g���b�N�@�̃�
extern const double c1;		// �o�b�N�g���b�N�@��c1
extern const double c2;		// �o�b�N�g���b�N�@��c2

extern const double eps_H; // H�����̃�

extern const int    ite_bfgs;	 // BFGS�̔�����
extern const double eps_bfgs;	 // BFGS�̎��������臒l

extern const int    ite_sdm; // �ŋ}�~���@�̔�����
extern const double eps_sdm; // �ŋ}�~���@�̎��������臒l

extern const int num_of_start; // BFGS�̏����_�̐�

extern       double mean;		// GP�̎��O���z
extern const double sigma_thre; // u=0�Ƃ���臒l

extern double maxf;			// �f�[�^�Z�b�g�̍ő�l

extern MatrixXd D_q;		// �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;			// �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��

extern MatrixXd** A; //���͍s��
extern MatrixXd*  B; //i�ɂ��Ă�sum
extern MatrixXd*  C; //j�ɂ��Ă�sum

extern const double beta; //����z�@�̃X�e�b�v�T�C�Y(alp0 := norm / beta)

extern const int    ite_sev; //����z�@�̔�����
extern const double eps_sev; //����z�@�̃�

extern const int    ite_local; //�Ǐ��œK�l�֐��̔�����
extern const double eps_local; //�Ǐ��œK�l�֐��̃�

extern const int    ite_relax; //�ɘa���̔�����
extern const double eps_relax; //�ɘa���̃�

extern const int    ite_bc;	//���}����@�̔�����
extern const double eps_bc; //���}����@�̃�

extern int t_find; // BO�ɂ����čő�l��T�������Ƃ��̔�����
extern int k_find; // BC�ɂ����čő�l��T�������Ƃ��̔�����

extern const int  select;		// �ړI�֐��̐؂�ւ�
extern const bool bfgs_or_rand;	// argmax�̐؂�ւ��iBFGS�̃����_���j