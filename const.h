#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <time.h>

using namespace std;
using namespace Eigen;

extern const int Inf; //�傫���l

extern		 int t;			// BO�̔�����
extern const int T;			// BO�̔����񐔂̏��

extern int d; // x�̎���
extern int m; // y�̎���
extern int n; // �Ώ̍s��̃T�C�Y

extern VectorXd Ux0;
extern VectorXd Lx0;
extern VectorXd Uy0;
extern VectorXd Ly0;

extern const double alp0;	// �o�b�N�g���b�N�@�̃��̏����l
extern const double eps_alp;// �o�b�N�g���b�N�@�̃���break����
extern const double rho;	// �o�b�N�g���b�N�@�̃�
extern const double c1;		// �o�b�N�g���b�N�@��c1
extern const double c2;		// �o�b�N�g���b�N�@��c2

extern const double eps_y; // H�����̃�

extern const int    ite_bfgs;	 // BFGS�̔�����
extern const double eps_bfgs;	 // BFGS�̎��������臒l

extern const int    ite_sdm; // �ŋ}�~���@�̔�����
extern const double eps_sdm; // �ŋ}�~���@�̎��������臒l

extern const int num_of_start; // BFGS�̏����_�̐�

extern       double mean;		// GP�̎��O���z
extern const double sigma_thre; // u=0�Ƃ���臒l
extern const double xi;			// �g���[�h�I�t�p�����[�^

extern double maxf_BO;			// �f�[�^�Z�b�g�̍ő�l�iBO�j
extern double maxf_lsBO;		// �f�[�^�Z�b�g�̍ő�l�ilsBO�j

extern MatrixXd D_q;		// �N�G��x�̃f�[�^�Z�b�g
extern VectorXd f;			// �lf(x)�̃f�[�^�Z�b�g

extern MatrixXd K;    // �����U�s��
extern MatrixXd Kinv; // ���x�s��

extern MatrixXd** A; //���͍s��
extern MatrixXd*  B; //i�ɂ��Ă�sum
extern MatrixXd*  C; //j�ɂ��Ă�sum

extern const int    ite_sev; //����z�@�̔�����
extern const double eps_sev; //����z�@�̃�

extern const int    ite_local; //�Ǐ��œK�l�֐��̔�����
extern const double eps_local; //�Ǐ��œK�l�֐��̃�

extern const int    ite_relax; //�ɘa���̔�����
extern const double eps_relax; //�ɘa���̃�

extern const int    ite_bc;	//���}����@�̔�����
extern const double eps_bc; //���}����@�̃�

extern       bool SDMorBFGS;	// argmax�̐؂�ւ��iSDM��BFGS�j