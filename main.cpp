#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	VectorXd x_next = VectorXd::Zero(d);

	for (t = 0; t < T; t++){
		//�yt�͂��̎��_�ɂ�����f�[�^�Z�b�g�̃T�C�Y�z//

		//u�̃e�X�g
		if (t > 0 && t%5 == 0){
			cout << "t=" << t << endl;
			cout << "���_ : " << u(VectorXd::Zero(d)) << endl;
			cout << "�Y�� : " << u(VectorXd::Constant(d, 0.5)) << endl;
		}


		//���̃T���v���_�̌���
		x_next = argmax_u();

		
		//�f�[�^�Z�b�g�̍X�V
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		maxf = max(f(t), maxf);
		//�y���̎��_�Ńf�[�^�Z�b�g�̃T�C�Y��t+1�z//
		

		//K��Kinv�̍X�V
		updateK(x_next);
		//�yupdateK���s��ꂽ��CK�̃T�C�Y��t+1�z//
	}

	//MatrixXd A = K*Kinv;
	//cout << A(T-1, T-2) << endl;

	//VectorXd v = VectorXd::Random(3);
	//cout << v << endl;
	//cout << v.head(2) << endl;

	return 0;
}