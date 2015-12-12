#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "debug.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	VectorXd x_next = VectorXd::Zero(d);

	for (t = 0; t < T; t++){
		//�yt�͂��̎��_�ɂ�����f�[�^�Z�b�g�̃T�C�Y�z//
		debug_inside();

		//���̃T���v���_(t+1�_��)�̌���
		x_next = argmax_u();
		
		//�f�[�^�Z�b�g�̍X�V
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		maxf = max(f(t), maxf);
		//�y���̎��_�Ńf�[�^�Z�b�g�̃T�C�Y��t+1�z//


		//K��Kinv�̍X�V
		update_K(x_next);
		//�yupdate_K���s��ꂽ��CK�̃T�C�Y��t+1�z//
	}

	debug_last();

	return 0;
}