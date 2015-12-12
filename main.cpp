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
		//【tはこの時点におけるデータセットのサイズ】//
		debug_inside();

		//次のサンプル点(t+1点目)の決定
		x_next = argmax_u();
		
		//データセットの更新
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		maxf = max(f(t), maxf);
		//【この時点でデータセットのサイズはt+1】//


		//KとKinvの更新
		update_K(x_next);
		//【update_Kが行われた後，Kのサイズはt+1】//
	}

	debug_last();

	return 0;
}