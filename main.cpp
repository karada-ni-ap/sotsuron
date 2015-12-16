#include <iostream>
#include <Eigen/Dense>
#include "const.h"
#include "BO.h"
#include "myfunc.h"
#include "obj.h"
#include "argmax.h"
#include "debug.h"
#include "sev.h"

using namespace std;
using namespace Eigen;

int main(void)
{	
	initA();

	for (t = 0; t < T; t++){
		//【tはこの時点におけるデータセットのサイズ】//

		//次のサンプル点(t+1点目)の決定
		x_next = argmax_u();
		debug_inside();
		
		//データセットの更新
		D_q.col(t) = x_next;
		f(t) = obj(x_next, select);

		if (f(t)>maxf){
			maxf = f(t);
			x_opt = x_next;
		}
		//【この時点でデータセットのサイズはt+1】//


		//KとKinvの更新
		update_K(x_next);
		//【update_Kが行われた後，Kのサイズはt+1】//
	}

	//このときt=Tで，データセットのサイズもT．
	debug_last();

	return 0;
}