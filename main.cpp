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
		//【tはこの時点におけるデータセットのサイズ】//

		//uのテスト
		if (t > 0 && t%5 == 0){
			cout << "t=" << t << endl;
			cout << "原点 : " << u(VectorXd::Zero(d)) << endl;
			cout << "ズレ : " << u(VectorXd::Constant(d, 0.5)) << endl;
		}


		//次のサンプル点の決定
		x_next = argmax_u();

		
		//データセットの更新
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		maxf = max(f(t), maxf);
		//【この時点でデータセットのサイズはt+1】//
		

		//KとKinvの更新
		updateK(x_next);
		//【updateKが行われた後，Kのサイズはt+1】//
	}

	//MatrixXd A = K*Kinv;
	//cout << A(T-1, T-2) << endl;

	//VectorXd v = VectorXd::Random(3);
	//cout << v << endl;
	//cout << v.head(2) << endl;

	return 0;
}