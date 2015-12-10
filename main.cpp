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

	for (t=0; t < T; t++){
		//mu,sigmaの確認
		if (t > 0){
			cout << Kinv.topLeftCorner(t, t) << endl;
			cout << "mu : " << mu(VectorXd::Constant(d,0.5)) << endl;
			cout << "---------------------------------------" << endl;
		}
		//【tは現在のデータセットのサイズ】//
		
		//次のサンプル点の決定
		if (t == 0)
			x_next = VectorXd::Random(d);
		else
			x_next = argmax_u();
		
		//データセットの更新
		D_q.col(t) = x_next;
		f(t) = obj(x_next);
		//【この時点でデータセットのサイズはt+1】//
		maxf = max(f(t), maxf);
		updateK(x_next);
		//【updateKが行われたとき，Kのサイズはt+1】//

	}

	//MatrixXd A = K*Kinv;
	//cout << A(T-1, T-2) << endl;

	//VectorXd v = VectorXd::Random(3);
	//cout << v << endl;
	//cout << v.head(2) << endl;

	return 0;
}