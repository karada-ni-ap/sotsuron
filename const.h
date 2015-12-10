#include <Eigen/Dense>
using namespace Eigen;

extern		 int t; // BOの反復回数
extern const int T; // BOの反復回数の上限

extern const int d; // xの次元
extern const int m; // yの次元
extern const int n; // 対称行列のサイズ

extern const int mean; // GPの事前分布

extern double maxf; // データセットの最大値

extern MatrixXd D_q; // クエリxのデータセット
extern VectorXd f;   // 値f(x)のデータセット

extern MatrixXd K;    // 共分散行列
extern MatrixXd Kinv; // 精度行列