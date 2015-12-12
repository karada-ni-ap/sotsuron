#include <Eigen/Dense>
using namespace Eigen;

extern		 int t; // BOの反復回数
extern const int T; // BOの反復回数の上限

extern const int d; // xの次元
extern const int m; // yの次元
extern const int n; // 対称行列のサイズ

extern const double Ub;  // 超立方体の上界
extern const double Lb;  // 超立方体の下界

extern const double Alp;		 // BFGSのステップサイズ
//extern const double conv_thre; // BFGSの収束判定の閾値
extern const int    ite_bfgs;	 // BFGSの反復回数
extern const int	num_bfgs;	 // BFGSの初期点の数

extern const int    mean;		// GPの事前分布
extern const double sigma_thre; // u=0とする閾値

extern double maxf; // データセットの最大値

extern MatrixXd D_q; // クエリxのデータセット
extern VectorXd f;   // 値f(x)のデータセット

extern MatrixXd K;    // 共分散行列
extern MatrixXd Kinv; // 精度行列

extern const bool bfgs_or_rand; // Debug用