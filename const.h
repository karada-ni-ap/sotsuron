#include <Eigen/Dense>
using namespace Eigen;

extern		 int t; // BOの反復回数
extern const int T; // BOの反復回数の上限

extern const int d; // xの次元
extern const int m; // yの次元
extern const int n; // 対称行列のサイズ

extern const double Ub;  // 超立方体の上界
extern const double Lb;  // 超立方体の下界

extern const double Alp0;		 // バックトップ法のαの初期値
extern const double rho;		 // バックトラック法のρ
extern const double c1_bfgs;	 // バックトラック法のc1

extern const double eps_bfgs;	 // BFGSの収束判定の閾値
extern const int    ite_bfgs;	 // BFGSの反復回数
extern const int	num_bfgs;	 // BFGSの初期点の数

extern const int    mean;		// GPの事前分布
extern const double sigma_thre; // u=0とする閾値

extern double maxf;			// データセットの最大値

extern VectorXd x_next;		// 次のサンプル点
extern VectorXd x_opt;		// 最適解

extern MatrixXd D_q;		// クエリxのデータセット
extern VectorXd f;			// 値f(x)のデータセット

extern MatrixXd K;    // 共分散行列
extern MatrixXd Kinv; // 精度行列

extern MatrixXd** A; //入力行列
extern MatrixXd*  B; //iについてのsum

extern const int    ite_subgrad; //最急降下法の反復回数
extern const int    Alp_subgrad; //最急降下法のステップサイズの初期値
extern const double rho_subgrad; //最急降下法のステップサイズの減少率

extern const int select;		// 目的関数の切り替え
extern const bool bfgs_or_rand;	// argmaxの切り替え（BFGS⇔ランダム）