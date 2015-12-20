#include <Eigen/Dense>
using namespace Eigen;

extern const int Inf; //大きい値

extern		 int t;			// BOの反復回数
extern       int t_find;	// maxfを探索したときの反復回数
extern const int T;			// BOの反復回数の上限

extern const int d; // xの次元
extern const int m; // yの次元
extern const int n; // 対称行列のサイズ

extern const VectorXd Ux0;
extern const VectorXd Lx0;
extern const VectorXd Uy0;
extern const VectorXd Ly0;

extern const double Alp0;		 // バックトラック法のαの初期値
extern const double rho;		 // バックトラック法のρ
extern const double c1_bfgs;	 // バックトラック法のc1

extern const double eps_bfgs;	 // BFGSの収束判定の閾値
extern const int    ite_bfgs;	 // BFGSの反復回数
extern const int	num_bfgs;	 // BFGSの初期点の数

extern const int    mean;		// GPの事前分布
extern const double sigma_thre; // u=0とする閾値

extern double maxf;			// データセットの最大値

extern MatrixXd D_q;		// クエリxのデータセット
extern VectorXd f;			// 値f(x)のデータセット

extern MatrixXd K;    // 共分散行列
extern MatrixXd Kinv; // 精度行列

extern MatrixXd** A; //入力行列
extern MatrixXd*  B; //iについてのsum
extern MatrixXd*  C; //jについてのsum

extern const double beta; //劣勾配法のステップサイズの倍率

extern const int    ite_sev; //劣勾配法の反復回数
extern const double eps_sev; //劣勾配法のε

extern const int    ite_local; //局所最適値関数の反復回数
extern const double eps_local; //局所最適値関数のε

extern const int    ite_relax; //緩和問題の反復回数
extern const double eps_relax; //緩和問題のε
extern const double alpW0;     //緩和問題のWの初期ステップサイズ

extern const int    ite_bc;	//分枝限定法の反復回数
extern const double eps_bc; //分枝限定法のε

extern const int  select;		// 目的関数の切り替え
extern const bool bfgs_or_rand;	// argmaxの切り替え（BFGS⇔ランダム）