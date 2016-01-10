#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <time.h>

using namespace std;
using namespace Eigen;

extern const int Inf; //大きい値

extern		 int t;			// BOの反復回数
extern const int T;			// BOの反復回数の上限

extern int d; // xの次元
extern int m; // yの次元
extern int n; // 対称行列のサイズ

extern VectorXd Ux0;
extern VectorXd Lx0;
extern VectorXd Uy0;
extern VectorXd Ly0;

extern const double alp0;	// バックトラック法のαの初期値
extern const double eps_alp;// バックトラック法のαのbreak条件
extern const double rho;	// バックトラック法のρ
extern const double c1;		// バックトラック法のc1
extern const double c2;		// バックトラック法のc2

extern const double eps_y; // H公式のε

extern const int    ite_bfgs;	 // BFGSの反復回数
extern const double eps_bfgs;	 // BFGSの収束判定の閾値

extern const int    ite_sdm; // 最急降下法の反復回数
extern const double eps_sdm; // 最急降下法の収束判定の閾値

extern const int num_of_start; // BFGSの初期点の数

extern       double mean;		// GPの事前分布
extern const double sigma_thre; // u=0とする閾値
extern const double xi;			// トレードオフパラメータ

extern double maxf_BO;			// データセットの最大値（BO）
extern double maxf_lsBO;		// データセットの最大値（lsBO）

extern MatrixXd D_q;		// クエリxのデータセット
extern VectorXd f;			// 値f(x)のデータセット

extern MatrixXd K;    // 共分散行列
extern MatrixXd Kinv; // 精度行列

extern MatrixXd** A; //入力行列
extern MatrixXd*  B; //iについてのsum
extern MatrixXd*  C; //jについてのsum

extern const int    ite_sev; //劣勾配法の反復回数
extern const double eps_sev; //劣勾配法のε

extern const int    ite_local; //局所最適値関数の反復回数
extern const double eps_local; //局所最適値関数のε

extern const int    ite_relax; //緩和問題の反復回数
extern const double eps_relax; //緩和問題のε

extern const int    ite_bc;	//分枝限定法の反復回数
extern const double eps_bc; //分枝限定法のε

extern       bool SDMorBFGS;	// argmaxの切り替え（SDM⇔BFGS）