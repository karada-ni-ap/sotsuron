#include <fstream>
#include "experiment.h"

void vs_BO(int d_, int m_, int n_, int num){
	d = d_;
	m = m_;
	n = n_;

	clock_t* sample_time_BO = new clock_t[T];
	clock_t* sample_time_BC = new clock_t[ite_bc];

	double* sample_val_BO = new double[T];
	double* sample_val_BC = new double[ite_bc];

	string bar = "_";
	string s = to_string(d) + bar + to_string(m) + bar + to_string(n) + ".csv";

	ofstream out;
	out.open(s, ios::trunc);

	out << d << ", " << m << ", " << n << endl;
	out << T << ", " << ite_bc << endl;

	for (int num_ = 1; num_ <= num; num_++){
		initA();
		initialize_for_BO();

		BO(sample_time_BO, sample_val_BO);
		branch_and_cut(sample_time_BC, sample_val_BC);

		deleting_for_BO();

		//以下ファイル出力
		for (int i = 0; i < T; i++){
			out << sample_time_BO[i] << ", ";
		}
		out << -1 << endl;

		for (int i = 0; i < T; i++){
			out << sample_val_BO[i] << ", ";
		}
		out << -1 << endl;

		for (int i = 0; i < ite_bc; i++){
			out << sample_time_BC[i] << ", ";
		}
		out << -1 << endl;

		for (int i = 0; i < ite_bc; i++){
			out << sample_val_BC[i] << ", ";
		}
		out << -1 << endl;
	}

	out.close();
}

void vs_lsBO(int d_, int m_, int n_){
	d = d_;
	m = m_;
	n = m_;

	clock_t* sample_time_lsBO = new clock_t[T];
	clock_t* sample_time_BC = new clock_t[ite_bc];

	double* sample_val_lsBO = new double[T];
	double* sample_val_BC = new double[ite_bc];

	initA();
	initialize_for_BO();

	lsBO(sample_time_lsBO, sample_val_lsBO);
	branch_and_cut(sample_time_BC, sample_val_BC);

	deleting_for_BO();

	for (int i = 0; i < T; i++){
		cout << sample_time_lsBO[i] << "  :  " << sample_val_lsBO[i] << endl;
	}

	cout << "------------------------------" << endl;

	for (int i = 0; i < ite_bc; i++){
		cout << sample_time_BC[i] << "  :  " << sample_val_BC[i] << endl;;
	}
}