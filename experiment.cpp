#include "experiment.h"

void vs_BO(int d_, int m_, int n_){
	d = d_;
	m = m_;
	n = m_;

	clock_t* sample_time_BO = new clock_t[T];
	clock_t* sample_time_BC = new clock_t[T];

	double* sample_val_BO = new double[T];
	double* sample_val_BC = new double[T];

	initA();
	initialize_for_BO();

	BO(sample_time_BO, sample_val_BO);
	//branch_and_cut(sample_time_BC, sample_val_BC);

	deleting_for_BO();

	for (int i = 0; i < T; i++){
		cout << sample_time_BO[i] << "  :  " << sample_val_BO[i] << endl;;
	}
}

void vs_lsBO(int d_, int m_, int n_){
	d = d_;
	m = m_;
	n = m_;

	clock_t* sample_time_lsBO = new clock_t[T];
	clock_t* sample_time_BC = new clock_t[T];

	double* sample_val_lsBO = new double[T];
	double* sample_val_BC = new double[T];

	initA();
	initialize_for_BO();

	lsBO(sample_time_lsBO, sample_val_lsBO);
	//branch_and_cut(sample_time_BC, sample_val_BC);

	deleting_for_BO();

	for (int i = 0; i < T; i++){
		cout << sample_time_lsBO[i] << "  :  " << sample_val_lsBO[i] << endl;;
	}
}