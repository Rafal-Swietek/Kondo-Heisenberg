#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <algorithm>
// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
//#define ARMA_NO_DEBUG // only after checking if no out_of_bounds, or memory outusage
//-------
#include <armadillo>
#include <cassert> // assert terminates program
#include <omp.h>
#include <thread>
#include <time.h>
#include <utility> // auto, etc. 
#include <memory> // smart ptr

using namespace std;
using namespace arma;

typedef unsigned long long int u64;
typedef std::complex<double> cpx;
typedef std::unique_ptr<std::vector<u64>> my_uniq_ptr;

// User makros
#define im cpx(0.0,1.0)
#define out std::cout << std::setprecision(16) << std::fixed
#define num_of_threads 24

#define memory_over_performance false // optimized by size --true-- (memory usage shortage) or performance --false--
#define show_system_size_parameters true // this parameter defines whether to print data such as system size for each object conctructor call
#define calculate_stand_dev false // enables the calculation of the standard deviation in the FTLM randomized method
#define calculate_Sq true // calculate static structure factor?
#define use_reorthonormalization true // enables in lanczos procedure full reorthogonalization - needs full krylov_space access


extern double pi;
extern double T;
extern double dT;
extern double T_end;
extern double domega;
extern double eta;
extern double E0;
constexpr double W = 2.1;
constexpr auto tolerance = 1e-2;


//----------------------------------------------------------------------------------------------
//--------------------------------------------------TOOLS---------------------------------------
//----------------------------------------------------------------------------------------------
inline std::vector<double> prepare_parameterVec(double _min, double _max, double step) {
	std::vector<double> parameter_vec;
	double temp = _min;
	while (temp <= _max) {
		parameter_vec.emplace_back(temp);
		temp += step;
	}
	return parameter_vec;
}
inline u64 binary_search(my_uniq_ptr& arr, int l_point, int r_point, u64 element) {
	if (r_point >= l_point) {
		int middle = l_point + (r_point - l_point) / 2;
		if (arr->at(middle) == element) return middle;
		else if (arr->at(middle) > element) return binary_search(arr, l_point, middle - 1, element);
		else return binary_search(arr, middle + 1, r_point, element);
	}
	//u64 j = findElement(arr, element);
	out << "Element " << element << " not present in the array" << endl;
	assert(false);
	return -1;
}
inline void int_to_binary(u64 idx, std::vector<int>& vec) {
	u64 temp = idx;
	for (int k = 0; k < vec.size(); k++) {
		vec[vec.size() - 1 - k] = static_cast<int>(temp % 8);
		temp = static_cast<u64>((double)temp / 8.);
	}
}
inline u64 binary_to_int(vector<int>& vec) {
	u64 val = 0;
	for (int k = 0; k < vec.size(); k++) {
		val += vec[vec.size() - 1 - k] * (u64)std::pow(8, k);
	}
	return val;
}
inline vec Create_Random_vec(unsigned int N) {
	vec random_vec(N, fill::zeros);
	double norm = 0;
	for (u64 j = 0; j < N; j++) {
		random_vec(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
		norm += random_vec(j) * random_vec(j);
	}
	return random_vec / norm;
}

