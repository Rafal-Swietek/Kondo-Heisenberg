#pragma once

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

using namespace std;
using namespace arma;

typedef unsigned long long int ull_int;
typedef std::complex<double> cpx;
typedef std::unique_ptr<std::vector<ull_int>> my_uniq_ptr;

// User makros
#define im cpx(0.0,1.0)
#define out std::cout << std::setprecision(16) << std::fixed
#define num_of_threads 16

#define memory_over_performance false // optimized by size --true-- (memory usage shortage) or performance --false--
#define show_system_size_parameters true // this parameter defines whether to print data such as system size for each object conctructor call
#define calculate_stand_dev false // enables the calculation of the standard deviation in the FTLM randomized method

#undef use_reorthonormalization // enables in lanczos procedure full reorthogonalization - needs full krylov_space access
#if defined(calculate_Sq)
	#define use_reorthonormalization true
#else
	#define use_reorthonormalization false
#endif
#define calculate_Sq false // calculate static structure factor?


extern double pi;
extern double T;
extern double dT;
extern double T_end;
extern double domega;
extern double eta;