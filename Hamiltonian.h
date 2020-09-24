#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>
//#include <Windows.h>
//#include <mpi.h>
//#include <clocale>
//include <thread> //multithreading
//#include <mutex>
//#include <omp.h>
//#include <execution>
//#include <atomic> // atomic variable type allows multithreading without mutex locking

using namespace std;
using namespace arma;

class HamiltonianKH {
private:
	int L; //chain length
	double t; //electron hopping
	double U; // electron repulsion
	double K; //exchange integral
	double J_H; // electron-spin interaction

	int num_of_electrons; //number of electrons

	int N; //number of states - if using blocks its binomial(n,i) - i-th block, otherwise N=2^L
	vector<int> base_vector;
	vector<int> mapping; //generates the mapping of the base_vector number to the number in the block
	vector<int> mapping_inv;


	mat H;
	mat eigenvectors;
	vec eigenvalues;

	/*mat Krylov_space;
	vec Lanczos_eigenVal;
	mat Lanczos_eigenVec;;
	mat H_Lanczos;*/

public:
	HamiltonianKH(int L, double t, double U, double K, double J_H); //Constructor for total Hamiltonian
	HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H); //Constructor for subblock Hamiltonian
	// num_of_electrons =0,...,2*L
	~HamiltonianKH();

	mat get_hamil();
	vec get_energy();
	cx_vec coeff;

	void Hamiltonian();
	void Diagonalization();

	void Density_of_states(int N_e);


	// Tools
		void generate_mapping_subblock();
		void generate_mapping_total();
	//-----------------------------

};

vector<int> int_to_binary(int idx, int L); //converges int to binary code of length N
int binary_to_int(vector<int> vec); //converges vector with binary code to decimal system

long long int factorial(int n);
long long int Binomial(int n, int k); //Binomial coefficient n over k

void Main_DOS_U(int L, int N_e, double t);


