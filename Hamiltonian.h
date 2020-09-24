#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

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

using namespace std;
using namespace arma;

#define out std::cout << std::setprecision(8) << std::fixed

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

public:
    mat H;
    mat eigenvectors;
    vec eigenvalues;
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
        void mappingAddIfNeeded(int &bSz, int &fSz, int& N_e, int &j, int &idx);
        void setHamiltonianElem(int k, double value, std::vector<int> temp);
	//-----------------------------

};

vector<int> int_to_binary(int idx, int L); //converges int to binary code of length N
int binary_to_int(vector<int> vec); //converges vector with binary code to decimal system
void Main_DOS_U(int L, int N_e, double t);

#endif

