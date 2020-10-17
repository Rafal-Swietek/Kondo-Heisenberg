#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <complex>
#include <armadillo>
#include <cassert> // assert terminates program
#include <omp.h>

using namespace std;
using namespace arma;

#define out std::cout << std::setprecision(16) << std::fixed

class HamiltonianKH {
private:
	int L; //chain length
	double t; //electron hopping
	double U; // electron repulsion
	double K; //exchange integral
	double J_H; // electron-spin interaction
	int Sz; // total spin for spin-sector selection
	int num_of_electrons; //number of electrons
	int N; //number of states

	vector<int> mapping; //generates the mapping of the base_vector number to the number in the block

public:
    mat H;
    mat eigenvectors;
    vec eigenvalues;

public:
	HamiltonianKH(); // default Constructor
	HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz); //Constructor for subblock Hamiltonian
	~HamiltonianKH();

	void Hamiltonian();
	void Diagonalization();
	vec Heat_Capacity();
	vec Total_Density_of_states(std::vector<double> omega_vec);

	void generate_mapping();
    void CreateMappingElem(int &bSz, int &fSz, int& N_e, int &j, int &idx);
    void setHamiltonianElem(int k, double value, std::vector<int> temp);

	void printEnergy();
	double spin_correlation_element(int site_i, int site_j, vec wavefunction);
};

//----------------------------------------------------------------------------------------------
//--------------------------------------------------TOOLS---------------------------------------
//----------------------------------------------------------------------------------------------
int binary_search(std::vector<int> arr, int l_point, int r_point, int element);

vector<int> int_to_binary(int idx, int L); //converges int to binary code of length N
int binary_to_int(vector<int> vec); //converges vector with binary code to decimal system
double FermiLevel(int L, int N_e, double t, double K, double U, double J_H);
std::vector<double> prepareOmegaVec(double omega_min, double omega_max, double dOmega);
void printDOS(vec resultDOS, double U, double N_e, int L, std::vector<double> omega_vec, double maximum, double E_fermi);
void print_Cv(vec Cv, double U, double N_e, int L);
void print_Cv_Lanczos(vec Cv, double U, double N_e, int L, int M);


void Main_U(int L, int N_e, double t);
void Main_Jh(int L, int N_e, double t, double K, double U);
void Main_Cv(int L, int N_e, double t, double K, double U, double J_H);
void Main_Cv_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M);
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H);


//----------------------------------------------------------------------------------------------
//--------------------------------------------------LANCZOS-------------------------------------
class Lanczos :public HamiltonianKH {
private:
	int L; //chain length
	double t; //electron hopping
	double U; // electron repulsion
	double K; //exchange integral
	double J_H; // electron-spin interaction
	int Sz; // total spin for spin-sector selection
	int num_of_electrons; //number of electrons
	int N; //number of states
	int lanczos_steps; // number of lanczos steps

	vector<int> mapping; //generates the mapping of the base_vector number to the number in the block

public:
	mat Krylov_space; // Krylov
	mat H_L; // lanczos hamiltonian
	vec eigenVal_L; // lanczos eigenvalues
	mat eigenVec_L; // lanczos eigenvect5ors in Krylov subspace
	vec Lanczos_GS;
	vec randVec_inKrylovSpace; // overlap of random vector and eigenvectors

	Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps);
	Lanczos();
	~Lanczos();
	void generate_mapping();
	void CreateMappingElem(int& bSz, int& fSz, int& N_e, int& j, int& idx);
	void setHamiltonianElem(int k, double value, std::vector<int> temp);


	vec Create_Random_vec();
	void Build_Lanczos_Hamil(vec initial_vec);
	void Build_Lanczos_Hamil_wKrylovSpace(vec initial_vec);
	void Lanczos_GroundState();

	void Lanczos_Diagonalization();
	vec Hamil_vector_multiply(vec initial_vec);

	vec RandVec_times_KrylovTranspose(vec randvec);
	double Cv_kernel(double T);
	vec Heat_Capacity_Lanczos();
};


#endif

