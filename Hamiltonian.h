#ifndef HAMILTONIAN_KH
#define HAMILTONIAN_KH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <complex>
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

using namespace std;
using namespace arma;

typedef unsigned long long int ull_int;
typedef std::complex<double> cpx;
typedef std::vector<double> vect;

#define out std::cout << std::setprecision(16) << std::fixed

class HamiltonianKH {
public:
	int L; //chain length
	double t; //electron hopping
	double U; // electron repulsion
	double K; //exchange integral
	double J_H; // electron-spin interaction
	int Sz; // total spin for spin-sector selection
	int num_of_electrons; //number of electrons

	vector<ull_int> mapping; //generates the mapping of the base_vector number to the number in the block
	ull_int N; //number of states
    mat H;
	sp_mat H_sparse;
    mat eigenvectors;
	vec ground_state;
    vec eigenvalues;

	sp_mat Sz_tot2; // square of total spin: Sz_tot^2 = sum_ij Sz^i Sz^j
	vec chi_0; // static susceptibility

public:
	HamiltonianKH(); // default Constructor
	HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz); //Constructor for subblock Hamiltonian
	~HamiltonianKH();

	void Hamiltonian();
	void Hamiltonian_sparse();
	void Diagonalization();

	void generate_mapping();
    void CreateMappingElem(int &bSz, int &fSz, int& N_e, ull_int& j, ull_int& idx);
    void setHamiltonianElem(ull_int k, double value, std::vector<int> temp);
	void setHamiltonianElem_sparse(ull_int k, double value, std::vector<int> temp);
	void printEnergy(double Ef);
	void show_ground_state();
	void total_spin_squared();

// Functions not usable by Lanczos
	vec Heat_Capacity();
	vec Total_Density_of_states(std::vector<double> omega_vec);
	double spin_correlation_element(int site_i, int site_j, vec wavefunction);
};

//----------------------------------------------------------------------------------------------
//--------------------------------------------------TOOLS---------------------------------------
//----------------------------------------------------------------------------------------------
ull_int binary_search(std::vector<ull_int> arr, int l_point, int r_point, int element);

vector<int> int_to_binary(ull_int idx, int L); //converges int to binary code of length N
int binary_to_int(vector<int> vec); //converges vector with binary code to decimal system
double FermiLevel(int L, int N_e, double t, double K, double U, double J_H);
std::vector<double> prepareOmegaVec(double omega_min, double omega_max, double dOmega);
void printDOS(vec resultDOS, double U, double N_e, int L, std::vector<double> omega_vec, double maximum, double E_fermi);
void print_Cv(vec Cv, double U, double N_e, int L);
void print_Cv_Lanczos(vec Cv, double U, double N_e, int L, int M, int random_steps);
void print_chi(vec chi, double U, double N_e, int L);


void Main_U(int L, int N_e, double t);
void Main_Jh(int L, int N_e, double t, double K, double U);
void Main_Cv(int L, int N_e, double t, double K, double U, double J_H);
void Main_Cv_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps);
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H);


//----------------------------------------------------------------------------------------------
//--------------------------------------------------LANCZOS-------------------------------------
class Lanczos : public HamiltonianKH {
public:
	int lanczos_steps; // number of lanczos steps

	mat Krylov_space; // Krylov
	mat H_L; // lanczos hamiltonian
	vec randVec_inKrylovSpace; // overlap of random vector and eigenvectors

	Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps);
	Lanczos();
	~Lanczos();

	vec Create_Random_vec();
	void Build_Lanczos_Hamil(vec initial_vec);
	void Build_Lanczos_Hamil_wKrylovSpace(vec initial_vec);
	void Lanczos_GroundState();

	void Lanczos_Diagonalization();
	vec Hamil_vector_multiply(vec initial_vec);

	vec RandVec_times_KrylovTranspose(vec randvec);
	double Cv_kernel(double T, int random_steps);
	vec Heat_Capacity_Lanczos(int random_steps);
};


#endif

