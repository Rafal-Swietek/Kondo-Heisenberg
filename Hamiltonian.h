#ifndef HAMILTONIAN_KH
#define HAMILTONIAN_KH
#include "headers.h"
//#include "sparse_mat.h"

class HamiltonianKH {
private:
	mat H;
public: // non-use in lanczos class
	virtual vec Total_Density_of_states(std::vector<double>&& omega_vec) final;

public:
	int L; //chain length
	double t; //electron hopping
	double U; // electron repulsion
	double K; //exchange integral
	double J_H; // electron-spin interaction
	int Sz; // total spin for spin-sector selection
	int num_of_electrons; //number of electrons

	my_uniq_ptr mapping;//generates the mapping of the base_vector number to the number in the block
	ull_int N; //number of states
    mat eigenvectors;
    vec eigenvalues;
	vec ground_state;

	HamiltonianKH(); // default Constructor
	HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz); //Constructor for subblock Hamiltonian
	~HamiltonianKH();

	// will be overridden in Lanczos class
	virtual void update_parameters(double t, double U, double K, double J_H, double Sz);
	virtual void Hamiltonian();
	virtual void setHamiltonianElem(ull_int& k, double value, std::vector<int>&& temp);
	virtual void Diagonalization();
	//--

	void generate_mapping();
	void mapping_kernel(ull_int start, ull_int stop, my_uniq_ptr& map_threaded, int _id);
	
	void printEnergy(double Ef);
	void show_ground_state();

	virtual vec static_structure_factor(double T);
	double partition_function(double T);
	mat correlation_matrix();
	void print_base_vector(std::vector<int>& base_vector, std::ofstream& out_str);
};

#endif

