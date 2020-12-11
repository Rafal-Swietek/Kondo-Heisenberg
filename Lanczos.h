#ifndef LANCZOS
#define LANCZOS
#include "Hamiltonian.h"

class Lanczos : public HamiltonianKH {
private:
	sp_mat H; // here the hamil matrix is sparse
public:
	int lanczos_steps; // number of lanczos steps
	double Z_constT; // partition function for given temperature;

	mat H_L; // lanczos hamiltonian
	vec randVec_inKrylovSpace; // overlap of random vector and eigenvectors
	vec partition_function; // partition function
	mat Krylov_space;

	Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps);
	Lanczos();
	~Lanczos();

	virtual void update_parameters(double t, double U, double K, double J_H, double Sz) override;
	virtual void Hamiltonian() override;
	virtual void setHamiltonianElem(ull_int& k, double value, std::vector<int>&& temp) override;
	virtual void Diagonalization() override;

	void Build_Lanczos_Hamil_wKrylov(vec& initial_vec);
	void Build_Lanczos_Hamil(vec& initial_vec);
	void Lanczos_convergence(vec& initial_vec);

	void Lanczos_GroundState();
	//void Lanczos_Diagonalization();
	void Hamil_vector_multiply_kernel(ull_int start, ull_int stop, vec& initial_vec, vec& result_vec_threaded);
	void Hamil_vector_multiply(vec& initial_vec, vec& result_vec);
	vec LDOS(int site, double T, std::vector<double>&& omega_vec);

	vec thermal_average_lanczos(vec&& quantity, int& random_steps);
	vec Heat_Capacity_Lanczos(int random_steps);
	vec static_spin_susceptibility(int random_steps);
	vec Sq_lanczos(int random_steps, double T);
};


#endif
