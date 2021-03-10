#ifndef LANCZOS
#define LANCZOS
#include "Hamiltonian.h"

class Lanczos : public HamiltonianKH {
private:
	sp_mat H; // here the hamil matrix is sparse
public:
	int lanczos_steps; // number of lanczos steps

	mat H_L; // lanczos hamiltonian
	vec randVec_inKrylovSpace; // overlap of random vector and eigenvectors
	vec partition_function; // partition function
	mat Krylov_space;

	Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps);
	Lanczos() = default;
	virtual ~Lanczos() = default;

	virtual void update_parameters(double t, double U, double K, double J_H, double Sz) override;
	virtual void Hamiltonian() override;
	virtual void setHamiltonianElem(u64& k, double value, std::vector<int>&& temp) override;
	virtual void Diagonalization() override;

	void Build_Lanczos_Hamil_wKrylov(vec& initial_vec);
	void Build_Lanczos_Hamil(vec& initial_vec);
	void Lanczos_convergence(vec& initial_vec);

	void Lanczos_GroundState(vec& initial_vec);
	void Hamil_vector_multiply_kernel(u64 start, u64 stop, vec& initial_vec, vec& result_vec_threaded);
	void Hamil_vector_multiply(vec& initial_vec, vec& result_vec);

	vec thermal_average_lanczos(vec&& quantity, int& random_steps);
	vec Heat_Capacity_Lanczos(int random_steps);
	vec entropy(int random_steps);
	vec static_spin_susceptibility(int random_steps);
	vec Sq_lanczos(int random_steps, double T);
	vec Sq_T0(int random_steps); //  zero-temperature static structure factor

	void SSF_T0(); // dynamical spin structure factor
};


#endif
