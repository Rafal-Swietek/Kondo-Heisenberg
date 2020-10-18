#include "Hamiltonian.h"

#define im std::complex<double>(0.0,1.0)
#define M_PI 3.14159265358979323846

double pi = M_PI;

typedef std::complex<double> cpx;
typedef std::vector<double> vect;

double dT = 0.002;
double T_end = 3.0;
double domega = 0.01;

//Destructor
HamiltonianKH::~HamiltonianKH() {}
// Constructor for Hamiltonian separated to blocks
HamiltonianKH::HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz) {
	this->L = L; //number of sites
	this->num_of_electrons = num_of_electrons; //number of electrons in lower orbital
	this->t = t; this->U = U; 
	this->K = K;
	this->J_H = J_H;
    this->Sz = Sz;

	this->mapping = vector<unsigned long int>();
	generate_mapping();
	this->N = mapping.size();
}
HamiltonianKH::HamiltonianKH() {}
//-------------------------

void HamiltonianKH::setHamiltonianElem(int k, double value, std::vector<int> temp){
    int idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
        H(idx, k) += value;
        H(k, idx) += value;
}
void HamiltonianKH::Hamiltonian() {
    this->H = sp_mat(N, N); //hamiltonian

    int s_i, s_j; //i=j, j=j+1
	bool PBC = 0; //allows periodic boundary conditions if =1
	int next_j;
    for (int k = 0; k < N; k++){
		std::vector<int> base_vector = int_to_binary(mapping[k], L);
		vector<int> temp(base_vector);
		for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;
                // Diagonal spin part
                if (base_vector[j] < 4) s_i = 1;
                else s_i = 0;

                // PBC = i+1 : (L-1)+1 = 0
                if (base_vector[next_j] < 4) s_j = 1;
                else s_j = 0;
                H(k, k) += K * (s_i - 0.5) * (s_j - 0.5);

                //Kinetic spin part: S+ S-
                temp = base_vector;
                if (s_i == 0 && s_j == 1) { // S_i^+ S_i+1^-
                    temp[j] = base_vector[j] - 4; //spin filp
                    temp[next_j] = base_vector[next_j] + 4;
                    setHamiltonianElem(k, K / 2., temp);
                }
                // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                    //spin up
                temp = base_vector;
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                    temp[next_j] -= 1; // anihilate spin-up electron
                    temp[j] += 1; // create spin-up electron
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        setHamiltonianElem(k, -t, temp);
                    }
                    else  setHamiltonianElem(k, +t, temp);
                }
                //spin down
                temp = base_vector;
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        temp[next_j] -= 2; // anihilate spin-down electron
                        temp[j] += 2; // create spin-down electron
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            setHamiltonianElem(k, -t, temp);
                        }
                        else  setHamiltonianElem(k, +t, temp);
                    }
                }
            }
			// electron repulsion
				if (base_vector[j] == 7 || base_vector[j] == 3) 
					H(k, k) += U;
			// electron-localised spin interaction ( interorbital electronic spin interaction)
				temp = base_vector;
				if (base_vector[j] == 5) {// S_i^+ s_i^-
					temp[j] = 2;
                    setHamiltonianElem(k,-J_H,temp);
                }
		    //Diagonal - z part
				if (base_vector[j] == 1 || base_vector[j] == 6)
					H(k, k) -= 2.0 * J_H * 0.25;
				if (base_vector[j] == 2 || base_vector[j] == 5)
					H(k, k) += 2.0 * J_H * 0.25;
        }
	}
}

//generates the vector, which maps the base_vector index to the index in given subblock
void HamiltonianKH::CreateMappingElem(int& bSz, int& fSz, int& N_e, int& j, unsigned long int& idx) {
    if ((bSz + fSz == Sz) && N_e == num_of_electrons) {
        mapping.push_back(j);
        idx++;
    }
}
std::tuple<int, int, int> calculateSpinElements(int L, int j) {
    int bSz = 0; //bosonic total spin - spin of upper orbital locked to n=1 filling
    int fSz = 0; //fermionic total spin
    int N_e = 0; // numer of electrons in given state
    vector<int> temp = int_to_binary(j, L);

    for (int k = 0; k < L; k++) {
        if (temp[k] < 4) bSz += 1;
        else bSz -= 1;
        if (temp[k] % 4 == 1) {
            fSz += 1;
            N_e += 1;
        }
        else if (temp[k] % 4 == 2) {
            fSz -= 1;
            N_e += 1;
        }
        else if (temp[k] % 4 == 3)
            N_e += 2;
    }

    return std::make_tuple(bSz, fSz, N_e);
}
void HamiltonianKH::generate_mapping() {
    unsigned long int idx = 0;
	for (int j = 0; j < std::pow(8, L); j++) {
        int bSz = 0, fSz = 0, N_e = 0;
        std::tie(bSz,fSz,N_e) = calculateSpinElements(L,j);
        CreateMappingElem(bSz,fSz, N_e,j,idx);
    }
    assert(idx>0 && "Not possible number of electrons - no. of states < 1");
}
//----------------------------------------------------

void HamiltonianKH::Diagonalization() {
	try {
		//arma::eig_sym(eigenvalues, eigenvectors, H);
        this->ground_state = eigenvectors.col(0);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
        out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
        	assert(false);
    }
    //out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
}

// Calculates the density of states using one-particle greens function
vec HamiltonianKH::Total_Density_of_states(std::vector<double> omega_vec) {
	double domega = 0.005;

    vec resultDOS(omega_vec.size());

    double maximum = 0;
//#pragma omp parallel for shared(omega_vec, resultDOS) num_threads(16)
    for (int w = 0; w < omega_vec.size(); w++) {
        double omega = omega_vec[w];
        double DOS = 0;
    //#pragma omp parallel for shared(omega_vec, resultDOS) reduction(+: DOS)
        for (int n = 0; n < N; n++)
                DOS += -1. / (double)L / pi * imag(1. / (omega + 2 * domega * 1i - eigenvalues(n)));

        if (DOS > maximum)
             maximum = DOS;

        resultDOS(w) = DOS;
    }
    //printDOS(resultDOS,U, num_of_electrons,L,omega_vec,maximum, E_F);
    return resultDOS;
}
vec HamiltonianKH::Heat_Capacity() {
    double energy_av; //average of energy E
    double energy2_av; //average of E^2
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->chi_0 = vec (static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    int k = 0;
    double T = dT;
    while (T <= T_end) {
        double Partition_Function = 0;
        energy_av = 0; energy2_av = 0;
        for (int j = 0; j < N; j++) {
            Partition_Function += std::exp(-(eigenvalues(j) - eigenvalues(0)) / T); //partition function(T)
            energy_av += eigenvalues(j) * std::exp(-(eigenvalues(j) - eigenvalues(0)) / T); //average energy(T)
            energy2_av += eigenvalues(j) * eigenvalues(j) * std::exp(-(eigenvalues(j) - eigenvalues(0)) / T); //average energy^2(T)
        }
        energy_av = energy_av / Partition_Function;
        energy2_av = energy2_av / Partition_Function;
        double heat_capacity = (energy2_av - energy_av * energy_av) / T / T / (L + 0.0);
        Cv(k) = heat_capacity;
        chi_0(k) = static_spin_susceptibility(T);
        T += dT; k++;
    }
    return Cv;
}
double HamiltonianKH::spin_correlation_element(int site_i, int site_j, vec wavefunction) {
    double result = 0;
#pragma omp parallel for reduction(+:result)
    for (int n = 0; n < N; n++) {
        vector<int> temp = int_to_binary(n, L);
        double fSz = 0, bSz = 0;

        if (temp[site_i] < 4) bSz += 0.5;
        else bSz -= 0.5;
        if (temp[site_i] % 4 == 1)
            fSz += 0.5;
        else if (temp[site_i] % 4 == 2)
            fSz -= 0.5;

        if (temp[site_j] < 4) bSz += 0.5;
        else bSz -= 0.5;
        if (temp[site_j] % 4 == 1)
            fSz += 0.5;
        else if (temp[site_j] % 4 == 2)
            fSz -= 0.5;

        result += wavefunction(n) * (fSz + bSz) * wavefunction(n);
    }
    return result;
}
double HamiltonianKH::static_spin_susceptibility(double T) {
    double Z = 0, chi_0 = 0;
    for (int n = 0; n < N; n++) {
        vector<int> temp = int_to_binary(n, L);
        for (int k = 0; k < L; k++) {
            for (int q = 0; q < L; q++) {
                double fSz1 = 0, fSz2 = 0, bSz1 = 0, bSz2 = 0;

                if (temp[k] < 4) bSz1 += 0.5;
                else bSz1 -= 0.5;
                if (temp[k] % 4 == 1)
                    fSz1 += 0.5;
                else if (temp[k] % 4 == 2)
                    fSz1 -= 0.5;

                if (temp[q] < 4) bSz2 += 0.5;
                else bSz2 -= 0.5;
                if (temp[q] % 4 == 1)
                    fSz2 += 0.5;
                else if (temp[q] % 4 == 2)
                    fSz2 -= 0.5;
                chi_0 += ground_state(n) * (fSz1 + bSz1) * (fSz2 + bSz2) * ground_state(n);
            }
        }
        Z += std::exp(-(eigenvalues(n) - eigenvalues(0)) / T);
    }
    return chi_0 / Z;
}

void HamiltonianKH::show_ground_state() {
    vec GS((int)std::pow(2, L), fill::zeros);
    for (int k = 0; k < N; k++) {
        std::vector<int> base_vector = int_to_binary(mapping[k], L);
        int val = 0;
        for (int j = 0; j < L; j++) 
            val += base_vector[base_vector.size() - 1 - j] / 4 * std::pow(2, j);
        GS(val) += ground_state(k);
    }
    GS = arma::abs(GS) / dot(GS, GS); //normalizing to be sure
    out << endl;
    double maximum = arma::max(GS);
    for (int k = 0; k < GS.size(); k++) {
        if (std::fabs(GS(k)) >= 0.5 * maximum) {
            vector<int> vec(L);
            int temp = k;
            for (int p = 0; p < L; p++) {
                vec[vec.size() - 1 - p] = temp % 2;
                temp = static_cast<int>((double)temp / 2.);
            }
            out << "Ground state:\t" << "|";
            for (int j = 0; j < L; j++) {
                out << vec[j];
            }
            out << ">\nwith probability\t p=" << GS(k) * GS(k) << endl << endl;
        }
    }
}

//----------------------------------------------------------------------------------------------
//--------------------------------------------------LANCZOS-------------------------------------
Lanczos::Lanczos() {}
Lanczos::Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps) {
    this->L = L; //number of sites
    this->num_of_electrons = num_of_electrons; //number of electrons in lower orbital
    this->t = t; this->U = U;
    this->K = K;
    this->J_H = J_H;
    this->Sz = Sz;
    this->lanczos_steps = lanczos_steps;

    this->mapping = vector<unsigned long int>();
    generate_mapping();
    this->N = mapping.size(); //out << "dim = " << N * sizeof(mapping[0]) << endl;
}
Lanczos::~Lanczos() {}

vec Lanczos::Create_Random_vec() {
    vec random_vec(N, fill::zeros);
    double norm = 0;
    for (int j = 0; j < N; j++) {
        random_vec(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
        norm += random_vec(j) * random_vec(j);
    }
    return random_vec / norm;
}
vec Lanczos::Hamil_vector_multiply(vec initial_vec) {
    vec result_vec(N, fill::zeros);
    for (int k = 0; k < N; k++) {
        std::vector<int> base_vector = int_to_binary(mapping[k], L);
        vector<int> temp(base_vector);
        int PBC = 0;
        int next_j, idx;
        int s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;

                if (base_vector[j] < 4) s_i = 1;
                else s_i = 0;
                if (base_vector[next_j] < 4) s_j = 1;
                else s_j = 0;

                result_vec(k) += K * (s_i - 0.5) * (s_j - 0.5) * initial_vec(k);
                //Kinetic spin part: S+ S-
                temp = base_vector;
                if (s_i == 0 && s_j == 1) {
                    temp[j] = base_vector[j] - 4; //spin filp
                    temp[next_j] = base_vector[next_j] + 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    result_vec(idx) += K / 2. * initial_vec(k); //   S_i^+ * S_i+1^-
                    result_vec(k) += K / 2. * initial_vec(idx); //   S_i^+ * S_i+1^-
                }
                //---------------------

                // electron hopping j+1 -> j
                    //spin up
                temp = base_vector;
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                    temp[next_j] -= 1; // anihilate spin-up electron
                    temp[j] += 1; // create spin-up electron
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        result_vec(idx) += -t * initial_vec(k);
                        result_vec(k) += -t * initial_vec(idx);
                    }
                    else {
                        result_vec(idx) += t * initial_vec(k);
                        result_vec(k) += t * initial_vec(idx);
                    }
                }
                // spin down
                temp = base_vector;
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        temp[next_j] -= 2; // anihilate spin-down electron
                        temp[j] += 2; // create spin-down electron
                        idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            result_vec(idx) += -t * initial_vec(k);
                            result_vec(k) += -t * initial_vec(idx);
                        }
                        else {
                            result_vec(idx) += t * initial_vec(k);
                            result_vec(k) += t * initial_vec(idx);
                        }
                    }
                }
            }
            // electron repulsion
            if (base_vector[j] == 7 || base_vector[j] == 3)
                result_vec(k) += U * initial_vec(k);
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            temp = base_vector;
            if (base_vector[j] == 5) {// S_i^+ s_i^-
                temp[j] = 2;
                idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                result_vec(idx) += -J_H * initial_vec(k);
                result_vec(k) += -J_H * initial_vec(idx);
            }
            //Diagonal - z part
            if (base_vector[j] == 1 || base_vector[j] == 6)
                result_vec(k) -= 2.0 * J_H * 0.25 * initial_vec(k);
            if (base_vector[j] == 2 || base_vector[j] == 5)
                result_vec(k) += 2.0 * J_H * 0.25 * initial_vec(k);
        }
    }
    return result_vec;
}

void Lanczos::Build_Lanczos_Hamil_wKrylovSpace(vec initial_vec) {
    this->H_L = mat(lanczos_steps, lanczos_steps, fill::zeros);
    //this->Lanczos_GS = vec(N);
    try { this->Krylov_space = mat(N, lanczos_steps, fill::zeros); }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        out << "dim(Krylov_space) = " << Krylov_space.size() * sizeof(Krylov_space(0, 0)) << "\n";
        assert(false);
    }
    //out << "dim(Krylov_space) = " << Krylov_space.size() * sizeof(H(0, 0)) << "\n";

    Krylov_space.col(0) = initial_vec;
    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); // normalized Krylov_space(j=0)

    vec tmp = Hamil_vector_multiply(Krylov_space.col(0)); // tmp = H * Krylov_space(0)
    double alfa = dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);

    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        Krylov_space.col(j) = tmp / beta;

        tmp = Hamil_vector_multiply(Krylov_space.col(j)); // tmp = H * tmp2
        alfa = dot(Krylov_space.col(j), tmp);
        tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;
    }
    tmp.~vec();
}
void Lanczos::Build_Lanczos_Hamil(vec initial_vec) {
    this->H_L = mat(lanczos_steps, lanczos_steps, fill::zeros);
    Hamiltonian();

    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    // already normalized input random vector

    this->randVec_inKrylovSpace = vec(lanczos_steps);
    randVec_inKrylovSpace(0) = dot(initial_vec, initial_vec); // =1

    //vec tmp = Hamil_vector_multiply(initial_vec); // tmp = H * Krylov_space(0)
    vec tmp = H * initial_vec;
    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    vec tmp2_prev = initial_vec;
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        double beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;
        randVec_inKrylovSpace(j) = dot(tmp2, initial_vec);
        tmp = H * tmp2;
        //tmp = Hamil_vector_multiply(tmp2); // tmp = H * tmp2
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * tmp2_prev;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;

        tmp2_prev = tmp2;
        //out << j << "lanczos" << endl;
    }
    tmp.~vec();
    tmp2_prev.~vec();
}
void Lanczos::Lanczos_Diagonalization() {
    srand(time(NULL));
    Build_Lanczos_Hamil(Create_Random_vec());
    eig_sym(eigenvalues, eigenvectors, H_L);
}

void Lanczos::Lanczos_GroundState() {

    srand(time(NULL));
    vec initial_vec = Create_Random_vec();
    Build_Lanczos_Hamil(initial_vec);

    eig_sym(eigenvalues, eigenvectors, H_L);
    vec GS = eigenvectors.col(0);

    Hamiltonian();
    this->ground_state = vec(N, fill::zeros);

    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    ground_state = GS(0) * initial_vec;

    //vec tmp = Hamil_vector_multiply(initial_vec); // tmp = H * Krylov_space(0)
    vec tmp = H * initial_vec;
    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    for (int j = 1; j < lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        ground_state += GS(j) * tmp2;
        tmp = H * tmp2;
        //tmp = Hamil_vector_multiply(tmp2); // tmp = H * tmp2
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        initial_vec = tmp2;
    }
    tmp.~vec();
}

vec Lanczos::RandVec_times_KrylovTranspose(vec randvec) {
    double beta = dot(randvec, randvec);
    randvec = randvec / sqrt(beta); // normalized Krylov_space(j=0)

    vec result_randvec(lanczos_steps, fill::zeros);
    result_randvec(0) = dot(randvec, randvec);

    vec tmp = Hamil_vector_multiply(randvec); // tmp = H * Krylov_space(0)
    double alfa = dot(randvec, tmp);
    tmp = tmp - alfa * randvec;

    vec tmp2_prev = randvec;
    for (int j = 1; j < lanczos_steps; j++) {
        double beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        result_randvec(j) = dot(randvec, tmp2); // overlap <Krylov.col(j)|ranvec>

        tmp = Hamil_vector_multiply(tmp2); // tmp = H * tmp2
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * tmp2_prev;

        tmp2_prev = tmp2;
    }
    tmp.~vec();
    tmp2_prev.~vec();

    return result_randvec;
}
double Lanczos::Cv_kernel(double T, int random_steps) {

    double Z = 0; //partition function
    double E_av = 0; // average energy
    double E2_av = 0; // average squared energy
    double overlap = 0; //overlap of random vectopr and lanczos eigenvectors

    for (int r = 0; r < random_steps; r++) {
        vec rand_vec = Create_Random_vec();
        Build_Lanczos_Hamil(rand_vec);
        eig_sym(eigenvalues, eigenvectors, H_L);

        for (int m = 0; m < lanczos_steps; m++) {
            overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
            overlap *= overlap;
            Z += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            E_av += eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            E2_av += eigenvalues(m) * eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
        }
    }
    E_av = E_av / Z * (double)N / (double)random_steps;
    E2_av = E2_av / Z * (double)N / (double)random_steps;

    return (E2_av - E_av * E_av) / T / T / (L + 0.0);
}
vec Lanczos::Heat_Capacity_Lanczos(int random_steps) {
    double T = dT;
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->chi_0 = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    int k = 0;
    srand(time(NULL));
    while (T <= T_end) {
        Cv(k) = Cv_kernel(T, random_steps);
        out << T << "\t\t" << Cv(k) << endl;
        T += dT; k++;
    }
    return Cv;
}



//----------------------------------------------------------------------------------------------
//----------------------------------------Main subroutines---------------------------------------
//----------------------------------------------------------------------------------------------

void Main_U(int L, int N_e, double t) {
//#pragma omp parallel for 
    for (int k = 2; k < 40; k+=2) {
        double U = (double)k / 10.0;
        double K, J_H;
        K = 4 * 0.15 * 0.15 / U;
        J_H = 0.25 * U;
        //Main_DOS(L, N_e, t, K, U, J_H);
        //Main_Cv(L, N_e, t, K, U, J_H);
        HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
        Hamil.Hamiltonian();
        Hamil.Diagonalization();

        double Ef = FermiLevel(L, N_e, t, K, U, J_H);
        Hamil.printEnergy(Ef);

        //out << "U = " << U << " done!" << endl;
    }
}
void Main_Jh(int L, int N_e, double t, double K, double U) {
    double J_H = 0.05 * U;
    while(J_H <= 0.3*U) {
        Main_DOS(L, N_e, t, K, U, J_H);
        Main_Cv(L, N_e, t, K, U, J_H);
        out << "J_H/U = " << J_H / U << " done!" << endl;
        J_H += 0.05 * U;
    }
}
void Main_Cv(int L, int N_e, double t, double K, double U, double J_H) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    for (int Sz = -N_e; Sz <= N_e; Sz +=2) {
        HamiltonianKH Hamil(L, N_e, t, U, K, J_H, Sz);
        Hamil.Hamiltonian();
        Hamil.Diagonalization();
        Cv += Hamil.Heat_Capacity() / (N_e + 1);
        chi_0 += Hamil.chi_0 / L;
    }
    print_Cv(Cv, U, N_e, L);
    print_chi(chi_0, U, N_e, L);
    Cv.~vec(); chi_0.~vec();
}
void Main_Cv_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
        Lanczos Hamil(L, N_e, t, U, K, J_H, Sz, M);
        Hamil.Lanczos_Diagonalization();
        if (random_steps > Hamil.N) random_steps = Hamil.N;
        Cv += Hamil.Heat_Capacity_Lanczos(random_steps) / double(N_e - 1.0);
        chi_0 += Hamil.chi_0 / std::pow(8, L);
        out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
    }
    print_Cv_Lanczos(Cv, U, N_e, L, M, random_steps);
    print_chi(chi_0, U, N_e, L);
    Cv.~vec(); chi_0.~vec();
}
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H) {
    HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
    Hamil.Hamiltonian();
    Hamil.Diagonalization();
    vector<double> omega_vec = prepareOmegaVec(Hamil.eigenvalues(0), Hamil.eigenvalues(Hamil.eigenvalues.size() - 1), domega);
    vec DOS = Hamil.Total_Density_of_states(omega_vec) / (N_e + 1);
    double maximum = max(DOS);
    double Ef = FermiLevel(L, N_e, t, K, U, J_H);
    Hamil.printEnergy(Ef);
    printDOS(DOS, U, N_e, L, omega_vec, maximum, Ef);
}


// Helpful tools
double FermiLevel(int L, int N_e, double t, double K, double U, double J_H) {
    double Ef;
    Lanczos Object(L, N_e + 1, t, U, K, J_H, ((N_e + 1) % 2 == 0) ? 0 : 1, 300);
    Object.Lanczos_Diagonalization();
    Ef = Object.eigenvalues(0);

    Lanczos Object2(L, N_e - 1, t, U, K, J_H, ((N_e - 1) % 2 == 0) ? 0 : 1, 300);
    Object2.Lanczos_Diagonalization();
    Ef = (Ef - Object2.eigenvalues(0)) / 2.0;

    return Ef;
}
std::vector<double> prepareOmegaVec(double omega_min, double omega_max, double dOmega) {
    vector<double> omega_vec;
    double omega = omega_min;
    while (omega <= omega_max) {
        omega_vec.push_back(omega);
        omega += dOmega;
    }
    return omega_vec;
}

void printDOS(vec resultDOS, double U, double N_e, int L, std::vector<double> omega_vec, double maximum, double E_fermi) {
    ofstream DOSfile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    DOSfile.open("DOS_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < omega_vec.size(); k++)
        DOSfile << omega_vec[k] - E_fermi << "\t\t" << resultDOS[k] / maximum << endl;

    DOSfile.close();
}
void print_Cv(vec Cv, double U, double N_e, int L) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("C_V_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << Cv(k) << endl; //save heat capacity to file
        T += dT; k++;
    }
    savefile.close();
}
void print_Cv_Lanczos(vec Cv, double U, double N_e, int L, int M, int random_steps) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    //savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + "_M=" + std::to_string(M) + "_R=" + std::to_string(random_steps) + ".txt");
    savefile.open("Cv_Lanczos_M=" + std::to_string(M) + "_R=" + std::to_string(random_steps) + ".txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << Cv(k) << endl; //save heat capacity to file
        T += dT; k++;
    }
    savefile.close();
}
void print_chi(vec chi, double U, double N_e, int L) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("chi_0_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << chi(k) << endl; //save heat capacity to file
        T += dT; k++;
    }
    savefile.close();
}
void HamiltonianKH::printEnergy(double Ef) {
    ofstream Efile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)num_of_electrons / (double)L;
    Efile.open("E_L=" + std::to_string(L) + ",N_e="+std::to_string(num_of_electrons)+"_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < 20; k++)
        Efile << U << "\t\t" << eigenvalues(k) << endl;

    Efile.close();

}

// Conversion of int to binary vector - using modulo operator
vector<int> int_to_binary(int idx, int L) {
    vector<int> vec(L);
    int temp = idx;
    for (int k = 0; k < L; k++) {
        vec[vec.size() - 1 - k] = temp % 8;
        temp = static_cast<int>((double)temp / 8.);
    }
    return vec;
}
int binary_to_int(vector<int> vec) {
    int val = 0;
    for (int k = 0; k < vec.size(); k++) {
        val += vec[vec.size() - 1 - k] * std::pow(8, k);
    }
    return val;
}

int binary_search(std::vector<unsigned long int> arr, int l_point, int r_point, int element) {
    if (r_point >= l_point) {
        int middle = l_point + (r_point - l_point) / 2;
        if (arr[middle] == element) return middle;
        else if (arr[middle] > element) return binary_search(arr, l_point, middle - 1, element);
        else return binary_search(arr, middle + 1, r_point, element);
    }
    assert(false && "Element not present in the array");
    return -1;
}
