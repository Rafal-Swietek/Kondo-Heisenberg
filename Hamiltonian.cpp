#include "Hamiltonian.h"

#define im std::complex<double>(0.0,1.0)
#define M_PI 3.14159265358979323846

double pi = M_PI;

double dT = 0.05;
double T_end = 10.0;
double domega = 0.01;

std::mutex my_mut;

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

	this->mapping = vector<ull_int>();
	generate_mapping();
	this->N = mapping.size();
}
HamiltonianKH::HamiltonianKH() {}
//-------------------------

void HamiltonianKH::setHamiltonianElem(ull_int k, double value, std::vector<int> temp){
    ull_int idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
        H(idx, k) += value;
        H(k, idx) += value;
}
void HamiltonianKH::Hamiltonian() {
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian

        int s_i, s_j; //i=j, j=j+1
	    bool PBC = 0; //allows periodic boundary conditions if =1
	    int next_j;
        for (ull_int k = 0; k < N; k++){
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
    }   catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
   //out << "dim(H) = " << (sizeof(H) + H.n_elem * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
}

void HamiltonianKH::setHamiltonianElem_sparse(ull_int k, double value, std::vector<int> temp) {
    ull_int idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
    H_sparse(idx, k) += value;
    H_sparse(k, idx) += value;
}
void HamiltonianKH::Hamiltonian_sparse() {
    try {
        this->H_sparse = sp_mat(N, N); //hamiltonian

        int s_i, s_j; //i=j, j=j+1
        bool PBC = 0; //allows periodic boundary conditions if =1
        int next_j;
        for (ull_int k = 0; k < N; k++) {
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
                    H_sparse(k, k) += K * (s_i - 0.5) * (s_j - 0.5);

                    //Kinetic spin part: S+ S-
                    temp = base_vector;
                    if (s_i == 0 && s_j == 1) { // S_i^+ S_i+1^-
                        temp[j] = base_vector[j] - 4; //spin filp
                        temp[next_j] = base_vector[next_j] + 4;
                        setHamiltonianElem_sparse(k, K / 2., temp);
                    }
                    // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                        //spin up
                    temp = base_vector;
                    //only odd numbers have up-electrons  //even numbers lack one up-electron
                    if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                        temp[next_j] -= 1; // anihilate spin-up electron
                        temp[j] += 1; // create spin-up electron
                        if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                            setHamiltonianElem_sparse(k, -t, temp);
                        }
                        else  setHamiltonianElem_sparse(k, +t, temp);
                    }
                    //spin down
                    temp = base_vector;
                    // the if statement contains every possible down-electron hopping: next_j->j
                    if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                        if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                            temp[next_j] -= 2; // anihilate spin-down electron
                            temp[j] += 2; // create spin-down electron
                            if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                                setHamiltonianElem_sparse(k, -t, temp);
                            }
                            else  setHamiltonianElem_sparse(k, +t, temp);
                        }
                    }
                }
                // electron repulsion
                if (base_vector[j] == 7 || base_vector[j] == 3)
                    H_sparse(k, k) += U;
                // electron-localised spin interaction ( interorbital electronic spin interaction)
                temp = base_vector;
                if (base_vector[j] == 5) {// S_i^+ s_i^-
                    temp[j] = 2;
                    setHamiltonianElem_sparse(k, -J_H, temp);
                }
                //Diagonal - z part
                if (base_vector[j] == 1 || base_vector[j] == 6)
                    H_sparse(k, k) -= 2.0 * J_H * 0.25;
                if (base_vector[j] == 2 || base_vector[j] == 5)
                    H_sparse(k, k) += 2.0 * J_H * 0.25;
            }
            //if (k % (ull_int)std::pow(10, 7) == 0) out << k << endl;
        }
    } catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    //out << "dim(H) = " << (sizeof(H_sparse) + H_sparse.n_elem * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
}

//generates the vector, which maps the base_vector index to the index in given subblock
std::tuple<int, int, int> calculateSpinElements(int L, ull_int& j) {
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
void HamiltonianKH::mapping_kernel(int& L, ull_int&start, ull_int& stop, std::vector<ull_int>& map_threaded, int& _id, int& S_z, int& no_electrons) {
    for (ull_int j = start; j < stop; j++) {
        int bSz = 0, fSz = 0, N_e = 0;
        std::tie(bSz, fSz, N_e) = calculateSpinElements(L, j);
        if ((bSz + fSz == S_z) && N_e == no_electrons)
            map_threaded.push_back(j);
    }
}
void HamiltonianKH::generate_mapping() {
    ull_int start = 0, stop = std::pow(8, L);
    int id = 0;
    int S_z = this->Sz;
    int no_electrons = this->num_of_electrons;
    mapping_kernel(L, start, stop, mapping, id, S_z, num_of_electrons);
    //Threaded
    /*std::vector<vector<ull_int>> map_threaded(num_of_threads);
    std::vector<thread> threads;
    threads.reserve(num_of_threads);
    ull_int size = 0;
    for (int t = 0; t < num_of_threads; t++) {
        start = t * (ull_int)std::pow(8, L) / num_of_threads;
        stop = ((t + 1) == num_of_threads ? (ull_int)std::pow(8, L) : (ull_int)std::pow(8, L) * (t + 1) / num_of_threads);
        map_threaded[t] = std::vector<ull_int>();
        threads.emplace_back(&HamiltonianKH::mapping_kernel, this, ref(L), ref(start), ref(stop), ref(map_threaded[t]), ref(t), ref(S_z), ref(no_electrons));
        out << "Thread " << t << " joined tha party! from " << start << " to " << stop << endl;
    }
    for (auto& t : threads) t.join(); 
    for (auto& t : threads) t.~thread();
    ull_int size;
    for (int t = num_of_threads - 1; t >= 0; t--) size += map_threaded[t].size();
    mapping = std::vector<ull_int>(size);
    for (int t = 0; t < num_of_threads; t++){
        mapping.insert(mapping.begin(), map_threaded[t].begin(), map_threaded[t].end());
    }*/
    out << "Mapping generated with  " << mapping.size() << "  elements" << endl;
    //out << mapping[0] << " " << mapping[mapping.size() - 1] << endl;
    assert(mapping.size() > 0 && "Not possible number of electrons - no. of states < 1");
}
//----------------------------------------------------

void HamiltonianKH::Diagonalization() {
	try {
		arma::eig_sym(eigenvalues, eigenvectors, H);
        this->ground_state = eigenvectors.col(0);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
        out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
        	assert(false);
    }
    //out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
}

void HamiltonianKH::print_base_vector(std::vector<int>& base_vector) {
    out << " |";
    for (int l = 0; l < L; l++)
        out << base_vector[l];
    out << ">" << endl;
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
        for (ull_int n = 0; n < N; n++)
            DOS += -1. / (double)L / pi * cpx(1. / (omega + 2 * domega * 1i - eigenvalues(n))).imag();

        if (DOS > maximum)
             maximum = DOS;

        resultDOS(w) = DOS;
    }
    //printDOS(resultDOS,U, num_of_electrons,L,omega_vec,maximum, E_F);
    return resultDOS;
}

void HamiltonianKH::show_ground_state() {
    vec GS((ull_int)std::pow(2, L), fill::zeros);
    for (ull_int k = 0; k < N; k++) {
        std::vector<int> base_vector = int_to_binary(mapping[k], L);
        int val = 0;
        for (int j = 0; j < L; j++)
            val += (1 - int(base_vector[base_vector.size() - 1 - j] / 4)) * std::pow(2, j);
        GS(val) += ground_state(k);
    }
    GS = arma::abs(GS);
    GS = GS / dot(GS, GS); //normalizing to be sure
    out << endl;
    double maximum = arma::max(GS);
    for (ull_int k = 0; k < GS.size(); k++) {
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

    this->mapping = vector<ull_int>();
    generate_mapping();
    this->N = mapping.size(); //out << "dim = " << N * sizeof(mapping[0]) << endl;
}
Lanczos::~Lanczos() {}

vec Lanczos::Create_Random_vec() {
    vec random_vec(N, fill::zeros);
    double norm = 0;
    for (ull_int j = 0; j < N; j++) {
        random_vec(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
        norm += random_vec(j) * random_vec(j);
    }
    return random_vec / norm;
}
vec Lanczos::Hamil_vector_multiply(vec initial_vec) {
    vec result_vec(N, fill::zeros);
    for (ull_int k = 0; k < N; k++) {
        std::vector<int> base_vector = int_to_binary(mapping[k], L);
        std::vector<int> temp(base_vector);
        int PBC = 0;
        int next_j;
        ull_int idx;
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

void Lanczos::Build_Lanczos_Hamil(vec initial_vec) {
    this->H_L = mat(lanczos_steps, lanczos_steps, fill::zeros);

    double beta = arma::dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    // already normalized input random vector

    this->randVec_inKrylovSpace = vec(lanczos_steps);
    randVec_inKrylovSpace(0) = arma::dot(initial_vec, initial_vec); // =1

    vec tmp = H_sparse * initial_vec;
    double alfa = arma::dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;
    //out << tmp.t();
    vec tmp2_prev = initial_vec;
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        double beta = sqrt(arma::dot(tmp, tmp));
        vec tmp2 = tmp / beta;
        randVec_inKrylovSpace(j) = dot(tmp2, initial_vec);
        tmp = H_sparse * tmp2;
        alfa = arma::dot(tmp2, tmp);
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

    Hamiltonian_sparse();
    srand(time(NULL));
    vec initial_vec = Create_Random_vec();
    Build_Lanczos_Hamil(initial_vec);

    eig_sym(eigenvalues, eigenvectors, H_L);
    vec GS = eigenvectors.col(0);

    this->ground_state = vec(N, fill::zeros);

    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    ground_state = GS(0) * initial_vec;

    //vec tmp = Hamil_vector_multiply(initial_vec); // tmp = H * Krylov_space(0)
    vec tmp = H_sparse * initial_vec;
    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    for (int j = 1; j < lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        ground_state += GS(j) * tmp2;
        tmp = H_sparse * tmp2;
        //tmp = Hamil_vector_multiply(tmp2); // tmp = H * tmp2
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        initial_vec = tmp2;
    }
    tmp.~vec();
}

vec Lanczos::Heat_Capacity_Lanczos(int random_steps) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec E_av(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec E_av2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->chi_0 = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    srand(time(NULL));
    Hamiltonian_sparse();
    for (int r = 0; r < random_steps; r++) {
        vec rand_vec = Create_Random_vec();
        Build_Lanczos_Hamil(rand_vec);
        eig_sym(eigenvalues, eigenvectors, H_L);
        double T = dT; int k = 0;
        while (T <= T_end) {
            double overlap;
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
                overlap *= overlap;
                Z(k) += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
                E_av(k) += eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
                E_av2(k) += eigenvalues(m) * eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            }
            T += dT; k++;
        }
    }
    double T = dT; int k = 0;
    while (T <= T_end) {
        E_av(k) = E_av(k) / Z(k) * (double)N / (double)random_steps;
        E_av2(k) = E_av2(k) / Z(k) * (double)N / (double)random_steps;
        Cv(k) = (E_av2(k) - E_av(k) * E_av(k)) / T / T / (double)L;
        T += dT; k++;
    }
    return Cv;
}

vec Lanczos::thermal_average_lanczos(vec& quantity, int& random_steps){
    assert(quantity.n_elem > 0 && "given object is empty");
    double T = dT;
    int t = 0;
    vec result(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    while(T <= T_end) {
        double Z = 0; //partition function
        double overlap = 0; //overlap of random vectopr and lanczos eigenvectors
        for (int m = 0; m < lanczos_steps; m++) {
            overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
            overlap *= overlap;
            Z += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            result(t) += quantity(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
        }
        result(t) = result(t) / Z * (double)N / (double)random_steps;
        T += dT; t++;
    }
    return result;
}
//----------------------------------------------------------------------------------------------
//----------------------------------------Main subroutines---------------------------------------
//----------------------------------------------------------------------------------------------
void Heat_Capacity(std::vector<arma::vec>& energies, vec& Cv) {
    std::vector temperature = prepare_parameterVec(dT, T_end, dT);
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k];
        double heat_capacity = 0;
        for (int l = 0; l < energies.size(); l++) {
            double Sz = -((double)energies.size() - 1.0) / 2. + (double)l;
            double E_av = 0;
            double E2_av = 0;
            double Z = 0;
            //out << Sz << endl;
            for (int n = 0; n < energies[l].size(); n++) {
                E_av += energies[l](n) * std::exp(-(energies[l](n) - energies[l](0)) / T);
                E2_av += energies[l](n) * energies[l](n) * std::exp(-(energies[l](n) - energies[l](0)) / T);
                Z += std::exp(-(energies[l](n) - energies[l](0)) / T);
            }
            E_av = E_av / Z;
            E2_av = E2_av / Z;
            heat_capacity += (E2_av - E_av * E_av) / T / T;
        }
        int L = (energies.size() - 1) * 2 / 3; // n=1.5 filling
        Cv(k) = heat_capacity / (double)L / (double)energies.size();

        //out << T << " " << Cv(k) << " " << L << endl;
    }
}
void static_spin_susceptibility(std::vector<arma::vec>& energies, vec& chi) {
    std::vector temperature = prepare_parameterVec(dT, T_end, dT);
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k]; 
        double Z = 0;
        double X_0 = 0;
        for (int l = 0; l < energies.size(); l++) {
            double Sz = -((double)energies.size() - 1.0) / 2. + (double)l;
            //out << Sz << endl;
            for (int n = 0; n < energies[l].size(); n++) {
                Z += std::exp(-energies[l](n) / T);
                X_0 += Sz * Sz * std::exp(-energies[l](n) / T);
            }
        }
        int L = (energies.size() - 1) * 2 / 3; // n=1.5 filling
        chi(k) = X_0 / Z / T / (double)L;

        out << T << " " << chi(k) << " " << L << endl;
    }
}

void Main_U(int L, int N_e, double t) {
//#pragma omp parallel for 
    for (int k = 2; k < 40; k+=2) {
        double U = (double)k / 10.0;
        double K, J_H;
        K = 4 * 0.15 * 0.15 / U;
        J_H = 0.25 * U;
        //Main_DOS(L, N_e, t, K, U, J_H);
        Main_Cv(L, N_e, t, K, U, J_H);
        /*HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
        Hamil.Hamiltonian();
        Hamil.Diagonalization();

        double Ef = FermiLevel(L, N_e, t, K, U, J_H);
        Hamil.printEnergy(Ef);*/

        out << "U = " << U << " done!" << endl;
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
    std::vector<arma::vec> energies;
    int idx = 0;
    for (int Sz = -N_e; Sz <= N_e; Sz +=2) {
        HamiltonianKH Hamil(L, N_e, t, U, K, J_H, Sz);
        Hamil.Hamiltonian();
        Hamil.Diagonalization();
        energies.push_back(Hamil.eigenvalues);
        out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
    }
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    Heat_Capacity(energies, chi_0);
    static_spin_susceptibility(energies, chi_0);
    print_Cv(Cv, U, N_e, L);
    print_chi(chi_0, U, N_e, L);
    Cv.~vec(); chi_0.~vec();
}
void Main_Cv_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
        Lanczos Hamil(L, N_e, t, U, K, J_H, Sz, M);
        if (M > Hamil.N) {
            Hamil.lanczos_steps = Hamil.N;
            out << "Lanczos steps greater than system size" << endl;
        }
        Cv += Hamil.Heat_Capacity_Lanczos(random_steps) / double(N_e + 1.0);
        //chi_0 += Hamil.chi_0 / double(N_e + 1.0);
        //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
    }
    print_Cv_Lanczos(Cv, U, N_e, L, M, random_steps);
    //print_chi(chi_0, U, N_e, L);
    Cv.~vec();
}
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H) {
    HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
    Hamil.Hamiltonian();
    Hamil.Diagonalization();
    vector<double> omega_vec = prepare_parameterVec(Hamil.eigenvalues(0), Hamil.eigenvalues(Hamil.eigenvalues.size() - 1), domega);
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
std::vector<double> prepare_parameterVec(double _min, double _max, double step) {
    vector<double> parameter_vec;
    double temp = _min;
    while (temp <= _max) {
        parameter_vec.push_back(temp);
        temp += step;
    }
    return parameter_vec;
}

void printDOS(vec& resultDOS, double U, double N_e, int L, std::vector<double>& omega_vec, double maximum, double E_fermi) {
    ofstream DOSfile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    DOSfile.open("DOS_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < omega_vec.size(); k++)
        DOSfile << omega_vec[k] - E_fermi << "\t\t" << resultDOS[k] / maximum << endl;

    DOSfile.close();
}
void print_Cv(vec& Cv, double U, double N_e, int L) {
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
void print_Cv_Lanczos(vec& Cv, double U, double N_e, int L, int M, int random_steps) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    //savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + "_M=" + std::to_string(M) + "_R=" + std::to_string(random_steps) + ".txt");
    savefile.open("Cv_Lanczos_M=" + std::to_string(M) + "_R=" + std::to_string(random_steps) + "_gKH.txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << Cv(k) << endl; //save heat capacity to file
        T += dT; k++;
    }
    savefile.close();
}
void print_chi(vec& chi, double U, double N_e, int L) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("chi_0_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << chi(k) << endl;
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
vector<int> int_to_binary(ull_int idx, int L) {
    vector<int> vec(L);
    ull_int temp = idx;
    for (int k = 0; k < L; k++) {
        vec[vec.size() - 1 - k] = static_cast<int>(temp % 8);
        temp = static_cast<ull_int>((double)temp / 8.);
    }
    return vec;
}
ull_int binary_to_int(vector<int> vec) {
    int val = 0;
    for (int k = 0; k < vec.size(); k++) {
        val += vec[vec.size() - 1 - k] * std::pow(8, k);
    }
    return val;
}

ull_int binary_search(std::vector<ull_int> arr, int l_point, int r_point, ull_int element) {
    if (r_point >= l_point) {
        int middle = l_point + (r_point - l_point) / 2;
        if (arr[middle] == element) return middle;
        else if (arr[middle] > element) return binary_search(arr, l_point, middle - 1, element);
        else return binary_search(arr, middle + 1, r_point, element);
    }
    assert(false && "Element not present in the array");
    return -1;
}
