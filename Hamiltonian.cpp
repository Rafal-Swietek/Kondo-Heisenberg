#include "Hamiltonian.h"

double pi = M_PI;

double T = 0.01;

double dT = 0.001;
double T_end = 3.0;
double domega = 0.005;
double eta = 0.02;

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

	this->mapping = new vector<ull_int>();
	generate_mapping();
	this->N = mapping->size();
    Hamiltonian();
}
HamiltonianKH::HamiltonianKH() {}
//-------------------------

void HamiltonianKH::update_parameters(double t, double U, double K, double J_H) {
    this->t = t; this->U = U, this->K = K; this->J_H = J_H;
}

void HamiltonianKH::setHamiltonianElem(ull_int& k, double value, std::vector<int>&& temp){
    ull_int idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
        H(idx, k) += value;
        H(k, idx) += value;
}
void HamiltonianKH::Hamiltonian() {
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    int s_i, s_j; //i=j, j=j+1
	bool PBC = 0; //allows periodic boundary conditions if =1
	int next_j;
    std::vector<int> base_vector(L);
    vector<int> temp = base_vector;
    for (ull_int k = 0; k < N; k++){
		int_to_binary(mapping->at(k), base_vector);
		temp = base_vector;
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
                    setHamiltonianElem(k, K / 2., std::move(temp));
                }
                // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                    //spin up
                temp = base_vector;
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                    temp[next_j] -= 1; // anihilate spin-up electron
                    temp[j] += 1; // create spin-up electron
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        setHamiltonianElem(k, -t, std::move(temp));
                    }
                    else  setHamiltonianElem(k, +t, std::move(temp));
                }
                //spin down
                temp = base_vector;
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        temp[next_j] -= 2; // anihilate spin-down electron
                        temp[j] += 2; // create spin-down electron
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            setHamiltonianElem(k, -t, std::move(temp));
                        }
                        else  setHamiltonianElem(k, +t, std::move(temp));
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
                    setHamiltonianElem(k,-J_H, std::move(temp));
                }
		    //Diagonal - z part
				if (base_vector[j] == 1 || base_vector[j] == 6)
					H(k, k) -= 2.0 * J_H * 0.25;
				if (base_vector[j] == 2 || base_vector[j] == 5)
					H(k, k) += 2.0 * J_H * 0.25;
        }
	}
   //out << "dim(H) = " << (sizeof(H) + H.n_elem * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
}

ull_int findElement(std::vector<ull_int>& vector, ull_int element) {
    std::vector<ull_int>::iterator it = find(vector.begin(), vector.end(), element);
    assert(it != vector.end() && "Element not present in the array");
    return (ull_int)std::distance(vector.begin(), it);
}
void HamiltonianKH::setHamiltonianElem_sparse(ull_int& k, double value, std::vector<int>&& temp) {
    ull_int idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));// findElement(std::move(mapping), binary_to_int(std::move(temp))); //
    assert(idx < N && "Somehow index out of scope, wtf?? Found element not possibly present in the array");
    H_sparse(idx, k) += value;
    H_sparse(k, idx) += value;
}
void HamiltonianKH::Hamiltonian_sparse() {
    try {
        this->H_sparse = sp_mat(N, N); //hamiltonian
    }catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    std::vector<int> base_vector(L);
    vector<int> temp(base_vector);
//#pragma omp parallel for num_threads(num_of_threads)
    for (ull_int k = 0; k < N; k++) {
        if (show_system_size_parameters == true && k % (mapping->at(N - 1) / 10) == 0)
            out << k << "done" << endl;
        int_to_binary(mapping->at(k), base_vector);
        temp = base_vector;
        int s_i, s_j; //i=j, j=j+1
        bool PBC = 0; //allows periodic boundary conditions if =1
        int next_j;
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
                    setHamiltonianElem_sparse(k, K / 2., std::move(temp));
                }
                // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                    //spin up
                temp = base_vector;
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                    temp[next_j] -= 1; // anihilate spin-up electron
                    temp[j] += 1; // create spin-up electron
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        setHamiltonianElem_sparse(k, -t, std::move(temp));
                    }
                    else  setHamiltonianElem_sparse(k, +t, std::move(temp));
                }
                //spin down
                temp = base_vector;
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        temp[next_j] -= 2; // anihilate spin-down electron
                        temp[j] += 2; // create spin-down electron
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            setHamiltonianElem_sparse(k, -t, std::move(temp));
                        }
                        else  setHamiltonianElem_sparse(k, +t, std::move(temp));
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
                setHamiltonianElem_sparse(k, -J_H, std::move(temp));
            }
            //Diagonal - z part
            if (base_vector[j] == 1 || base_vector[j] == 6)
                H_sparse(k, k) -= 2.0 * J_H * 0.25;
            if (base_vector[j] == 2 || base_vector[j] == 5)
                H_sparse(k, k) += 2.0 * J_H * 0.25;
        }
        //if (k % (ull_int)std::pow(10, 2) == 0) out << k << endl;
    }
    if (show_system_size_parameters)
        out << "Hamiltonian complete" << endl;
    //out << "dim(H) = " << (sizeof(H_sparse) + H_sparse.n_elem * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
}

//generates the vector, which maps the base_vector index to the index in given subblock
std::tuple<int, int, int> calculateSpinElements(int L, ull_int& j, std::vector<int>& temp) {
    int bSz = 0; //bosonic total spin - spin of upper orbital locked to n=1 filling
    int fSz = 0; //fermionic total spin
    int N_e = 0; // numer of electrons in given state
    int_to_binary(j, temp);

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
void HamiltonianKH::mapping_kernel(ull_int start, ull_int stop, std::vector<ull_int>* map_threaded, int _id) {
    int n = 1;
    //out << "A new thread joined tha party! from " << start << " to " << stop << endl;
    std::vector<int> temp(L);
    for (ull_int j = start; j < stop; j++) {
        int bSz = 0, fSz = 0, N_e = 0;
        std::tie(bSz, fSz, N_e) = calculateSpinElements(this->L, j, temp);
        if ((bSz + fSz == this->Sz) && N_e == this->num_of_electrons) 
             map_threaded->push_back(j);
        if (show_system_size_parameters && (j - start) % ull_int((stop - start) * n / 4) == 0 && j > 0) { 
            out << n << "-th quarter of " << _id << endl; 
            n++; 
        }
    }
}
void HamiltonianKH::generate_mapping() {
    ull_int start = 0, stop = std::pow(8, L);
    //mapping_kernel(start, stop, mapping, L, Sz, num_of_electrons);
    //Threaded
    std::vector<std::vector<ull_int>*> map_threaded(num_of_threads);
    std::vector<std::thread> threads;
    threads.reserve(num_of_threads);
    for (int t = 0; t < num_of_threads; t++) {
        start = t * (ull_int)std::pow(8, L) / num_of_threads;
        stop = ((t + 1) == num_of_threads ? (ull_int)std::pow(8, L) : (ull_int)std::pow(8, L) * (t + 1) / num_of_threads);
        map_threaded[t] = new std::vector<ull_int>();
        threads.emplace_back(&HamiltonianKH::mapping_kernel, this, start, stop, map_threaded[t], t);
    }
    for (auto& t : threads) t.join();

    for (auto & t : map_threaded)
        mapping->insert(mapping->end(), t->begin(), t->end());
    //sort(mapping->begin(), mapping->end());
    if (show_system_size_parameters) {
        out << "Mapping generated with  " << mapping->size() << "  elements" << endl;
        out << "Last element = " << mapping->at(mapping->size() - 1) << endl;
    }
    //out << mapping[0] << " " << mapping[mapping.size() - 1] << endl;
    assert(mapping->size() > 0 && "Not possible number of electrons - no. of states < 1");
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

void HamiltonianKH::print_base_vector(std::vector<int>& base_vector, std::ofstream& out_str) {
    out_str << " |";
    for (auto it = base_vector.begin(); it != base_vector.end(); ++it)
        out_str << *it;
    out_str << ">" << endl;
}
vec HamiltonianKH::Total_Density_of_states(std::vector<double>&& omega_vec) {
    vec resultDOS(omega_vec.size());

    double maximum = 0;
#pragma omp parallel for shared(omega_vec, resultDOS) num_threads(16)
    for (int w = 0; w < omega_vec.size(); w++) {
        double omega = omega_vec[w];
        double DOS = 0;
    //#pragma omp parallel for shared(omega_vec, resultDOS) reduction(+: DOS)
        for (ull_int n = 0; n < N; n++)
            DOS += -1. / (double)L / pi * cpx(1. / (omega + eta * 1i - eigenvalues(n))).imag();

        resultDOS(w) = DOS;
    }
    return resultDOS;
}

vec HamiltonianKH::static_structure_factor(double T) { /// here sth wrong?., check
    vec static_structure_factor(L + 1, fill::zeros);
    int number_of_thr = L + 1; // if - elseif - else statement
#pragma omp parallel for num_threads(L+1)
    for (int l = 0; l <= L; l++) {
        double q = (double)l * pi / ((double)L + 1.0);
        cx_vec cpx_vec, a;
        std::vector<int> vect(L);
        sp_cx_mat Sq(sp_mat(N, N), sp_mat(N, N));
        for (int p = 0; p < N; p++) {
            int_to_binary(mapping->at(p), vect);
            Sq(p, p) = 0;
            for (int m = 0; m < L; m++) {
                double Szm = 0;
                if (vect[m] < 4) Szm += 0.5;
                else Szm -= 0.5;
                if (vect[m] % 4 == 1) Szm += 0.5;
                else if (vect[m] % 4 == 2) Szm -= 0.5;
                Sq(p, p) += std::exp(1i * q * (m + 0.0)) * Szm;
            }
        }
        cpx Sq0 = 0;
        if (T > 0) {
            for (int n = 0; n < N; n++) {
                cpx_vec = cx_vec(eigenvectors.col(n), vec(N, fill::zeros));
                a = Sq * cpx_vec;
                Sq0 += std::exp(-eigenvalues(n) / T) * cdot(a, a);
            }
        }
        else {
            cpx_vec = cx_vec(eigenvectors.col(0), vec(N, fill::zeros));
            a = Sq * cpx_vec;
            Sq0 = cdot(a,a);
        }
        /*Sq0 = 0;
        mat cor_mat = correlation_matrix();
        for (int m = 0; m < L; m++) {
            for (int k = 0; k < L; k++) {
                Sq0 += std::exp(1i * q * (m - k + 0.0)) * cor_mat(m,k);
            }
        }*/
        static_structure_factor(l) = real(2.0 * Sq0 / pi / (L + 1.0));
    }
    return static_structure_factor;
}
double HamiltonianKH::partition_function(double T) {
    double Z = 0;
    if (T > 0) {
        for (int n = 0; n < eigenvalues.size(); n++)
            Z += std::exp(-eigenvalues(n) / T);
    }
    else
        Z = 1;
    return Z;
}

void HamiltonianKH::show_ground_state() {
    ofstream GSfile("Ground state for L = " + std::to_string(L) + ".txt");
    vec GS((ull_int)std::pow(2, L), fill::zeros);
    std::vector<int> base_vector(L);
    for (ull_int k = 0; k < N; k++) {
        int_to_binary(mapping->at(k), base_vector);
        int val = 0;
        for (int j = 0; j < L; j++)
            val += (1 - int(base_vector[base_vector.size() - 1 - j] / 4)) * std::pow(2, j);
        GS(val) += ground_state(k) * ground_state(k);
    }
    //GS = arma::abs(GS);
    //GS = GS / dot(GS, GS); //normalizing to be sure
    GSfile << endl;
    double maximum = arma::max(GS);
    std::vector<int> vec(L);
    for (ull_int k = 0; k < GS.size(); k++) {
        if (std::fabs(GS(k)) >= 0.2 * maximum) {
            ull_int temp = k;
            for (int p = 0; p < L; p++) {
                vec[vec.size() - 1 - p] = temp % 2;
                temp = static_cast<int>((double)temp / 2.);
            }
            GSfile << "Ground state:\t";
            print_base_vector(vec, GSfile);
            GSfile << " with probability\t p=" << GS(k) * GS(k) << endl << endl;
        }
    }
    GSfile.close();
}

mat HamiltonianKH::correlation_matrix() {
    mat cor_mat(L, L, fill::zeros);
    vector<int> vect(L);
    for (int p = 0; p < N; p++) {
        int_to_binary(mapping->at(p), vect);
        for (int m = 0; m < L; m++) {
            double Szm = 0;
            if (vect[m] < 4) Szm += 0.5;
            else Szm -= 0.5;
            if (vect[m] % 4 == 1) Szm += 0.5;
            else if (vect[m] % 4 == 2) Szm -= 0.5;

            double Szk = 0;
            for (int k = 0; k < L; k++) {
                double Szk = 0;
                if (vect[k] < 4) Szk += 0.5;
                else Szk -= 0.5;
                if (vect[k] % 4 == 1) Szk += 0.5;
                else if (vect[k] % 4 == 2) Szk -= 0.5;
                cor_mat(m, k) += ground_state(p) * Szm * Szk * ground_state(p);
            }
        }
    }
    return cor_mat;
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

    this->mapping = new vector<ull_int>();
    generate_mapping();
    this->N = mapping->size(); 
    if(show_system_size_parameters) 
        out << "dim = " << N << endl;
}
Lanczos::Lanczos(std::unique_ptr<Lanczos>& obj) {
    this->L = obj->L; //number of sites
    this->num_of_electrons = obj->num_of_electrons; //number of electrons in lower orbital
    this->t = obj->t; 
    this->U = obj->U;
    this->K = obj->K;
    this->J_H = obj->J_H;
    this->Sz = obj->Sz;
    this->lanczos_steps = obj->lanczos_steps;

    this->mapping = obj->mapping;
    this->N = this->mapping->size();
    if (show_system_size_parameters == true)
        out << "dim = " << N << endl;
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
void Lanczos::Hamil_vector_multiply_kernel(int start, int stop, vec& initial_vec, vec& result_vec_threaded) {
    std::vector<int> base_vector(L);
    int PBC = 0;
    int next_j;
    ull_int idx;
    int s_i, s_j;
    for (int k = start; k < stop; k++) {
        int_to_binary(mapping->at(k), base_vector);
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;

                if (base_vector[j] < 4) s_i = 1;
                else s_i = 0;
                if (base_vector[next_j] < 4) s_j = 1;
                else s_j = 0;

                //localised spin diagonal part
                result_vec_threaded(k) += K * (s_i - 0.5) * (s_j - 0.5) * initial_vec(k);
                //Kinetic localised spin part: S+ S-
                if (s_i == 0 && s_j == 1) {
                    base_vector[j] -= 4; //spin filp
                    base_vector[next_j] += 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    result_vec_threaded(k) += K / 2. * initial_vec(idx); //   S_i^+ * S_i+1^-
                    base_vector[j] += 4; //spin filp
                    base_vector[next_j] -= 4;
                }
                if (s_i == 1 && s_j == 0) {
                    base_vector[j] += 4; //spin filp
                     base_vector[next_j] -= 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    result_vec_threaded(k) += K / 2. * initial_vec(idx); //   S_i^+ * S_i+1^-
                    base_vector[j] -= 4; //spin filp
                    base_vector[next_j] += 4;
                }
                //---------------------

                // electron hopping
                //spin up
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) { //  j+1 -> j
                    base_vector[next_j] -= 1; // anihilate spin-up electron
                    base_vector[j] += 1; // create spin-up electron
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    base_vector[next_j] += 1;
                    base_vector[j] -= 1;
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        result_vec_threaded(k) += -t * initial_vec(idx);
                    }
                    else
                        result_vec_threaded(k) += t * initial_vec(idx);
                }
                if (base_vector[j] % 2 == 1 && base_vector[next_j] % 2 == 0) { // j -> j+1
                    base_vector[next_j] += 1; // anihilate spin-up electron
                    base_vector[j] -= 1; // create spin-up electron
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    base_vector[next_j] -= 1;
                    base_vector[j] += 1;
                    if (base_vector[j] % 2 == 1 && base_vector[next_j] % 4 == 2) { // sign change when up-electron hop over down-electron (if 3=du)
                        result_vec_threaded(k) += -t * initial_vec(idx);
                    }
                    else
                        result_vec_threaded(k) += t * initial_vec(idx);
                }
                // spin down
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        base_vector[next_j] -= 2; // anihilate spin-down electron
                        base_vector[j] += 2; // create spin-down electron
                        idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                        base_vector[next_j] += 2;
                        base_vector[j] -= 2;
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            result_vec_threaded(k) += -t * initial_vec(idx);
                        }
                        else
                            result_vec_threaded(k) += t * initial_vec(idx);
                    }
                }
                // the if statement contains every possible down-electron hopping: j->next_j
                if (base_vector[j] % 4 == 2 || base_vector[j] % 4 == 3) {
                    if (base_vector[next_j] % 4 == 0 || base_vector[next_j] % 4 == 1) {
                        base_vector[j] -= 2; // anihilate spin-down electron
                        base_vector[next_j] += 2; // create spin-down electron
                        idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                        base_vector[next_j] -= 2;
                        base_vector[j] += 2;
                        if ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 1) || ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 1)) || ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 0))) {
                            result_vec_threaded(k) += -t * initial_vec(idx);
                        }
                        else
                            result_vec_threaded(k) += t * initial_vec(idx);
                    }
                }
            }
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            if (base_vector[j] == 5) {// S_i^+ s_i^-
                base_vector[j] = 2;
                idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                base_vector[j] = 5;
                result_vec_threaded(k) += -J_H * initial_vec(idx);
            }
            if (base_vector[j] == 2) {// S_i^+ s_i^-
                base_vector[j] = 5;
                idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                base_vector[j] = 2;
                result_vec_threaded(k) += -J_H * initial_vec(idx);
            }
            // electron repulsion
            if (base_vector[j] == 7 || base_vector[j] == 3)
                result_vec_threaded(k) += U * initial_vec(k);
            //Diagonal - z part
            if (base_vector[j] == 1 || base_vector[j] == 6)
                result_vec_threaded(k) -= 2.0 * J_H * 0.25 * initial_vec(k);
            if (base_vector[j] == 2 || base_vector[j] == 5)
                result_vec_threaded(k) += 2.0 * J_H * 0.25 * initial_vec(k);
        }
    }
}
void Lanczos::Hamil_vector_multiply(vec& initial_vec, vec& result_vec) {
    result_vec.zeros();
    std::vector<arma::vec> result_threaded(num_of_threads);
    std::vector<std::thread> threads;
    //Hamil_vector_multiply_kernel(0, N, initial_vec, result_vec);
    threads.reserve(num_of_threads);
    for (int t = 0; t < num_of_threads; t++) {
        int start = t * N / num_of_threads;
        int stop = ((t + 1) == num_of_threads ? N : N * (t + 1) / num_of_threads);
        result_threaded[t] = arma::vec(stop - start, fill::zeros);
        threads.emplace_back(&Lanczos::Hamil_vector_multiply_kernel, this, start, stop, ref(initial_vec), ref(result_vec));
    }for (auto& t : threads) t.join();
}
void Lanczos::Lanczos_convergence(vec& initial_vec) {
    this->Krylov_space = mat(N, 1);
    this->H_L = mat(1, 1, fill::zeros);
    vec e_prev(5, fill::zeros);
    ofstream convergence("Lanczos convergence.txt");
    ofstream energy("Energy Lanczos.txt");
    convergence << std::setprecision(20) << std::fixed;

    Krylov_space.col(0) = initial_vec;

    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H_sparse * Krylov_space.col(0);

    double alfa = arma::dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        Krylov_space.resize(N, j + 1);
        H_L.resize(j + 1, j + 1);

        beta = sqrt(dot(tmp, tmp));
        Krylov_space.col(j) = tmp / beta;
        vec tmp2 = Krylov_space.col(j);
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H_sparse * tmp2;

        alfa = arma::dot(Krylov_space.col(j), tmp);
        //tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);
        vec temporary(N, fill::zeros);
        for (int k = 0; k <= j; k++)
            temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
        tmp = tmp - temporary;
        temporary.zeros();
        for (int k = 0; k <= j; k++)
            temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
        tmp = tmp - temporary;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;
        eig_sym(eigenvalues, eigenvectors, H_L);

        for (int k = 0; k < j; k++) {
            energy << j << "\t" << eigenvalues(k) << endl;
        }
        if (j >= 5) {
            vec e(e_prev.size());
            e(0) = eigenvalues(0);
            e(1) = eigenvalues(1);
            e(2) = eigenvalues(2);
            e(3) = eigenvalues(3);
            e(4) = eigenvalues(4);
            convergence << j << "\t";
            for (int k = 0; k < e.size(); k++) {
                convergence << fabs((e_prev(k) - e(k)) / e(k)) << "\t";
            }
            convergence << endl;
            e_prev = e;
        }
        out << j << endl;
    }
    energy.close();
    convergence.close();
}

void Lanczos::Build_Lanczos_Hamil_wKrylov(vec& initial_vec) {
    this->Krylov_space = mat(N, 1);
    this->H_L = mat(1, 1, fill::zeros);

    Krylov_space.col(0) = initial_vec;

    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H_sparse * Krylov_space.col(0);

    double alfa = arma::dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        try {
            Krylov_space.resize(N, j + 1);
            H_L.resize(j + 1, j + 1);
        }
        catch (const bad_alloc& e) {
            std::cout << "Memory exceeded " << e.what() << "\n";
            out << "dim(H) = " << Krylov_space.size() * sizeof(Krylov_space(0, 0)) * std::pow(10, -9) << " GB" << "\n";
            assert(false);
        }

        beta = sqrt(dot(tmp, tmp));
        Krylov_space.col(j) = tmp / beta;
        vec tmp2 = Krylov_space.col(j);
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H_sparse * tmp2;

        alfa = arma::dot(Krylov_space.col(j), tmp);
        if (use_reorthonormalization) {
            vec temporary(N, fill::zeros);
            for (int k = 0; k <= j; k++)
                temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
            tmp = tmp - temporary;
            temporary.zeros();
            for (int k = 0; k <= j; k++)
                temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
            tmp = tmp - temporary;
        }else
            tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;
    }
}
void Lanczos::Build_Lanczos_Hamil(vec& initial_vec) {
    this->H_L = mat(lanczos_steps, lanczos_steps, fill::zeros);

    double beta = arma::dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    // already normalized input random vector

    this->randVec_inKrylovSpace = vec(lanczos_steps);
    randVec_inKrylovSpace(0) = arma::dot(initial_vec, initial_vec); // =1

    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H_sparse * initial_vec;

    double alfa = arma::dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;
    //out << tmp.t();
    vec tmp2_prev = initial_vec;
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        double beta = sqrt(arma::dot(tmp, tmp));
        vec tmp2 = tmp / beta;
        randVec_inKrylovSpace(j) = dot(tmp2, initial_vec);

        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H_sparse * tmp2;

        alfa = arma::dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * tmp2_prev;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;

        tmp2_prev = tmp2;
        //out << j << "lanczos" << endl;
    }
}
void Lanczos::Lanczos_Diagonalization() {
    //Hamiltonian_sparse();
    //vec rand = randu<vec>(N);
    int averaged_iterations = 10;
    this->eigenvalues = vec(lanczos_steps, fill::zeros);
    this->eigenvectors = mat(lanczos_steps, lanczos_steps, fill::zeros);
    for (int r = 0; r < averaged_iterations; r++) {
        vec rand = Create_Random_vec();
        Build_Lanczos_Hamil(rand);
        vec E; 
        mat V;
        eig_sym(E, V, H_L);
        eigenvalues += E / averaged_iterations;
        eigenvectors += V / averaged_iterations;
    }
}
void Lanczos::Lanczos_GroundState() {

    if(!memory_over_performance) 
        Hamiltonian_sparse();

    auto initial_vec = Create_Random_vec();
    Build_Lanczos_Hamil(initial_vec);
    try {
        eig_sym(eigenvalues, eigenvectors, H_L);
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded " << e.what() << "\n";
        out << "dim(H) = " << H_L.size() * sizeof(H_L(0, 0)) << "\n";
        assert(false);
    }
    out << eigenvalues(0) << endl << endl;
    vec GS = eigenvectors.col(0);

    this->ground_state = vec(N, fill::zeros);

    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)

    vec tmp(N, fill::zeros);
    if(memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H_sparse * initial_vec;

    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    for (int j = 1; j < lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        ground_state += GS(j) * tmp2;
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H_sparse * tmp2;
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        initial_vec = tmp2;
        //out << j << "lanczos" << endl;
    }
}

vec Lanczos::Heat_Capacity_Lanczos(int random_steps) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec E_av(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec E_av2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->chi_0 = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->partition_function = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    //this->chi_0 = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    auto temperature = prepare_parameterVec(dT, T_end, dT);
    Hamiltonian_sparse();
    for (int r = 0; r < random_steps; r++) {
        auto rand_vec = Create_Random_vec();
        Build_Lanczos_Hamil(rand_vec);
        eig_sym(eigenvalues, eigenvectors, H_L);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
                overlap *= overlap;
                Z(k) += (double)N / (double)random_steps * overlap * std::exp(-eigenvalues(m) / T);
                E_av(k) += eigenvalues(m) * overlap * std::exp(-eigenvalues(m)/ T);
                E_av2(k) += eigenvalues(m) * eigenvalues(m) * overlap * std::exp(-eigenvalues(m) / T);
                partition_function(k) += Z(k);
                chi_0(k) += Sz * Sz * overlap * std::exp(-eigenvalues(m) / T) * (double)N / (double)random_steps / T / (double)L;
            }
        }
    }
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k];
        E_av(k) = E_av(k) / Z(k) * (double)N / (double)random_steps;
        E_av2(k) = E_av2(k) / Z(k) * (double)N / (double)random_steps;
        Cv(k) = (E_av2(k) - E_av(k) * E_av(k)) / T / T / (double)L;
    }
    return Cv;
}

vec Lanczos::thermal_average_lanczos(vec&& quantity, int& random_steps){
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
void Heat_Capacity(std::vector<arma::vec>&& energies, vec&& Cv) {
    auto temperature = prepare_parameterVec(dT, T_end, dT);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k];
        double heat_capacity = 0;
        int l = 0;
        for (auto E : energies) {
            double Sz = -((double)E.size() - 1.0) / 2. + (double)l;
            double E_av = 0;
            double E2_av = 0;
            double Z = 0;
            //out << Sz << endl;
            for (int n = 0; n < E.size(); n++) {
                E_av += E(n) * std::exp(-E(n) / T);
                E2_av += E(n) * E(n) * std::exp(-E(n) / T);
                Z += std::exp(-E(n) / T);
            }
            E_av = E_av / Z;
            E2_av = E2_av / Z;
            heat_capacity += (E2_av - E_av * E_av) / T / T;
            l++;
        }
        int L = (energies.size() - 1) * 2 / 3; // n=1.5 filling
        Cv(k) = heat_capacity / (double)L / (double)energies.size();

        //out << T << " " << Cv(k) << " " << L << endl;
    }
}
void static_spin_susceptibility(std::vector<arma::vec>&& energies, vec&& chi) {
    auto temperature = prepare_parameterVec(dT, T_end, dT);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k]; 
        double Z = 0;
        double X_0 = 0;
        int l = 0;
        for (auto E : energies) {
            double Sz = -((double)E.size() - 1.0) / 2. + (double)l;
            for (int n = 0; n < E.size(); n++) {
                Z += std::exp(-E(n)/ T);
                X_0 += Sz * Sz * std::exp(-E(n) / T);
            }
            l++;
        }
        int L = (energies.size() - 1) * 2 / 3; // n=1.5 filling
        chi(k) = X_0 / Z / T / (double)L;
        //out << T << " " << chi(k) << " " << L << endl;
    }
}
vec Sq_lanczos(int random_steps, double T, std::unique_ptr<Lanczos>& obj) {
    vec Sq(obj->L + 1, fill::zeros);
    if(memory_over_performance)
        obj->Hamiltonian_sparse();
    int number_of_thr = (obj->L > 10) ? 1 : ((obj->L < 10) ? obj->L + 1 : 2); // if - elseif - else statement
#pragma omp parallel for num_threads(number_of_thr)
    for (int l = 0; l <= obj->L; l++) {
        vector<int> vect(obj->L);
        double q = (double)l * pi / ((double)obj->L + 1.0);
        sp_cx_mat Sq_mat(sp_mat(obj->N, obj->N), sp_mat(obj->N, obj->N));
        for (int p = 0; p < obj->N; p++) {
            int_to_binary(obj->mapping->at(p), vect);
            Sq_mat(p, p) = 0;
            for (int m = 0; m < obj->L; m++) {
                double Szm = 0;
                if (vect[m] < 4) Szm += 0.5;
                else Szm -= 0.5;
                if (vect[m] % 4 == 1) Szm += 0.5;
                else if (vect[m] % 4 == 2) Szm -= 0.5;
                Sq_mat(p, p) += std::exp(1i * q * (m + 0.0)) * Szm;
            }
        }
        double Sq0 = 0, Z = 0;
        std::unique_ptr<Lanczos> obj_copy(new Lanczos(obj));
        obj_copy->H_sparse = obj->H_sparse;
        for (int r = 0; r < random_steps; r++) {
            auto rand_vec = obj_copy->Create_Random_vec();
            obj_copy->Build_Lanczos_Hamil_wKrylov(rand_vec);
            mat V;
            eig_sym(obj_copy->eigenvalues, V, obj_copy->H_L);
            obj_copy->eigenvectors = obj_copy->Krylov_space * V;
            cx_vec temp = Sq_mat * Sq_mat.t() * rand_vec;
            for (int m = 0; m < obj_copy->lanczos_steps; m++) {
                double overlap = dot(rand_vec, obj_copy->eigenvectors.col(m));
                Z += (double)obj_copy->N / (double)random_steps * overlap * overlap * std::exp(-obj_copy->eigenvalues(m) / T);
                overlap = real(overlap * cdot(cx_vec(obj_copy->eigenvectors.col(m), vec(obj_copy->N, fill::zeros)), temp));
                Sq0 += (double)obj_copy->N / (double)random_steps * overlap * std::exp(-obj_copy->eigenvalues(m) / T);
            }
        }
        Sq(l) = real(2.0 * Sq0 / pi / (obj_copy->L + 1.0));
        obj->Z_constT = Z;
    }
    return Sq;
}

void Cv_Umap(int L, int N_e, double t) {
    ofstream savefile;
    stringstream Nstr;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("C_V_L=" + std::to_string(L) + "_Umap.txt");
    for (int k = 2; k < 40; k += 1) {
        double U = (double)k / 10.0;
        double K, J_H;
        K = 4 * 0.15 * 0.15 / U;
        J_H = 0.25 * U;
        std::vector<arma::vec> energies;
        //std::unique_ptr<std::vector<arma::vec>> energies(new std::vector<arma::vec>);
        for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
            //HamiltonianKH Hamil(L, N_e, t, U, K, J_H, Sz);
            std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, Sz));
            Hamil->Hamiltonian();
            Hamil->Diagonalization();
            energies.emplace_back(Hamil->eigenvalues);
        }
        vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        Heat_Capacity(std::move(energies), std::move(Cv));
        int q = 0;
        double T = dT;
        while (T <= T_end) {
            savefile << T << "\t\t" << U << "\t\t" << Cv(q) << endl; //save heat capacity to file
            T += dT; q++;
        }
        //Main_DOS(L, N_e, t, K, U, J_H);
        //Main_Cv(L, N_e, t, K, U, J_H);

        out << "U = " << U << " done!" << endl;
    }
    savefile.close();
}
void DOS_Umap(int L, int N_e, double t) {
    ofstream savefile;
    stringstream Nstr;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("DOS_L=" + std::to_string(L) + "_Umap.txt");
    for (int k = 1; k < 80; k += 1) {
        double U = (double)k / 10.0;
        double K, J_H;
        K = 4 * 0.15 * 0.15 / U;
        J_H = 0.25 * U;

        std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1));
        Hamil->Hamiltonian();
        Hamil->Diagonalization();
        vector<double> omega_vec = prepare_parameterVec(Hamil->eigenvalues(0) - 100 * domega, Hamil->eigenvalues(Hamil->eigenvalues.size() - 1) + 100 * domega, domega);
        vec DOS = Hamil->Total_Density_of_states(std::move(omega_vec));
        int q = 0;
        double Ef = FermiLevel(L, N_e, t, K, U, J_H);
        for(int q = 0; q < omega_vec.size(); q++){
            savefile << omega_vec[q] - Ef << "\t\t" << U << "\t\t" << DOS(q) << endl; //save heat capacity to file
        }

        out << "U = " << U << " done!" << endl;
    }
    savefile.close();
}

void Main_U(int L, int N_e, double t) {
//#pragma omp parallel for 
    for (int k = 1; k <= 50; k += 1) {
        double U = (double)k / 10.0;
        double K, J_H;
        K = 4 * 0.15 * 0.15 / U;
        J_H = 0.25 * U;
        //K = 0; J_H = 0;
        Main_DOS(L, N_e, t, K, U, J_H);
        Main_Cv(L, N_e, t, K, U, J_H);

        out << "U = " << U << " done!" << endl;
    }
}
void Main_Jh(int L, int N_e, double t, double K, double U) {
    double J_H = 0.01 * U;
    while(J_H <= 0.3*U) {
        Main_DOS(L, N_e, t, K, U, J_H);
        Main_Cv(L, N_e, t, K, U, J_H);
        out << "J_H/U = " << J_H / U << " done!" << endl;
        J_H += 0.01 * U;
    }
}

void Main_Cv(int L, int N_e, double t, double K, double U, double J_H) {
    std::vector<arma::vec> energies;
    vec Sq(L + 1, fill::zeros);
    double Z = 0;
    for (int Sz = -N_e; Sz <= N_e; Sz +=2) {
        std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, Sz));
        Hamil->Diagonalization();
        Sq += Hamil->static_structure_factor(T);
        Z += Hamil->partition_function(T);
        energies.emplace_back(std::move(Hamil->eigenvalues));
        //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
    }
    Sq = Sq / Z;
    print_Sq(std::move(Sq), U, N_e, L, T);

    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    Heat_Capacity(std::move(energies),  std::move(Cv));
    print_Cv(std::move(Cv), U, N_e, L);

    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    static_spin_susceptibility(std::move(energies), std::move(chi_0));
    print_chi(std::move(chi_0), U, N_e, L);
}
void Main_Sq(int L, int N_e, double t, double K, double U, double J_H) {
    vec Sq(L + 1, fill::zeros);
    double Z = 0;
    double T = 1.0;
    for (double T = 0.01; T <= 1; T += 0.02) {
        for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
            std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, Sz));
            Hamil->Diagonalization();
            Sq += Hamil->static_structure_factor(T);
            Z += Hamil->partition_function(T);
            //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
        }
        Sq = Sq / Z;
        print_Sq(std::move(Sq), U, N_e, L, T);
    }
}
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H) {
    //HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
    std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1));
    Hamil->Diagonalization();
    vector<double> omega_vec = prepare_parameterVec(Hamil->eigenvalues(0) - 100 * domega, Hamil->eigenvalues(Hamil->eigenvalues.size() - 1) + 100 * domega, domega);
    vec DOS = Hamil->Total_Density_of_states(std::move(omega_vec)); //rvalue of DOS
    double maximum = max(DOS);
    double Ef = FermiLevel(L, N_e, t, K, U, J_H);
    //Hamil->printEnergy(Ef);
    printDOS(std::move(DOS), U, N_e, L, std::move(omega_vec), maximum, Ef);
}
void Main_Cv_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps) {
    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Sq(L + 1, fill::zeros);
    double Z_constT = 0;
    for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
        //Lanczos Hamil(L, N_e, t, U, K, J_H, Sz, M);
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, Sz, M));
        if (M > Hamil->N) {
            Hamil->lanczos_steps = Hamil->N;
            //out << "Lanczos steps greater than system size" << endl;
        }
        Cv += Hamil->Heat_Capacity_Lanczos(random_steps) / double(N_e + 1.0);
        chi_0 += Hamil->chi_0;
        Z += Hamil->partition_function;
        Sq += Sq_lanczos(random_steps, T, Hamil);
        Z_constT += Hamil->Z_constT;
        //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
    }
    /*std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0)? 0:1, M));
    Sq = Sq_lanczos(random_steps, T, Hamil);
    Z_constT = Hamil->Z_constT;*/
    Sq = Sq / Z_constT;
    print_Sq(std::move(Sq), U, N_e, L, T);

    chi_0 = chi_0 / Z; // elementwise division
    print_Cv_Lanczos(std::move(Cv), U, N_e, L, M, random_steps);
    print_chi(std::move(chi_0), U, N_e, L);
}


// Helpful tools
double FermiLevel(int L, int N_e, double t, double K, double U, double J_H) {
    double Ef;
    std::unique_ptr<HamiltonianKH> Object(new HamiltonianKH(L, N_e + 1, t, U, K, J_H, ((N_e + 1) % 2 == 0) ? 0 : 1));
    Object->Hamiltonian();
    Object->Diagonalization();
    Ef = Object->eigenvalues(0);

    Object.reset(new HamiltonianKH(L, N_e - 1, t, U, K, J_H, ((N_e - 1) % 2 == 0) ? 0 : 1));
    Object->Hamiltonian();
    Object->Diagonalization();
    Ef = (Ef + Object->eigenvalues(0)) / 2.0;

    return Ef;
}
std::vector<double> prepare_parameterVec(double _min, double _max, double step) {
    vector<double> parameter_vec;
    double temp = _min;
    while (temp <= _max) {
        parameter_vec.emplace_back(temp);
        temp += step;
    }
    return parameter_vec;
}

void printDOS(vec&& resultDOS, double U, double N_e, int L, std::vector<double>&& omega_vec, double maximum, double E_fermi) {
    ofstream DOSfile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    DOSfile.open("DOS_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
    //DOSfile.open("DOS_dw=" + std::to_string(domega) + "_eta=" + std::to_string(eta) + ".txt");
    for (int k = 0; k < omega_vec.size(); k++)
        DOSfile << omega_vec[k] - E_fermi << "\t\t" << resultDOS[k] / maximum << endl;

    DOSfile.close();
}
void print_Cv(vec&& Cv, double U, double N_e, int L) {
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
void print_Cv_Lanczos(vec&& Cv, double U, double N_e, int L, int M, int random_steps) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    //savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + "_M=" + std::to_string(M) + "_R=" + std::to_string(random_steps) + ".txt");
    savefile.open("Cv_Lanczos_L=" + std::to_string(L) + "_gKH.txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << Cv(k) << endl; //save heat capacity to file
        T += dT; k++;
    }
    savefile.close();
}
void print_chi(vec&& chi, double U, double N_e, int L) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("chi_0_L=" + std::to_string(L) + "_U=" + Ustr.str() + "_L.txt");
    int k = 0;
    double T = dT;
    while (T <= T_end) {
        savefile << T << "\t\t" << chi(k) << endl;
        T += dT; k++;
    }
    savefile.close();
}
void print_Sq(vec&& Sq, double U, double N_e, int L, double T) {
    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    savefile.open("Sq_L=" + std::to_string(L) + "_U=" + Ustr.str() + "_T=" + std::to_string(T) + ".txt");
    savefile << std::setprecision(8) << fixed;
    for (int l = 0; l <= L; l++) {
        double q = (double)l * pi / ((double)L + 1.0);
        savefile << q << "\t\t" << Sq(l) << endl;
        out << q / pi << "\t" << Sq(l) << endl;
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
        Efile << U << "\t\t" << eigenvalues(k) - Ef << endl;

    Efile.close();

}

// Conversion of int to binary vector - using modulo operator
void int_to_binary(ull_int idx, std::vector<int> & vec) {
    ull_int temp = idx;
    for (int k = 0; k < vec.size(); k++) {
        vec[vec.size() - 1 - k] = static_cast<int>(temp % 8);
        temp = static_cast<ull_int>((double)temp / 8.);
    }
}
ull_int binary_to_int(vector<int>& vec) {
    ull_int val = 0;
    for (int k = 0; k < vec.size(); k++) {
        val += vec[vec.size() - 1 - k] * std::pow(8, k);
    }
    return val;
}

ull_int binary_search(std::vector<ull_int>* arr, int l_point, int r_point, ull_int element) {
    //out << element << endl;
    if (r_point >= l_point) {
        int middle = l_point + (r_point - l_point) / 2;
        if (arr->at(middle) == element) return middle;
        else if (arr->at(middle) > element) return binary_search(arr, l_point, middle - 1, element);
        else return binary_search(arr, middle + 1, r_point, element);
    }
    out << "Element " << element << " not present in the array" << endl;
    assert(false);
    return -1;
}

void Lanczos_convergence(int L, int N_e) {
    double t = 0.5, U = 2.1, K = 4 * 0.15 * 0.15 / U, J_H = 0.25 * U;
    std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, 500));
    Hamil->Hamiltonian_sparse();
    vec random_vec = Hamil->Create_Random_vec();
    Hamil->Lanczos_convergence(random_vec);
}