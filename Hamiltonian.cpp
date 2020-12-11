#include "Hamiltonian.h"

double pi = 3.14159265358979323846;
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

    this->mapping = std::make_unique<std::vector<ull_int>>();
    //this->mapping = std::vector<ull_int>();
	generate_mapping();
	this->N = mapping->size();
    Hamiltonian();
}
HamiltonianKH::HamiltonianKH() {}
//-------------------------

void HamiltonianKH::update_parameters(double t, double U, double K, double J_H, double Sz) {
    this->t = t; this->U = U, this->K = K; this->J_H = J_H;
    this->Sz = Sz;
    this->mapping.reset(new vector<ull_int>());
    //this->mapping = std::vector<ull_int>();
    generate_mapping();
    this->N = mapping->size();
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

ull_int findElement(std::vector<ull_int>* vector, ull_int element) {
    std::vector<ull_int>::iterator it = find(vector->begin(), vector->end(), element);
    assert(it != vector->end() && "Element not present in the array");
    return (ull_int)std::distance(vector->begin(), it);
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
void HamiltonianKH::mapping_kernel(ull_int start, ull_int stop, my_uniq_ptr& map_threaded, int _id) {
    int n = 1;
    //out << "A new thread joined tha party! from " << start << " to " << stop << endl;
    std::vector<int> temp(L);
    for (ull_int j = start; j < stop; j++) {
        int bSz, fSz, N_e;
        std::tie(bSz, fSz, N_e) = calculateSpinElements(this->L, j, temp);
        if ((bSz + fSz == this->Sz) && N_e == this->num_of_electrons) 
             map_threaded->push_back(j);
    }
}
void HamiltonianKH::generate_mapping() {
    ull_int start = 0, stop = (ull_int)std::pow(8, L);
    //mapping_kernel(start, stop, mapping, L, Sz, num_of_electrons);
    //Threaded
    //std::vector<my_uniq_ptr> map_threaded(num_of_threads);
    std::vector<my_uniq_ptr> map_threaded(num_of_threads);
    std::vector<std::thread> threads;
    threads.reserve(num_of_threads);
    for (int t = 0; t < num_of_threads; t++) {
        start = t * (ull_int)std::pow(8, L) / num_of_threads;
        stop = ((t + 1) == num_of_threads ? (ull_int)std::pow(8, L) : ull_int(std::pow(8, L) / (double)num_of_threads * (double)(t + 1) ));
        //map_threaded[t] = my_uniq_ptr(new std::vector<ull_int>());
        map_threaded[t] = std::make_unique<std::vector<ull_int>>();
        //map_threaded[t]->reserve(std::pow(2, L)); ??
        threads.emplace_back(&HamiltonianKH::mapping_kernel, this, start, stop, ref(map_threaded[t]), t);
    }
    for (auto& t : threads) t.join();

    for (auto& t : map_threaded)
        mapping->insert(mapping->end(), std::make_move_iterator(t->begin()), std::make_move_iterator(t->end()));
    mapping->shrink_to_fit();
    //sort(mapping->begin(), mapping->end());
    if (show_system_size_parameters) {
        out << "Mapping generated with  " << mapping->size() << "  elements" << endl;
        out << "Last element = " << mapping->at(mapping->size() - 1) << endl;
    }
    //out << mapping[0] << " " << mapping[mapping->size() - 1] << endl;
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
#pragma omp parallel for shared(omega_vec, resultDOS) num_threads(num_of_threads)
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
                Sq(p, p) += std::exp(cpx(1i * q * (m + 0.0))) * Szm;
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
            val += (1 - int(base_vector[base_vector.size() - 1 - j] / 4)) * (int)std::pow(2, j);
        GS(val) += ground_state(k) * ground_state(k);
    }
    //GS = arma::abs(GS);
    GS = GS / arma::dot(GS, GS); //normalizing to be sure
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

