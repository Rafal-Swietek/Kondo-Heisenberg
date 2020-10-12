#include "Hamiltonian.h"

#define im std::complex<double>(0.0,1.0)
#define M_PI 3.14159265358979323846

double pi = M_PI;

typedef std::complex<double> cpx;
typedef std::vector<double> vect;

void myassert(std::string text){
    std::cerr << text << std::endl;
    assert(false);
}

HamiltonianKH::~HamiltonianKH() {
	H.~Mat();
	eigenvectors.~Mat();
	eigenvalues.~vec();

    H_L.~Mat();
    eigenVal_L.~vec();
    Krylov_space.~Mat();
}
//-----------------------

//--------------------------------------------------------------------------------------------------
//----------------------------------------SUBBLOCKS HAMILTONIAN-------------------------------------
//------------------------------------------- S_z = 0 block-----------------------------------------
//--------------------------------------------------------------------------------------------------
// Constructor for Hamiltonian separated to blocks
HamiltonianKH::HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H) {
	this->L = L; //number of sites
	this->num_of_electrons = num_of_electrons; //number of electrons in lower orbital
	this->t = t; this->U = U; 
	this->K = K;
	this->J_H = J_H;

	this->mapping = vector<int>(std::pow(8, L));
	this->mapping_inv = vector<int>(); 
	generate_mapping_subblock();
	this->N = mapping_inv.size();

	this->H = mat(N, N, arma::fill::zeros); //hamiltonian
}
//------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------TOTAL HAMILTONIAN-------------------------------------
//----------------------------------------------------------------------------------------------
// Constructor for total Hamiltonian (no subblocks)
HamiltonianKH::HamiltonianKH(int L, double t, double U, double K, double J_H) {
	this->L = L; //number of sites
	this->N = std::pow(8, L);
	this->t = t; this->K = K;
	this->J_H = J_H; this->U = U;

	this->H = mat(N, N, arma::fill::zeros); //hamiltonian

	this->mapping = vector<int>(N);
	this->mapping_inv = vector<int>(N);
	generate_mapping_total();
}
//-------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------BUILD HAMILTONIAN-------------------------------------
//----------------------------------------------------------------------------------------------

void HamiltonianKH::setHamiltonianElem(int k, double value, std::vector<int> temp){
    int idx = mapping[binary_to_int(temp)];
        H(idx, k) += value;
        H(k, idx) += value;
}
void HamiltonianKH::Hamiltonian() {
    int s_i, s_j; //i=j, j=j+1
	bool PBC = 0; //allows periodic boundary conditions if =1
	int next_j;
    for (int k = 0; k < N; k++){
		std::vector<int> base_vector = int_to_binary(mapping_inv[k], L);
		vector<int> temp(base_vector);
		for (int j = 0; j <= L - 1; j++) {
            if (PBC == 1 && j == L - 1) next_j = 0;
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
//----------------------------------------------------

//generates the vector, which maps the base_vector index to the index in given subblock
void HamiltonianKH::CreateMappingElem(int& bSz, int& fSz, int& N_e, int& j, int& idx) {
    mapping[j] = -1;
    if ((bSz + fSz == 0) && N_e == num_of_electrons) {
        mapping_inv.push_back(j);
        mapping[j] = idx;
        idx++;
    }
}
std::tuple<double, double, int> calculateSpinElements(int L, int j) {
    double bSz = 0; //bosonic total spin - spin of upper orbital locked to n=1 filling
    double fSz = 0; //fermionic total spin
    int N_e = 0; // numer of electrons in given state
    vector<int> temp = int_to_binary(j, L);

    for (int k = 0; k < L; k++) {
        if (temp[k] < 4) bSz += 0.5;
        else bSz -= 0.5;
        if (temp[k] % 4 == 1) {
            fSz += 0.5;
            N_e += 1;
        }
        else if (temp[k] % 4 == 2) {
            fSz -= 0.5;
            N_e += 1;
        }
        else if (temp[k] % 4 == 3)
            N_e += 2;
    }

    return std::make_tuple(bSz, fSz, N_e);
}
void HamiltonianKH::generate_mapping_subblock() {
	int idx = 0;
	for (int j = 0; j < std::pow(8, L); j++) {
        int bSz, fSz, N_e;
        std::tie(bSz,fSz,N_e) = calculateSpinElements(L,j);
        CreateMappingElem(bSz,fSz, N_e,j,idx);
    }
    assert(idx>0 && "Not possible number of electrons - no. of states < 1");
}

void HamiltonianKH::generate_mapping_total() {
	for (int j = 0; j < N; j++) { mapping[j] = j; mapping_inv[j] = j; }
}
//----------------------------------------------------

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
//----------------------------------------------------

//----------------------------------------------------------------------------------------------
//--------------------------------------Diagonalization & Quantities-----------------------------
//----------------------------------------------------------------------------------------------

void HamiltonianKH::Diagonalization() {
	try {
		arma::eig_sym(eigenvalues, H);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		std::cout << H.size() * sizeof(H(0, 0)) << "\n";
        	assert(false);
    	}
}
//----------------------------------------------------

// Helpful tools
void HamiltonianKH::printEnergy() {
    ofstream Efile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)num_of_electrons / (double)L;
    Efile.open("E_n=" + Nstr.str() + "_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < N; k++)
        Efile << U << "\t\t" << eigenvalues(k) - eigenvalues(0) << endl;

    Efile.close();

}

std::vector<double> prepareOmegaVec(arma::vec &eigenvalues, double dOmega){
    vector<double> omega_vec;
    double omega = eigenvalues(0);
    while (omega <= eigenvalues(eigenvalues.size() - 1)) {
        omega_vec.push_back(omega);
        omega += dOmega;
    }
    return omega_vec;
}
//-----------

// Calculates the density of states using one-particle greens function
void printDOS(vect resultDOS, double U, double N_e, int L, vect omega_vec, double maximum) {
    ofstream DOSfile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
    DOSfile.open("DOS_n=" + Nstr.str() + "_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < omega_vec.size(); k++)
        DOSfile << omega_vec[k] << "\t\t" << resultDOS[k] / maximum << endl;

    DOSfile.close();
}
void HamiltonianKH::Density_of_states(int N_e) {
	double domega = 0.005;

    vector<double> omega_vec = prepareOmegaVec(eigenvalues,domega);
    vector<double> resultDOS(omega_vec.size());

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

        resultDOS[w] = DOS;
    }
    printDOS(resultDOS,U,N_e,L,omega_vec,maximum);
}

void HamiltonianKH::Heat_Capacity() {
    double dT = 0.05;
    double T = 0.05;
    double energy_av; //average of energy E
    double energy2_av; //average of E^2

    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(2) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)num_of_electrons / (double)L;
    //savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + ".txt");
    savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + ".txt");

    while (T <= 5.0) {
        double Partition_Function = 0;
        energy_av = 0; energy2_av = 0;
        for (int j = 0; j < N; j++) {
            Partition_Function += std::exp(-eigenvalues(j) / T); //partition function(T)
            energy_av += eigenvalues(j) * std::exp(-eigenvalues(j) / T); //average energy(T)
            energy2_av += eigenvalues(j) * eigenvalues(j) * std::exp(-eigenvalues(j) / T); //average energy^2(T)
        }
        energy_av = energy_av / Partition_Function;
        energy2_av = energy2_av / Partition_Function;
        double heat_capacity = (energy2_av - energy_av * energy_av) / T / T / (L + 0.0);
        savefile << T << "\t\t" << heat_capacity << endl; //save heat capacity to file
        T += dT;
    }
    savefile.close();
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

//----------------------------------------------------------------------------------------------
//--------------------------------------------------LANCZOS-------------------------------------
//----------------------------------------------------------------------------------------------
vec HamiltonianKH::Create_Random_vec() {
    vec random_vec(N, fill::zeros);
    for (int j = 0; j < N; j++)
        random_vec(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
    return random_vec;
}


vec HamiltonianKH::Hamil_vector_multiply(vec initial_vec) {
    vec result_vec(N, fill::zeros);
    bool PBC = 0;
    int next_j, idx;
    int s_i, s_j;
    for (int k = 0; k < N; k++) {
        std::vector<int> base_vector = int_to_binary(mapping_inv[k], L);
        vector<int> temp(base_vector);
        for (int j = 0; j <= L - 1; j++) {
            if (PBC == 1 && j == L - 1) next_j = 0;
            else if (PBC == 0 && j == L - 1) goto kinetic_term_omitted2;
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
                idx = mapping[binary_to_int(temp)];
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
                idx = mapping[binary_to_int(temp)];
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
                    idx = mapping[binary_to_int(temp)];
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
        kinetic_term_omitted2:
            // electron repulsion
            if (base_vector[j] == 7 || base_vector[j] == 3)
                result_vec(k) += U * initial_vec(k);
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            temp = base_vector;
            if (base_vector[j] == 5) {// S_i^+ s_i^-
                temp[j] = 2;
                idx = mapping[binary_to_int(temp)];
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
void HamiltonianKH::Build_Lanczos_Hamil_wKrylovSpace(vec initial_vec, int Lanczos_steps) {
    this->H_L = mat(Lanczos_steps, Lanczos_steps, fill::zeros);
    //this->Lanczos_GS = vec(N);
    this->Krylov_space = mat(N, Lanczos_steps);

    Krylov_space.col(0) = initial_vec;
    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); // normalized Krylov_space(j=0)

    vec tmp = Hamil_vector_multiply(Krylov_space.col(0)); // tmp = H * Krylov_space(0)
    double alfa = dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);

    H_L(0, 0) = alfa;
    for (int j = 1; j < Lanczos_steps; j++) {
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
void HamiltonianKH::Build_Lanczos_Hamil(vec initial_vec, int Lanczos_steps) {
    this->H_L = mat(Lanczos_steps, Lanczos_steps, fill::zeros);
    
    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    
    vec tmp = Hamil_vector_multiply(initial_vec); // tmp = H * Krylov_space(0)
    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    H_L(0, 0) = alfa;
    for (int j = 1; j < Lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        tmp = Hamil_vector_multiply(tmp2); // tmp = H * tmp2
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;

        initial_vec = tmp2;
    }
    tmp.~vec();
}
void HamiltonianKH::Lanczos_Diagonalization(int lanczos_steps) {
    srand(time(NULL));
    Build_Lanczos_Hamil(Create_Random_vec(), lanczos_steps);
    this->eigenVal_L = vec(lanczos_steps);
    eig_sym(eigenVal_L, H_L);
}

double HamiltonianKH::Cv_kernel(double T) {
    int random_cycles = 25; // number of random cycles to compute heat capacity
    int Lancz_steps = 100;
    double Z = 0; //partition function
    double E_av = 0; // average energy
    double E2_av = 0; // average squared energy
    double overlap = 0; //overlap of random vectopr and lanczos eigenvectors

    for (int r = 0; r < random_cycles; r++) {
        vec rand_vec = Create_Random_vec();
        Build_Lanczos_Hamil_wKrylovSpace(rand_vec, Lancz_steps);
        mat eigenVec;
        eig_sym(eigenVal_L, eigenVec, H_L);
        mat Lanczos_eigenvec(N, Lancz_steps);
        Lanczos_eigenvec = Krylov_space * eigenVec;
        eigenVec.~Mat();

        for (int m = 0; m < Lancz_steps; m++) {
            overlap = dot(rand_vec, Lanczos_eigenvec.col(m));
            overlap *= overlap;
            Z += (double)N / (double)random_cycles * overlap * std::exp(-eigenVal_L(m) / T);
            E_av += eigenVal_L(m) * overlap * std::exp(-eigenVal_L(m) / T);
            E2_av += eigenVal_L(m) * eigenVal_L(m) * overlap * std::exp(-eigenVal_L(m) / T);
        }
        Lanczos_eigenvec.~Mat();
    }
    E_av = E_av / Z * (double)N / (double)random_cycles;
    E2_av = E2_av / Z * (double)N / (double)random_cycles;

    return (E2_av - E_av * E_av) / T / T / (L + 0.0);
}
void HamiltonianKH::Heat_Capacity_Lanczos() {
    double dT = 0.05, T = dT;

    ofstream savefile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)num_of_electrons / (double)L;
    savefile.open("C_V_n=" + Nstr.str() + "_U=" + Ustr.str() + "_M=.txt");

    srand(time(NULL));
    while (T <= 4.0) {
        double Cv = Cv_kernel(T);
        savefile << T << "\t\t" << Cv << endl;
        T += dT;
        out << T << "\t\t" << Cv << endl;
    }
    savefile.close();
}



//----------------------------------------------------------------------------------------------
//----------------------------------------Main subroutines---------------------------------------
//----------------------------------------------------------------------------------------------

void Main_U(int L, int N_e, double t) {
	for (double U = 0.2; U < 3.0; U += 0.2) {
		double K, J_H;
		K = 4 * 0.15 * 0.15 / U;
		J_H = 0.25 * U;
		HamiltonianKH Object(L, N_e, t, U, K, J_H);
		Object.Hamiltonian();
		Object.Diagonalization();
		Object.Density_of_states(N_e);
        Object.printEnergy();
        //Object.Heat_Capacity();

        out << "U = " << U << " done!" << endl;
	}
}

void Main_Jh(int L, int N_e, double t, double K, double U) {
    double J_H = 0.05 * U;
    while(J_H <= 0.3*U) {
        HamiltonianKH Object(L, N_e, t, U, K, J_H);
        Object.Hamiltonian();
        Object.Diagonalization();
        //Object.Density_of_states(N_e);
        Object.Heat_Capacity();

        out << "J_H/U = " << J_H / U << " done!" << endl;
        J_H += 0.05 * U;
    }
}

