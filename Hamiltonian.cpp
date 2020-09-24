#include "Hamiltonian.h"


mutex my_mutex;
complex<double> i = 1i;
double pi = 3.141592653;

//Destructor
HamiltonianKH::~HamiltonianKH() {
	H.~Mat();
	eigenvectors.~Mat();
	/*H_Lanczos.~Mat();
	Krylov_space.~Mat();
	Lanczos_eigenVal.~vec();*/

	eigenvalues.~vec();
    //mapping.~vector();
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

	this->base_vector = vector<int>(L);
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

	this->base_vector = vector<int>(L); //spin on each site, for instance |0100> for L=4
	this->H = mat(N, N, arma::fill::zeros); //hamiltonian

	this->mapping = vector<int>(N);
	this->mapping_inv = vector<int>(N);
	generate_mapping_total();
}
//-------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------BUILD HAMILTONIAN-------------------------------------
//----------------------------------------------------------------------------------------------

// Generating Hamiltonian
void HamiltonianKH::Hamiltonian() {
	int s_i, s_j; //i=j, j=j+1
	int idx = 0; //indices equivalent to spin-filp due to kinetic term
	bool PBC = 0; //allows periodic boundary conditions if =1
	int next_j;
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(mapping_inv[k], L);
		vector<int> temp(base_vector);
		for (int j = 0; j <= L - 1; j++) {
			if (PBC == 1 && j == L - 1) next_j = 0;
			else if (PBC == 0 && j == L - 1) goto kinetic_term_omitted;
			else next_j = j + 1;
			// Diagonal spin part
			// i
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
				idx = mapping[binary_to_int(temp)];
				H(idx, k) += K / 2.;
			}
			temp = base_vector;
			if (s_i == 1 && s_j == 0) { // S_i^- S_i+1^+
				temp[j] = base_vector[j] + 4; //spin filp
				temp[next_j] = base_vector[next_j] - 4;
				idx = mapping[binary_to_int(temp)];
				H(idx, k) += K / 2.;
			}
			//---------------------
			// electron hopping
				//spin up
						//j+1 -> j
				temp = base_vector;
				//only odd numbers have up-electrons  //even numbers lack one up-electron
				if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
					temp[next_j] -= 1; // anihilate spin-up electron
					temp[j] += 1; // create spin-up electron
					idx = mapping[binary_to_int(temp)];
					H(idx, k) += t;
				}
				//j -> j+1
				temp = base_vector;
				if (base_vector[j] % 2 == 1 && base_vector[next_j] % 2 == 0) {
					temp[j] -= 1; // anihilate spin-up electron
					temp[next_j] += 1; // create spin-up electron
					idx = mapping[binary_to_int(temp)];
					H(idx, k) += t;
				}
			//spin down
				temp = base_vector;
				// the if statement contains every possible down-electron hopping: next_j->j
				if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
					if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
						temp[next_j] -= 2; // anihilate spin-down electron
						temp[j] += 2; // create spin-down electron 
						idx = mapping[binary_to_int(temp)];
						if (base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 0)
							H(idx, k) -= t;
						else if (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)
							H(idx, k) -= t;
						else
							H(idx, k) += t;
					}
				}
				//j -> j+1
				temp = base_vector;
				if (base_vector[j] % 4 == 2 || base_vector[j] % 4 == 3) {
					if (base_vector[next_j] % 4 == 0 || base_vector[next_j] % 4 == 1) {
						temp[j] -= 2; // anihilate spin-down electron
						temp[next_j] += 2; // create spin-down electron 
						idx = mapping[binary_to_int(temp)];
						if (base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 0)
							H(idx, k) -= t;
						else if (base_vector[next_j] % 4 == 1 && base_vector[j] % 4 == 2)
							H(idx, k) -= t;
						else
							H(idx, k) += t;
					}
				}
			kinetic_term_omitted:
				//---------------------
			// electron repulsion
				if (base_vector[j] == 7 || base_vector[j] == 3) 
					H(k, k) += U;
			//--------------------
			// electron-localised spin interaction ( interorbital electronic spin interaction)
				temp = base_vector;
				if (base_vector[j] == 5) {// S_i^+ s_i^-
					temp[j] = 2;
					idx = mapping[binary_to_int(temp)];
					H(idx, k) -= J_H;
				}
				if (base_vector[j] == 2) {// S_i^-s_i^+
					temp[j] = 5;
					idx = mapping[binary_to_int(temp)];
					H(idx, k) -= J_H;
				}
				//Diagonal - z part
				if (base_vector[j] == 1 || base_vector[j] == 6)
					H(k, k) -= 2.0 * J_H * 0.25;
				if (base_vector[j] == 2 || base_vector[j] == 5)
					H(k, k) += 2.0 * J_H * 0.25;
				//--------------------
		}
        //temp.~vector();
	}
}
//----------------------------------------------------

//generates the vector, which maps the base_vector index to the index in given subblock
void HamiltonianKH::generate_mapping_subblock() {
	int idx = 0;
	for (int j = 0; j < std::pow(8, L); j++) {
		vector<int> temp = int_to_binary(j, L);
		int bSz = 0; //bosonic total spin - spin of upper orbital locked to n=1 filling
		int fSz = 0; //fermionic total spin
		int N_e = 0; // numer of electrons in given state
		mapping[j] = -1;
		for (int k = 0; k < L; k++) {
			if (temp[k] < 4) bSz += 1;
			else bSz -= 1;
			// if temp[k] % 4 == 0 then fSz += 0 and N_e += 0
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
		if ((bSz + fSz == 0) && N_e == num_of_electrons) {
			mapping_inv.push_back(j);
            mapping[j] = idx;
			idx++;
        }
    }
	if (idx < 1) {
		cout << "Not possible number of electrons - no. of states < 1" << endl;
		exit(-3);
	}
}
//----------------------------------------------------

//generates the vector, which maps for the total hamiltonian: map(i) = i
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
// Conversion of binary vector to int
int binary_to_int(vector<int> vec) {
	int val = 0;
	for (int k = 0; k < vec.size(); k++) {
		val += vec[vec.size() - 1 - k] * std::pow(8, k);
	}
	return val;
}
//----------------------------------------------------

//----------------------------------------------------------------------------------------------
//--------------------------------------Diagonalization & exercises-----------------------------
//----------------------------------------------------------------------------------------------

//Diagonalizes the hamiltonian
void HamiltonianKH::Diagonalization() {
	this->eigenvalues = vec(N); //eigenvalues
	this->eigenvectors = mat(N, N); //eigenvectors
	try {
		arma::eig_sym(eigenvalues, H);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		std::cout << H.size() * sizeof(H(0, 0)) << "\n";
		exit(123);
    }
}
//----------------------------------------------------

// Calculates the density of states using one-particle greens function
void HamiltonianKH::Density_of_states(int N_e) {
	double omega = eigenvalues(0); // eigenvalues(0) - eigenvalues(N - 1);
	double domega = 0.001;

	vector<double> omega_vec;
	while (omega <= eigenvalues(N - 1)) {
		omega_vec.push_back(omega);
		omega += domega;
	}
	int nloop = omega_vec.size();

	ofstream DOSfile;
	stringstream Ustr, Nstr; 
	Ustr << setprecision(1) << fixed << U;
	Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
	DOSfile.open("DOS_n=" + Nstr.str() + "_U=" + Ustr.str() + ".txt");
	//DOSfile.open("DOS_2_U=" + Ustr.str() + ".txt");

	vector<double> resultSF(nloop);
		double DOS = 0;
		double maximum = 0;
		for (int w = 0; w < nloop; w++) {
			omega = omega_vec[w];
			DOS = 0;
			for (int n = 0; n < N; n++)
//				DOS += -1. / (double)L / pi * imag(1. / (omega + 2*domega*1i - eigenvalues(n)));
			if (DOS > maximum) 
				maximum = DOS;
			resultSF[w] = DOS;
		}
		for (int k = 0; k < nloop; k++)
			DOSfile << omega_vec[k] << "\t\t" << resultSF[k]/maximum + 5.1*U << endl;

	resultSF.~vector();
	omega_vec.~vector();
	DOSfile.close();
}

//----------------------------------------------------------------------------------------------
//----------------------------------------Rest methods & functions------------------------------
//----------------------------------------------------------------------------------------------

//Getting private fields for usage outside the class
mat HamiltonianKH::get_hamil() {
	return this->H;
}
//----------------------------------------------------

// Getting the energy eigenvalues from the private field
vec HamiltonianKH::get_energy() {
	return this->eigenvalues;
}
//----------------------------------------------------

//Factorial!!! - screw tgamma function, is shit
long long int factorial(int n) {
	if (n > 1) return n * factorial(n - 1);
	else return 1;
}

// Binomial: n po k = n!/( k!*(n-k)! )
long long int Binomial(int n, int k) {
	return factorial(n) / factorial(k) / factorial(n - k);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------


void Main_DOS_U(int L, int N_e, double t) {
	// Changing U
	for (double U = 0.0; U <= 2.8; U += 0.2) {
		double K, J_H;
		if (U != 0) {
			K = 4 * 0.15 * 0.15 / U;
			J_H = 0.25 * U;
		}
		else { K = 1; J_H = 0; }
		K = 0; J_H = 0;
		HamiltonianKH Object(L, N_e, t, U, K, J_H);
		Object.Hamiltonian();
		Object.Diagonalization();
		Object.Density_of_states(N_e);

		Object.~HamiltonianKH();
		cout << "U = " << U << " done!" << endl;
	}
	cout << "\n DOS for various U calculated" << endl;
}
