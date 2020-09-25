#include "Hamiltonian.h"

using namespace std;

int main() {

	int L = 6; //chain length
	double t = 1; //electron hopping = t_00
    double U = 1; // electron repulsion
	double K = 1; // spin exchange integral*.
	double J_H = 1; // electron-spin interaction
	int N_e = 6; // numer of electrons - parity same as number of sites

	// Main Program---------------------------------------
		//Main_DOS_U(L, N_e, t);
        	cout << "Ground state for:" << endl;
		HamiltonianKH Object(L, N_e, t, U, K, J_H);
		Object.Hamiltonian();
		Object.Diagonalization();
		//Object.Density_of_states(N_e);
            	out << "E = " << Object.eigenvalues(0) << endl;

	//------------------------------------------------------------



	return 0;
}
