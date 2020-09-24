#include "Hamiltonian.h"

using namespace std;

int main() {

	int L = 4; //chain length
	double t = 1; //electron hopping = t_00
	double U = 1; // electron repulsion
	double K = 1; // spin exchange integral*.
	double J_H = 1; // electron-spin interaction
	int N_e = 4; // numer of electrons - parity same as number of sites

	// Main Program---------------------------------------
		//Main_DOS_U(L, N_e, t);
        cout << "Ground state for:" << endl;
			HamiltonianKH Object(L, N_e, t, U, K, J_H);
			Object.Hamiltonian();
			Object.Diagonalization();
			//Object.Density_of_states(N_e);
            cout << "E = " << setprecision(8) << Object.get_energy()(0) << endl;

	//------------------------------------------------------------



	return 0;
}
