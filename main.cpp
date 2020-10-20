#include "Hamiltonian.h"

using namespace std;

int main() {
    clock_t tStart = clock();
    int L = 4; //chain length
	double t = 0; //electron hopping = t_00
    double U = 0; // electron repulsion
    double K = 1.0; // spin exchange integral*.
    double J_H = 0; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites
	// Main Program---------------------------------------

    //Main_U(L, N_e, t);

    Main_Cv(L, N_e, t, K, U, J_H);

	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}