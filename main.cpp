#include "Hamiltonian.h"

using namespace std;

int main() {
    clock_t tStart = clock();
    int L = 4; //chain length
    double t = 0; //electron hopping = t_00
    double U = 0; // electron repulsion
    double K = 1; // spin exchange integral*.
    double J_H = 0; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites
    // Main Program---------------------------------------

    //Main_U(L, N_e, t);
    Main_Cv(L, N_e, t, K, U, J_H);
    //Main_Cv_Lanczos(L, N_e, t, K, U, J_H, 150, 10);
    //Main_Cv(L, N_e, 0.5, 0.0429, 2.1, 2.1 / 4);

    //out << sizeof(intmax_t);
   /* Lanczos Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, 150);
    Hamil.Lanczos_GroundState();
    Hamil.show_ground_state();*/

    //Main_Cv_Lanczos(L, N_e, t, K, U, J_H, 150, 100);

	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}