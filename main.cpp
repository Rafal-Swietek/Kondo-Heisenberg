#include "Hamiltonian.h"

using namespace std;

int main() {
    clock_t tStart = clock();
    int L = 4; //chain length
    double t = 0.5; //electron hopping = t_00
    double U = 2.1; // electron repulsion
    double K = 0.0492; // spin exchange integral*.
    double J_H = U / 4.; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites
    int M = 150; // number of Lanczos steps
    int R = 10; // number of randomized steps for thermal average
    // Main Program---------------------------------------

    Main_Cv(L, N_e, 0, 1, 0, 0);
    /*std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
    Hamil->Lanczos_GroundState();
    Hamil->show_ground_state();*/
    /*Main_Cv_Lanczos(4, 6, t, K, U, J_H, M, R);
    Main_Cv_Lanczos(6, 9, t, K, U, J_H, M, R);
    Main_Cv_Lanczos(8, 12, t, K, U, J_H, M, R);*/
	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}