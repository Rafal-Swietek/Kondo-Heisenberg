#include "Hamiltonian.h"

using namespace std;

int main() {
    clock_t tStart = clock();
    int L = 8; //chain length
	double t = 0.5; //electron hopping = t_00
    double U = 2.1; // electron repulsion
    double K = 0.0429; // spin exchange integral*.
    double J_H = U / 4.; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites

	// Main Program---------------------------------------
    /*out << "L=6:" << endl;
    Main_U(L, N_e, t);
    Main_U(L, N_e - 1, t);
    Main_U(L, N_e + 1, t);*/
    /*HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
    Hamil.Hamiltonian();
    Hamil.Diagonalization();
    Hamil.show_ground_state();*/
    int M = 200;
    Lanczos Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M);
    Hamil.Lanczos_Diagonalization();
    out << Hamil.eigenvalues(0) << endl;
    Hamil.Lanczos_GroundState();
    Hamil.show_ground_state();
    

    /*out << "Cv Lanczos for different M & R" << endl;
#pragma omp parallel for collapse(2)
    for (int R = 10; R <= 20; R += 5) {
        for (int M = 50; M <= 350; M += 100) {
            Main_Cv_Lanczos(L, N_e, t, K, U, J_H, M, R);
        }
    }*/
	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}