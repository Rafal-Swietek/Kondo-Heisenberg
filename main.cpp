#include "Hamiltonian.h"

using namespace std;

int main() {
    srand(time(NULL));

    clock_t tStart = clock();
    int L = 8; //chain length
    double t = 0.5; //electron hopping = t_00
    double U = 2.1; // electron repulsion
    double K = 4 * 0.15 * 0.15 / U; // spin exchange integral*.
    double J_H = U / 4.; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites
    int M = (L == 4)? 100 : 150; // number of Lanczos steps
    int R = 15; // number of randomized steps for thermal average
    // Main Program---------------------------------------


    std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
    Hamil->Lanczos_GroundState();
    Hamil->show_ground_state();
    //out << Hamil->correlation_matrix();

    //out << endl << " Lanczos:\n q \t\t S(q) " << endl;
    //Main_Cv_Lanczos(L, N_e, t, K, U, J_H, M, R);
    
    /*std::unique_ptr<HamiltonianKH> Hamil2(new HamiltonianKH(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1));
    Hamil2->Diagonalization();
    out << Hamil2->correlation_matrix();*/
    //Main_U(L, N_e, t);
    //out << endl << " ED:\n q \t\t S(q) " << endl;
    //Main_Cv(L, N_e, t, K, U, J_H);
    //Lanczos_convergence(L, N_e);
	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}