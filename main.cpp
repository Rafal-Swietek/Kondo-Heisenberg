#include "Hamiltonian.h"

using namespace std;

int main() {

    int L = 6; //chain length
	double t = 0; //electron hopping = t_00
    double U = 0; // electron repulsion
    double K = 1;// 0.0429; // spin exchange integral*.
    double J_H = 0;// U / 4.; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites

	// Main Program---------------------------------------
		    //Main_U(L, N_e, t);
            //Main_Jh(L, N_e, t, K, U);
    Main_Cv(L, N_e, t, K, U, J_H);
    for (int M = 50; M <= 300; M += 50) {
        Main_Cv_Lanczos(L, N_e, t, K, U, J_H, M);
        out << M << endl;
    }

    /*Lanczos Ham(L, N_e, t, U, K, J_H, -N_e + 6, 200);
    Ham.Lanczos_Diagonalization();
    out << Ham.eigenVal_L.t();*/
    //Main_DOS(L, N_e, t, K, U, J_H);
	//------------------------------------------------------------



	return 0;
}
