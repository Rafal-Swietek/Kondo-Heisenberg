#include "Hamiltonian.h"

using namespace std;

int main() {

    int L = 4; //chain length
	double t = 0.5; //electron hopping = t_00
    double U = 2.1; // electron repulsion
	double K = 0.0429; // spin exchange integral*.
    double J_H = U / 4.; // electron-spin interaction
    int N_e = 6; // numer of electrons - parity same as number of sites

	// Main Program---------------------------------------
		    //Main_U(L, N_e, t);
            //Main_Jh(L, N_e, t, K, U);
            HamiltonianKH Object(L, N_e, t, U, K, J_H);
            Object.Hamiltonian();
            Object.Diagonalization();
            srand(time(NULL));
            Object.Build_Lanczos_Hamil(Object.Create_Random_vec(), 100);
            out << "E   = " << Object.eigenvalues(0) << endl;
            out << "E_L = " << Object.eigenVal_L(0) << endl;
            Object.Heat_Capacity();
            Object.Heat_Capacity_Lanczos();

	//------------------------------------------------------------



	return 0;
}
