#include "subroutines.h"
using namespace Routines;

//externs
double domega = 0.001; // freq step for DOS
double eta = 0.02; // ~(peak width) in DOS, eta->0
double T = 0.01; // temperature at which Sq is calculated
double dT = 0.0001; // temperature step
double T_end = 10.0; // end temperature, saturated entropy?
bool PBC = 0; // allow periodic boundary conditions?
int model = 0;
// model = 0 - gKH
// model = 1 - Heisenberg
// model = 2 - Hubbard


int main() {
    srand(time(NULL));

    clock_t tStart = clock();
    int L = 8; //chain length
    double t = 0.5; //electron hopping = t_00
    double U = W; // electron repulsion
    double K = 4 * 0.15 * 0.15 / U; // spin exchange integral*.
    double J_H = U / 4.; // electron-spin interaction
    int N_e = 3 * L / 2; // numer of electrons - parity same as number of sites
    int M = 300; // number of Lanczos steps
    int R = 50; // number of randomized steps for thermal average

    out << "t=" << t << "\tU/W=" << U / W << "\tK=" << K << "\tJh/U=" << J_H / U << endl;
    // Main Program---------------------------------------
    /*arma::arma_version ver;
    std::cout << "ARMA version: " << ver.as_string() << std::endl << std::endl;*/

    //Main_Cv(L, N_e, t, K, U, J_H);
    //Main_X(L, N_e, t, K, W);

    //Main_Lanczos(L, N_e, t, K, U, J_H, M, R);

    std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, 0, M));
    Hamil->SSF_T0();
    //Find_Cv_peaks(L);
    /*
    model = 2;
    for (U = 0; U <= W; U += 0.01*W) {
        //Main_DOS(L, N_e, t, (U == 0) ? 0 : 4 * 0.15 * 0.15 / U, U, U / 4.);
        Main_DOS(L, L, 0.2, 0, U, 0);
    }*/
    //Main_U(L, L, 1, R);
    //Energu_U_sweep(L, N_e, t);
    //Main_Jh(L, N_e, t, K, U, R);
    //Main_X(L, N_e, t, K, U, R);
    //Main_U(L, N_e, t, R, M);
    //Sq_max_map(L, M, R, t);
  
    //Lanczos_convergence(L, N_e);
	//------------------------------------------------------------
    out << "Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}
// TODO: fsosffodmfkmsfsfsif
// ! this is important
// ? question
//// WRONG

/*ofstream files1("S_M=" + std::to_string(M) + ".txt");
    ofstream files2("Cv_M=" + std::to_string(M) + ".txt");
    std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
    files1 << Hamil->entropy(R);
    files2 << Hamil->Heat_Capacity_Lanczos(R);
    files1.close();
    files2.close();*/

    /*for (M = 500; M >= 50; M -= 50) {
        ofstream files1("S_M=" + std::to_string(M) + ".txt");
        ofstream files2("Cv_M=" + std::to_string(M) + ".txt");
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
        files1 << Hamil->entropy(R);
        files2 << Hamil->Heat_Capacity_Lanczos(R);
        files1.close();
        files2.close();
    }*/
/*
for (U = W; U <= 3 * W; U += W) {
    for (int L = 2; L <= 14; L += 2) {
        N_e = 3 * L / 2;
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
        vec Sq(L + 1, fill::zeros);
        vec zereo(L + 1, fill::zeros);
        Sq = Hamil->Sq_T0(1);
        print_Sq(std::move(Sq), std::move(zereo), U, J_H, N_e, L, T);
        Hamil->show_ground_state();
    }
}*/


/*
// n = 2;
for (U = W; U <= 3 * W; U += W) {
    K = 4 * 0.15 * 0.15 / U; J_H = U / 4.;
    double n = 0;
    for (int i = 0; i <= L; i++) {
        n = 2 + (double)i / (double)L;
        N_e = L + i;
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
        vec Sq(L + 1, fill::zeros);
        vec zereo(L + 1, fill::zeros);
        Sq = Hamil->Sq_T0(1);
        print_Sq(std::move(Sq), std::move(zereo), U, J_H, N_e, L, T);
        Hamil->show_ground_state();
        out << "Filling n= " << n << " at U/W= " << U / W << " done" << endl;
    }
}*/

/*U = 0.1 * W;
   out << "U/W" << "\t\t\t\t T*" << endl;
   while (U <= 3.05 * W) {
       double K, J_H;
       K = 4 * 0.15 * 0.15 / U;
       J_H = 0.25 * U;
       std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
       Hamil->Diagonalization();
       out << U / W << "\t\t" << (Hamil->eigenvalues(M - 1) - Hamil->eigenvalues(0)) / (double)M << endl;
       U += 0.1 * W;
   }*/