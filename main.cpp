#include "Hamiltonian.h"

using namespace std;

int main() {

	/*time_t t_start, t_end;
	time(&t_start);*/

	int L = 4; //chain length
	double t = 1; //electron hopping = t_00
	double U = 1; // electron repulsion
	double K = 1; // spin exchange integral*.
	double J_H = 1; // electron-spin interaction
	int N_e = 4; // numer of electrons - parity same as number of sites

	/*TCHAR Npath[MAX_PATH];
	GetCurrentDirectory(MAX_PATH, Npath);
	CreateDirectory(L"./Results/", NULL);
	SetCurrentDirectory(L"./Results/");*/

	// Main Program---------------------------------------
		//Main_DOS_U(L, N_e, t);
		cout << "Ground state for:" << endl;
		/*for (int N_e = L; N_e <= 2*L; N_e += 2) {
			HamiltonianKH Object(L, N_e, t, U, K, J_H);
			Object.Hamiltonian();
			Object.Diagonalization();
			cout << "n = N/L = " << (double)N_e / (double)L << "\t->\tE = " << setprecision(8) << Object.get_energy()(0) << endl;
			Object.static_StructureFactor();
			Object.~HamiltonianKH();
		}*/
			HamiltonianKH Object(L, N_e, t, U, K, J_H);
			Object.Hamiltonian();
			Object.Diagonalization();
			//Object.Density_of_states(N_e);
			cout << "E = " << setprecision(8) << Object.get_energy()(0) << endl;
            //Object.~HamiltonianKH();

	//------------------------------------------------------------


	/*SetCurrentDirectory(Npath);
	//Calculate execution time
	time(&t_end);
	int hours, minutes, seconds;
	hours = static_cast<int>((t_end - t_start) / 3600);
	minutes = static_cast<int>((t_end - t_start - 3600 * hours) / 60);
	seconds = t_end - t_start - 3600 * hours - 60 * minutes;
	cout << "\nProgram executed in:\t" << fixed << hours << "h " << minutes << "min " << seconds << "sec" << setprecision(2) << endl;*/

	return 0;
}
