#ifndef ROUTINES
#define ROUTINES
#include "Lanczos.h"

/*
//----------------------------------------------------------------------------------------------
//---------------------------------------------DECLARATION---------------------------------------
//----------------------------------------------------------------------------------------------
double FermiLevel(int L, int N_e, double t, double K, double U, double J_H);

//Printing data
void printDOS(vec&& resultDOS, double U, double N_e, int L, std::vector<double>&& omega_vec, double maximum, double E_fermi);
void print_Cv(vec&& Cv, vec&& Cv_stand_dev, double U, double N_e, int L);
void print_chi(vec&& chi, vec&& chi_stand_dev, double U, double N_e, int L);
void print_Sq(vec&& Sq, vec&& Sq_stand_dev, double U, double N_e, int L, double T);


//Quantities averaged over spin blocks
void Heat_Capacity(std::vector<arma::vec>&& energies, vec& Cv);
void static_spin_susceptibility(std::vector<arma::vec>&& energies, vec& chi);


//Main routines
void Main_U(int L, int N_e, double t);
void Main_Jh(int L, int N_e, double t, double K, double U);
void Main_Cv(int L, int N_e, double t, double K, double U, double J_H);
void Main_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps);
void Main_DOS(int L, int N_e, double t, double K, double U, double J_H);
void Main_Sq(int L, int N_e, double t, double K, double U, double J_H);

void Cv_Umap(int L, int N_e, double t);
void DOS_Umap(int L, int N_e, double t);
void Lanczos_convergence(int L, int N_e);

*/

//----------------------------------------------------------------------------------------------
//---------------------------------------------DEFINITION---------------------------------------
//----------------------------------------------------------------------------------------------
namespace Routines{
    inline double FermiLevel(int L, int N_e, double t, double K, double U, double J_H) {
        double Ef;
        std::unique_ptr<HamiltonianKH> Object(new HamiltonianKH(L, N_e + 1, t, U, K, J_H, ((N_e + 1) % 2 == 0) ? 0 : 1));
        Object->Hamiltonian();
        Object->Diagonalization();
        Ef = Object->eigenvalues(0);

        Object.reset(new HamiltonianKH(L, N_e - 1, t, U, K, J_H, ((N_e - 1) % 2 == 0) ? 0 : 1));
        Object->Hamiltonian();
        Object->Diagonalization();
        Ef = (Ef + Object->eigenvalues(0)) / 2.0;

        return Ef;
    }

    inline void printDOS(vec&& resultDOS, double U, double N_e, int L, std::vector<double>&& omega_vec, double maximum, double E_fermi) {
        ofstream DOSfile;
        stringstream Ustr, Nstr;
        Ustr << setprecision(1) << fixed << U;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        DOSfile.open("DOS_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
        //DOSfile.open("DOS_dw=" + std::to_string(domega) + "_eta=" + std::to_string(eta) + ".txt");
        for (int k = 0; k < omega_vec.size(); k++)
            DOSfile << omega_vec[k] - E_fermi << "\t\t" << resultDOS[k] / maximum << endl;

        DOSfile.close();
    }
    inline void print_Cv(vec&& Cv, vec&& Cv_stand_dev, double U, double N_e, int L) {
        ofstream savefile;
        stringstream Ustr, Nstr;
        Ustr << setprecision(2) << fixed << U;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        savefile.open("Cv_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
        int k = 0;
        double T = dT;
        //out << "Cv:" << endl;
        while (T <= T_end) {
            savefile << T << "\t\t" << Cv(k) << "\t\t" << Cv_stand_dev(k) << endl; //save heat capacity to file
            //out << T << "\t\t" << Cv(k) << "\t\t" << Cv_stand_dev(k) << endl;
            T += dT; k++;
        }
        savefile.close();
    }
    inline void print_chi(vec&& chi, vec&& chi_stand_dev, double U, double N_e, int L) {
        ofstream savefile;
        stringstream Ustr, Nstr;
        Ustr << setprecision(2) << fixed << U;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        savefile.open("chi_0_L=" + std::to_string(L) + "_U=" + Ustr.str() + "_L.txt");
        int k = 0;
        double T = dT;
        //out << "chi:" << endl;
        while (T <= T_end) {
            savefile << T << "\t\t" << chi(k) << "\t\t" << chi_stand_dev(k) << endl;
            //out << T << "\t\t" << chi(k) << "\t\t" << chi_stand_dev(k) << endl;
            T += dT; k++;
        }
        savefile.close();
    }
    inline void print_Sq(vec&& Sq, vec&& Sq_stand_dev, double U, double N_e, int L, double T) {
        ofstream savefile;
        stringstream Ustr, Nstr;
        Ustr << setprecision(2) << fixed << U;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        savefile.open("Sq_L=" + std::to_string(L) + "_U=" + Ustr.str() + "_T=" + std::to_string(T) + ".txt");
        savefile << std::setprecision(8) << fixed;
        out << "Sq : " << endl;
        for (int l = 0; l <= L; l++) {
            double q = (double)l * pi / ((double)L + 1.0);
            savefile << q << "\t\t" << Sq(l) << "\t\t" << Sq_stand_dev(l) << endl;
            out << q / pi << "\t" << Sq(l) << "\t\t" << Sq_stand_dev(l) << endl;
        }
        savefile.close();
    }

    // ED quantities
    inline void Heat_Capacity(std::vector<arma::vec>&& energies, vec& Cv) {
        auto temperature = prepare_parameterVec(dT, T_end, dT);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double T = temperature[k];
            double heat_capacity = 0;
            int l = 0;
            for (auto E : energies) {
                double Sz = -((double)E.size() - 1.0) / 2. + (double)l;
                double E_av = 0;
                double E2_av = 0;
                double Z = 0;
                //out << Sz << endl;
                for (int n = 0; n < E.size(); n++) {
                    E_av += E(n) * std::exp(-E(n) / T);
                    E2_av += E(n) * E(n) * std::exp(-E(n) / T);
                    Z += std::exp(-E(n) / T);
                }
                E_av = E_av / Z;
                E2_av = E2_av / Z;
                heat_capacity += (E2_av - E_av * E_av) / T / T;
                l++;
            }
            int L = (energies.size() - 1) * 2 / 3; // n=1.5 filling
            Cv(k) = heat_capacity / (double)L / (double)energies.size();

            //out << T << " " << Cv(k) << " " << L << endl;
        }
    }
    inline void static_spin_susceptibility(std::vector<arma::vec>&& energies, vec& chi) {
        auto temperature = prepare_parameterVec(dT, T_end, dT);
        int L = (energies.size() - 1) * 2. / 3.; // n = 1.5 filling
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double T = temperature[k];
            double Z = 0;
            double X_0 = 0;
            int l = 0;
            for (auto E : energies) {
                double Sz = -((double)energies.size() - 1.0) / 2. + (double)l; // --gKH
                //double Sz = -(double)L / 2. + (double)l; // --heisenberg
                for (int n = 0; n < E.size(); n++) {
                    Z += std::exp(-E(n) / T);
                    X_0 += Sz * Sz * std::exp(-E(n) / T);
                }
                l++;
            }
            chi(k) = X_0 / Z / T / (double)L;
        }
    }
    //--

    inline void Cv_Umap(int L, int N_e, double t) {
        ofstream savefile;
        stringstream Nstr;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        savefile.open("C_V_L=" + std::to_string(L) + "_Umap.txt");
        for (int k = 2; k < 40; k += 1) {
            double U = (double)k / 10.0;
            double K, J_H;
            K = 4 * 0.15 * 0.15 / U;
            J_H = 0.25 * U;
            std::vector<arma::vec> energies;
            //std::unique_ptr<std::vector<arma::vec>> energies(new std::vector<arma::vec>);
            for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
                //HamiltonianKH Hamil(L, N_e, t, U, K, J_H, Sz);
                std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, Sz));
                Hamil->Hamiltonian();
                Hamil->Diagonalization();
                energies.emplace_back(Hamil->eigenvalues);
            }
            vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
            Heat_Capacity(std::move(energies), Cv);
            int q = 0;
            double T = dT;
            while (T <= T_end) {
                savefile << T << "\t\t" << U << "\t\t" << Cv(q) << endl; //save heat capacity to file
                T += dT; q++;
            }
            //Main_DOS(L, N_e, t, K, U, J_H);
            //Main_Cv(L, N_e, t, K, U, J_H);

            out << "U = " << U << " done!" << endl;
        }
        savefile.close();
    }
    inline void DOS_Umap(int L, int N_e, double t) {
        ofstream savefile;
        stringstream Nstr;
        Nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        savefile.open("DOS_L=" + std::to_string(L) + "_Umap.txt");
        for (int k = 1; k < 80; k += 1) {
            double U = (double)k / 10.0;
            double K, J_H;
            K = 4 * 0.15 * 0.15 / U;
            J_H = 0.25 * U;

            std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1));
            Hamil->Hamiltonian();
            Hamil->Diagonalization();
            vector<double> omega_vec = prepare_parameterVec(Hamil->eigenvalues(0) - 100 * domega, Hamil->eigenvalues(Hamil->eigenvalues.size() - 1) + 100 * domega, domega);
            vec DOS = Hamil->Total_Density_of_states(std::move(omega_vec));
            int q = 0;
            double Ef = FermiLevel(L, N_e, t, K, U, J_H);
            for (int q = 0; q < omega_vec.size(); q++) {
                savefile << omega_vec[q] - Ef << "\t\t" << U << "\t\t" << DOS(q) << endl; //save heat capacity to file
            }

            out << "U = " << U << " done!" << endl;
        }
        savefile.close();
    }

    inline void Main_Cv(int L, int N_e, double t, double K, double U, double J_H) {
        std::vector<arma::vec> energies;
        vec Sq(L + 1, fill::zeros);
        vec Sq_stand_dev(L + 1, fill::zeros);
        vec chi0_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        double Z = 0;
        std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, -N_e));
        for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
            Hamil->update_parameters(t, U, K, J_H, Sz);
            Hamil->Diagonalization();
            Sq += Hamil->static_structure_factor(T);
            Z += Hamil->partition_function(T);
            energies.emplace_back(std::move(Hamil->eigenvalues));
            //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
        }
        Sq = Sq / Z;
        print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, N_e, L, T);

        vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        Heat_Capacity(std::move(energies), Cv);
        print_Cv(std::move(Cv), std::move(Cv_stand_dev), U, N_e, L);

        vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        static_spin_susceptibility(std::move(energies), chi_0);
        print_chi(std::move(chi_0), std::move(chi0_stand_dev), U, N_e, L);
    }
    inline void Main_Sq(int L, int N_e, double t, double K, double U, double J_H) {
        vec Sq(L + 1, fill::zeros);
        vec Sq_stand_dev(L + 1, fill::zeros);
        double Z = 0;
        double T = 1.0;
        for (double T = 0.01; T <= 1; T += 0.02) {
            for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
                std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, Sz));
                Hamil->Diagonalization();
                Sq += Hamil->static_structure_factor(T);
                Z += Hamil->partition_function(T);
                //out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
            }
            Sq = Sq / Z;
            print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, N_e, L, T);
        }
    }
    inline void Main_DOS(int L, int N_e, double t, double K, double U, double J_H) {
        //HamiltonianKH Hamil(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1);
        std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1));
        Hamil->Diagonalization();
        vector<double> omega_vec = prepare_parameterVec(Hamil->eigenvalues(0) - 100 * domega, Hamil->eigenvalues(Hamil->eigenvalues.size() - 1) + 100 * domega, domega);
        vec DOS = Hamil->Total_Density_of_states(std::move(omega_vec)); //rvalue of DOS
        double maximum = max(DOS);
        double Ef = FermiLevel(L, N_e, t, K, U, J_H);
        //Hamil->printEnergy(Ef);
        printDOS(std::move(DOS), U, N_e, L, std::move(omega_vec), maximum, Ef);
    }
    inline void Main_Lanczos(int L, int N_e, double t, double K, double U, double J_H, int M, int random_steps) {
        auto temperature = prepare_parameterVec(dT, T_end, dT);
        vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_R = Cv;
        vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec chi_0_2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec chi_R = chi_0;
        vec Sq(L + 1, fill::zeros);
        vec Sq2(L + 1, fill::zeros);
        vec Sq_R = Sq;
        vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

        vec Sq_stand_dev(L + 1, fill::zeros);
        vec chi0_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        double Z_constT = 0;
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, M));
        Hamil->Diagonalization();
        E0 = Hamil->eigenvalues(0); // set global GS energy
        out << E0 << endl;
        int R = (calculate_stand_dev) ? 50 : 1; // for bigger systems lesser averaging to decrease computation time for now  
        // ((L >= 8) ? 20 : 50)
        for (int r = 0; r < R; r++) {
            Cv_R.zeros(); Sq_R.zeros();  chi_R.zeros(); Z_constT = 0; Z.zeros();
            for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
                Hamil->update_parameters(t, U, K, J_H, Sz);
                if (M > Hamil->N)  Hamil->lanczos_steps = Hamil->N;
                else Hamil->lanczos_steps = M;

                Cv_R += Hamil->Heat_Capacity_Lanczos(random_steps) / double(N_e + 1.0);
                chi_R += Hamil->static_spin_susceptibility(random_steps);
                Z += Hamil->partition_function;
                if (calculate_Sq)
                    Sq_R += Hamil->Sq_lanczos(random_steps, T);
                Z_constT += Hamil->Z_constT;
                out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
            }
            Sq_R = Sq_R / Z_constT; chi_R = chi_R / Z;
            Sq += Sq_R / (double)R; Sq2 += arma::square(Sq_R) / (double)R;
            chi_0 += chi_R / (double)R; chi_0_2 += arma::square(chi_R) / (double)R;
            Cv += Cv_R / (double)R; Cv_2 += arma::square(Cv_R) / (double)R;
            out << "r = " << r << endl;
            if (calculate_stand_dev && r % 1 == 0) {
                std::ofstream fileR("Cv_L=" + std::to_string(L) + "_R=" + std::to_string(r) + ".txt");
                //fileR << (Cv * (double)R / (double)r);
                fileR << Cv_R;
                fileR.close();
                fileR.open("Sq_L=" + std::to_string(L) + "_R=" + std::to_string(r) + ".txt");
                fileR << Sq_R;
                fileR.close();
                fileR.open("X_L=" + std::to_string(L) + "_R=" + std::to_string(r) + ".txt");
                fileR << chi_R;
                fileR.close();
            }
        }
        Sq_stand_dev = arma::sqrt(Sq2 - square(Sq));
        print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, N_e, L, T);

        Cv_stand_dev = arma::sqrt(Cv_2 - square(Cv));
        print_Cv(std::move(Cv), std::move(Cv_stand_dev), U, N_e, L);

        chi0_stand_dev = arma::sqrt(chi_0_2 - square(chi_0));
        print_chi(std::move(chi_0), std::move(chi0_stand_dev), U, N_e, L);
        E0 = 0;
    }

    inline void Main_Jh(int L, int N_e, double t, double K, double U) {
        double J_H = 0.0 * U;
        while (J_H <= 0.4 * U) {
            //Main_DOS(L, N_e, t, K, U, J_H);
            //Main_Cv(L, N_e, t, K, U, J_H);
            Main_Lanczos(L, N_e, t, K, U, J_H, (L == 4) ? 100 : 150, 20);
            out << "J_H/U = " << J_H / U << " done!" << endl;
            J_H += 0.05 * U;
        }
    }
    inline void Main_U(int L, int N_e, double t) {
        double U = 0.1 * W;
        while (U <= 3*W) {
            double K, J_H;
            K = 4 * 0.15 * 0.15 / U;
            J_H = 0.25 * U;
            Main_Lanczos(L, N_e, t, K, U, J_H, (L == 4) ? 100 : 150, 15);

            out << "U = " << U << " done!" << endl;
            U += 0.1 * W;
        }
    }

    inline void Lanczos_convergence(int L, int N_e) {
        double t = 0.5, U = 2.1, K = 4 * 0.15 * 0.15 / U, J_H = 0.25 * U;
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, 500));
        vec random_vec = Create_Random_vec(Hamil->N);
        Hamil->Lanczos_convergence(random_vec);
    }

}
#endif
