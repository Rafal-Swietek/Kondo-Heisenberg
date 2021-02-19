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
        DOSfile << std::setprecision(16) << std::fixed;
        for (int k = 0; k < omega_vec.size(); k++)
            DOSfile << omega_vec[k] - E_fermi << "\t\t" << resultDOS[k] / maximum << endl;

        DOSfile.close();
    }
    inline void print_Cv(vec&& Cv, vec&& Cv_stand_dev, double U, double J_H, double N_e, int L) {
        stringstream Ustr, Jhstr, nstr;
        nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        Ustr << setprecision(2) << fixed << U / W;
        Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
        std::ofstream savefile("Cv_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
        int k = 0;
        double T = dT;
        //out << "Cv:" << endl;
        savefile << std::setprecision(16) << std::fixed;
        while (T <= T_end) {
            savefile << T << "\t\t" << Cv(k) << "\t\t" << Cv_stand_dev(k) << endl; //save heat capacity to file
            //out << T << "\t\t" << Cv(k) << "\t\t" << Cv_stand_dev(k) << endl;
            T += dT; k++;
        }
        savefile.close();
    }
    inline void print_chi(vec&& chi, vec&& chi_stand_dev, double U, double J_H, double N_e, int L) {
        stringstream Ustr, Jhstr, nstr;
        nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        Ustr << setprecision(2) << fixed << U / W;
        Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
        std::ofstream savefile("X0_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
        int k = 0;
        double T = dT;
        savefile << std::setprecision(16) << std::fixed;
        //out << "chi:" << endl;
        while (T <= T_end) {
            savefile << T << "\t\t" << chi(k) << "\t\t" << chi_stand_dev(k) << endl;
            //out << T << "\t\t" << chi(k) << "\t\t" << chi_stand_dev(k) << endl;
            T += dT; k++;
        }
        savefile.close();
    }
    inline void print_Sq(vec&& Sq, vec&& Sq_stand_dev, double U, double J_H,  double N_e, int L, double T) {
        stringstream Ustr, Jhstr, nstr;
        nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        Ustr << setprecision(2) << fixed << U / W;
        Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
        std::ofstream savefile("Sq_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str()\
            + "U__" + "_T=" + std::to_string(T) + "_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
        savefile << std::setprecision(16) << std::fixed;
        out << "Sq : " << endl;
        for (int l = 0; l <= L; l++) {
            double q = (double)l * pi / ((double)L + 1.0);
            savefile << q << "\t\t" << Sq(l) << "\t\t" << Sq_stand_dev(l) << endl;
            out << q / pi << "\t" << Sq(l) << "\t\t" << Sq_stand_dev(l) << endl;
        }
        savefile.close();
    }
    inline void print_S(vec&& S, vec&& S_stand_dev, double U, double J_H, double N_e, int L) {
        stringstream Ustr, Jhstr, nstr;
        nstr << setprecision(2) << fixed << (double)N_e / (double)L;
        Ustr << setprecision(2) << fixed << U / W;
        Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
        std::ofstream savefile("S_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
        int k = 0;
        double T = dT;
        savefile << std::setprecision(16) << std::fixed;
        //out << "Cv:" << endl;
        while (T <= T_end) {
            savefile << T << "\t\t" << S(k) << "\t\t" << S_stand_dev(k) << endl; //save entropy to file
            //out << T << "\t\t" << Cv(k) << "\t\t" << Cv_stand_dev(k) << endl;
            T += dT; k++;
        }
        savefile.close();
    }
    // ED quantities
    inline void Heat_Capacity(std::vector<arma::vec>& energies, vec& Cv) {
        auto temperature = prepare_parameterVec(dT, T_end, dT);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double T = temperature[k];
            double heat_capacity = 0;
            int l = 0;
            //out << "size: " << energies.size() << endl << endl;
            for (auto E : energies) {
                double Sz = -((double)energies.size() - 1.0) / 2. + (double)l;
                double E_av = 0;
                double E2_av = 0;
                double Z = 0;
                //out << Sz << endl;
                for (int n = 0; n < E.size(); n++) {
                    E_av += E(n) * std::exp(-(E(n) - E0) / T);
                    E2_av += E(n) * E(n) * std::exp(-(E(n) - E0) / T);
                    Z += std::exp(-(E(n) - E0) / T);
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
    inline void static_spin_susceptibility(std::vector<arma::vec>& energies, vec& chi) {
        auto temperature = prepare_parameterVec(dT, T_end, dT);
        int L = (energies.size() - 1) * 2. / 3.; // n = 1.5 filling gKH
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double T = temperature[k];
            double Z = 0;
            double X_0 = 0;
            int l = 0;
            for (auto E : energies) {
                double Sz = -((double)energies.size() - 1.0) / 2. + (double)l; // --gKH
                //out << Sz << endl;
                //double Sz = -(double)L / 2. + (double)l; // --heisenberg
                for (int n = 0; n < E.size(); n++) {
                    Z += std::exp(-(E(n) - E0) / T);
                    X_0 += Sz * Sz * std::exp(-(E(n) - E0) / T);
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
            Heat_Capacity(energies, Cv);
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

    inline void GS_map_to_Heis(int L, int N_e, double t, double U, double K, double J_H, int M, int R) {
        int Sz;
        switch (model) {
        case 0: //gKH
            Sz = (N_e % 2 == 0) ? 0 : 1;
            break;
        case 1: //Heisenberg
            Sz = (L % 2 == 0) ? 0 : 1;
            t = 0; U = 0; J_H = 0; N_e = 0;
            break;
        case 2: // Hubbard
            Sz = (N_e % 2 == 0) ? 0 : 1;
            K = 0; J_H = 0;
            break;
        }
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, Sz, M));
        Hamil->Lanczos_GroundState();
        Hamil->show_ground_state();
        //out << Hamil->correlation_matrix() << endl;
        //out << Hamil->Sq_T0();
    }
    inline void Main_Cv(int L, int N_e, double t, double K, double U, double J_H) {
        std::vector<arma::vec> energies;
        vec Sq(L + 1, fill::zeros);
        vec Sq_stand_dev(L + 1, fill::zeros);
        vec chi0_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_stand_dev(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        double Z = 0;
        std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (L % 4 == 0) ? 0 : 1));
        Hamil->Diagonalization();
        E0 = Hamil->eigenvalues(0); // set global GS energy
        out << "Ground state E0 = " << E0 << endl;
        for (int Sz = -N_e; Sz <= N_e; Sz += 2) {
            Hamil->update_parameters(t, U, K, J_H, Sz);
            Hamil->Diagonalization();
            if (calculate_Sq) {
                Sq += Hamil->static_structure_factor(T);
                Z += Hamil->partition_function(T);
            }
            energies.emplace_back(std::move(Hamil->eigenvalues));
            out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
        }
        if (calculate_Sq) {
            Sq = Sq / Z;
            print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, J_H, N_e, L, T);
        }

        vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        Heat_Capacity(energies, Cv);
        print_Cv(std::move(Cv), std::move(Cv_stand_dev), U, J_H, N_e, L);

        vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        static_spin_susceptibility(energies, chi_0);
        print_chi(std::move(chi_0), std::move(chi0_stand_dev), U, J_H, N_e, L);
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
            print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, J_H, N_e, L, T);
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
        int Sz_tot;
        switch (model) {
        case 0: //gKH
            Sz_tot = (N_e % 2 == 0) ? 0 : 1;
            break;
        case 1: //Heisenberg
            Sz_tot = (L % 2 == 0) ? 0 : 1;
            t = 0; U = 0; J_H = 0; N_e = 0;
            break;
        case 2: // Hubbard
            Sz_tot = (N_e % 2 == 0) ? 0 : 1;
            K = 0; J_H = 0;
            break;
        }
        //GS_map_to_Heis(L, N_e, t, U, K, J_H, M, random_steps);
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, Sz_tot, M));
        Hamil->Diagonalization();
        E0 = Hamil->eigenvalues(0); // set global GS energy
        out << "Ground state E0 = " << E0 << endl;

        vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Cv_2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec chi_0_2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec Sq(L + 1, fill::zeros);
        vec Sq2(L + 1, fill::zeros);
        vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec S(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
        vec S2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

        double Z_constT = 0;
        int R = (calculate_stand_dev) ? 20 : 1; // for bigger systems lesser averaging to decrease computation time for now  
        // ((L >= 8) ? 20 : 50)
        int Sz_max = (model == 0) ? N_e : ((model == 1) ? L : 2*L-N_e);
        out << "dim(H) = " << Hamil->N << endl;
        for (int r = 0; r < R; r++) {
            Z_constT = 0; Z.zeros();
            vec Cv_R, Sq_R, chi_R, S_R;
            if (calulate_X0) {
                chi_R = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
                for (int Sz = -Sz_max; Sz <= Sz_max; Sz += 2) {
                    std::unique_ptr<Lanczos> Hamil2(new Lanczos(L, N_e, t, U, K, J_H, Sz, M));
                    chi_R += Hamil2->static_spin_susceptibility(random_steps);// / (Sz_max + 1);
                    Z += Hamil2->partition_function;
                    out << "Sector Sz = " << double(Sz) / 2.0 << "done" << endl;
                }
                chi_R = chi_R / Z;
                chi_0 += chi_R / (double)R; 
                chi_0_2 += arma::square(chi_R) / (double)R;
            }
            if (calculate_Sq) {
                Sq_R;
                if (T == 0) Sq_R = Hamil->Sq_T0(random_steps);
                else Sq_R = Hamil->Sq_lanczos(random_steps, T) / Hamil->Z_constT;
                Sq += Sq_R / (double)R;
                Sq2 += arma::square(Sq_R) / (double)R;
            }
            if (calculate_Cv) {
                Cv_R = Hamil->Heat_Capacity_Lanczos(random_steps);// / double(N_e + 1.0);
                Cv += Cv_R / (double)R;
                Cv_2 += arma::square(Cv_R) / (double)R;
            }
            if (calculate_entropy) {
                S_R = Hamil->entropy(random_steps);
                S += S_R / (double)R;
                S2 += arma::square(S_R) / (double)R;
            }
            if (calculate_stand_dev && r % 1 == 0) {
                stringstream Ustr, Jhstr, nstr;
                nstr << setprecision(2) << fixed << (double)N_e / (double)L;
                Ustr << setprecision(2) << fixed << U / W;
                Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
                std::ofstream fileR;
                fileR << std::setprecision(16) << std::fixed;
                if (calculate_Cv) {
                    fileR.open("Cv_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + "__R=" + std::to_string(r) + ".txt");
                    //fileR << (Cv * (double)R / (double)r);
                    fileR << Cv_R;
                    fileR.close();
                }
                if (calculate_entropy) {
                    fileR.open("S_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + "__R=" + std::to_string(r) + ".txt");
                    fileR << S_R;
                    fileR.close();
                }
                if (calculate_Sq) {
                    fileR.open("Sq_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + "_T=" + std::to_string(T) + "__R=" + std::to_string(r) + ".txt");
                    fileR << Sq_R;
                    fileR.close();
                }
                if (calulate_X0) {
                    fileR.open("X0_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + "__R=" + std::to_string(r) + ".txt");
                    fileR << chi_R;
                    fileR.close();
                }
            }
            out << "r = " << r << " fine" << endl;
        }
        if (calculate_Sq) {
            vec Sq_stand_dev = arma::sqrt(Sq2 - square(Sq));
            print_Sq(std::move(Sq), std::move(Sq_stand_dev), U, J_H, N_e, L, T);
        }
        if (calulate_X0) {
            vec chi0_stand_dev = arma::sqrt(chi_0_2 - square(chi_0));
            print_chi(std::move(chi_0), std::move(chi0_stand_dev), U, J_H, N_e, L);
        }

        if (calculate_Cv) {
            vec Cv_stand_dev = arma::sqrt(Cv_2 - square(Cv));
            print_Cv(std::move(Cv), std::move(Cv_stand_dev), U, J_H, N_e, L);
        }

        if (calculate_entropy) {
            vec S_stand_dev = arma::sqrt(S2 - square(S));
            print_Cv(std::move(S), std::move(S_stand_dev), U, J_H, N_e, L);
        }
        E0 = 0;
    }

    inline void Main_Jh(int L, int N_e, double t, double K, double U, int R) {
        double J_H = 0.0 * U;
        while (J_H <= 0.4 * U) {
            //Main_DOS(L, N_e, t, K, U, J_H);
            //Main_Cv(L, N_e, t, K, U, J_H);
            Main_Lanczos(L, N_e, t, K, U, J_H, 150, R);
            out << "J_H/U = " << J_H / U << " done!" << endl;
            J_H += (J_H / U < 0.1) ? 0.01 * U : 0.05 * U;
        }
    }
    inline void Main_U(int L, int N_e, double t, int R) {
        double U = 0.1 * W;
        while (U <= 3.05*W) {
            double K, J_H;
            K = 4 * 0.15 * 0.15 / U;
            J_H = 0.25 * U;
            Main_Lanczos(L, N_e, t, K, U, J_H, 150, R);

            /*std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (L % 4 == 0) ? 0 : 1));
            Hamil->Diagonalization();
            Hamil->show_ground_state();*/
            out << "U/W = " << U / W << " done!" << endl;
            U += 0.1 * W;
        }
    }
    inline void Main_X(int L, int N_e, double t, double K, double X, int R) {
        double J_H = 0.0 * X;
        while (J_H <= 0.4 * X) {
            Main_Lanczos(L, N_e, t, K, 0, J_H, (L == 4) ? 100 : 150, R);
            out << "J_H/X = " << J_H / X << " done!" << endl;
            J_H += 0.01 * X;
        }
    }

    inline void Lanczos_convergence(int L, int N_e) {
        double t = 0.5, U = 2.1, K = 4 * 0.15 * 0.15 / U, J_H = 0.25 * U;
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, 500));
        vec random_vec = Create_Random_vec(Hamil->N);
        Hamil->Lanczos_convergence(random_vec);
    }
    inline void size_scaling(double t, double K, double U, double J_H, int M, int random_steps) {
        Main_Lanczos(2, 3, t, K, U, J_H, M, random_steps);
        Main_Lanczos(4, 6, t, K, U, J_H, M, random_steps);
        Main_Lanczos(6, 9, t, K, U, J_H, M, random_steps);
        Main_Lanczos(8, 12, t, K, U, J_H, M, random_steps);
        Main_Lanczos(10, 15, t, K, U, J_H, M, random_steps);
        Main_Lanczos(12, 18, t, K, U, J_H, M, random_steps);
        //Main_Lanczos(14, 21, t, K, U, J_H, M, random_steps);
    }
    inline void Energu_U_sweep(int L, int N_e, double t) {
        double U = 0.1 * W;
        ofstream file("Energy U sweep_L=" + std::to_string(L) + ".txt");
        while (U <= 3 * W) {
            double K = 4 * 0.15 * 0.15 / U;
            double J_H = 0.25 * U;

            std::unique_ptr<HamiltonianKH> Hamil(new HamiltonianKH(L, N_e, t, U, K, J_H, (L % 4 == 0) ? 0 : 1));
            Hamil->Diagonalization();
            double E0 = Hamil->eigenvalues(0);
            int degeneracy = 0;
            /*for (int k = 0; k < 10; k++) {
                out << Hamil->eigenvalues(k) << ' ';
            }*/
            while (abs(Hamil->eigenvalues(0) - Hamil->eigenvalues(degeneracy + 1)) < 1e-3) degeneracy++;
            for (int k = 0; k < 10; k++) {
                file << U / W << "\t\t" << Hamil->eigenvalues(k) - Hamil->eigenvalues(0) << endl;
            }
            out << "U/W = " << U / W << " has " << degeneracy << " ground state degeneracy" << endl;
            U += 0.1 * W;
        }
        file.close();

    }

    inline void Find_Cv_peaks(int L) {
        double U = 0.1 * W;
        double Jh = 0.0 * W;
        std::ofstream peaks("Cv_peaks3.txt");
        peaks << "U / W\t\t T1\t\tCv(T1)\t\t\tT2\t\tCv(T2)\t\t\t..." << endl;
        //for (U = 0.1 * W; U <= 3 * W; U += 0.1*W) {
        for (U = 0.15; U <= 0.4; U += 0.05) {
            stringstream Ustr;
            Ustr << setprecision(2) << fixed << U;
            std::ifstream file("Cv_L=" + std::to_string(L) + "_U=" + Ustr.str() + ".txt");
            //std::ifstream file("Cv_L=" + std::to_string(L) + "_U=3W_Jh=" + Ustr.str() + ".txt");
            std::vector<float> temperature, Cv;
            //std::ifstream file("Cv_L=10_U=2.10.txt");
            double a, b;
            std::string s;
            file >> a >> b;
            while (std::getline(file, s)) {
                temperature.push_back(a);
                Cv.push_back(b);
                file >> a >> b;
            }
            file.close();
            std::vector<int> peak_index;
            Peaks::findPeaks(Cv, peak_index);
            //peaks << U / W << "\t\t";
            peaks << U << "\t\t";
            for (int i = 0; i < peak_index.size(); i++) {
                peaks << temperature[peak_index[i]] << "\t\t" << Cv[peak_index[i]] << "\t\t";
            }
            peaks << endl;
        }
        peaks.close();
    }

    inline void Sq_max_map(int L, int lanczos_steps, int random_steps, double t) {
        std::ofstream Sqmax_map("Sq_maximum_L=" + std::to_string(L) + ".txt");
        double U;
        for (int k = L; k >= 0; k--) {
            int N_e = L + k;
            for (U = 0.05 * W; U <= 3.05 * W; U += 0.05 * W) {
                double K, J_H;
                K = 4 * 0.15 * 0.15 / U;
                J_H = 0.25 * U;
                std::unique_ptr<Lanczos> Hamil(new Lanczos(L, N_e, t, U, K, J_H, (N_e % 2 == 0) ? 0 : 1, lanczos_steps));
                if (T > 0) {
                    Hamil->Diagonalization();
                    E0 = Hamil->eigenvalues(0); // set global GS energy
                    //out << "Ground state E0 = " << E0 << endl;
                }
                vec Sq(L + 1, fill::zeros);
                if (T == 0)
                    Sq = Hamil->Sq_T0(random_steps);
                else {
                    // THIS SHIT IS ONLY FOR  T > 0
                    out << "Doing T > 0" << endl;
                    Sq = Hamil->Sq_lanczos(random_steps, T) / Hamil->Z_constT;
                }
                int l = arma::index_max(Sq);
                double q_max = (double)l * pi / ((double)L + 1.0);
                Sqmax_map << 2 - (double)N_e / (double)L << "\t\t" << U / W << "\t\t" << q_max << std::endl;
                Sqmax_map.flush();
                out << " U = " << U / W << "W done!\n";
            }
            out << "\n filling n=" << 2 - (double)N_e / (double)L << " done" << endl << endl;
        }
        Sqmax_map.close();
    }


}
#endif
