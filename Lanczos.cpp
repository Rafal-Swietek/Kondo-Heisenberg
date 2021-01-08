#include "Lanczos.h"


//----------------------------------------------------------------------------------------------
//--------------------------------------------------LANCZOS-------------------------------------
Lanczos::Lanczos(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz, int lanczos_steps) {
    this->L = L; //number of sites
    this->num_of_electrons = num_of_electrons; //number of electrons in lower orbital
    this->t = t; this->U = U;
    this->K = K;
    this->J_H = J_H;
    this->Sz = Sz;
    this->lanczos_steps = lanczos_steps;

    this->mapping = std::make_unique<std::vector<u64>>();
    //this->mapping = std::vector<u64>();
    generate_mapping();
    this->N = mapping->size();
    if (show_system_size_parameters)
        out << "dim = " << N << endl;
    try {
        this->H = sp_mat(N, N); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    if (!memory_over_performance)
        Hamiltonian();
}

void Lanczos::setHamiltonianElem(u64& k, double value, std::vector<int>&& temp) {
    u64 idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));// findElement(std::move(mapping), binary_to_int(std::move(temp))); //
    assert(idx < N && "Somehow index out of scope, wtf?? Found element not possibly present in the array");
    H(idx, k) += value;
    H(k, idx) += value;
}
void Lanczos::Hamiltonian() {
    std::vector<int> base_vector(L);
    std::vector<int> temp(base_vector);
    //#pragma omp parallel for num_threads(num_of_threads)
    for (u64 k = 0; k < N; k++) {
        int_to_binary(mapping->at(k), base_vector);
        temp = base_vector;
        int s_i, s_j; //i=j, j=j+1
        bool PBC = 0; //allows periodic boundary conditions if =1
        int next_j;
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;
                // Diagonal spin part
                if (base_vector[j] < 4) s_i = 1;
                else s_i = 0;

                // PBC = i+1 : (L-1)+1 = 0
                if (base_vector[next_j] < 4) s_j = 1;
                else s_j = 0;
                H(k, k) += K * (s_i - 0.5) * (s_j - 0.5);

                //Kinetic spin part: S+ S-
                temp = base_vector;
                if (s_i == 0 && s_j == 1) { // S_i^+ S_i+1^-
                    temp[j] = base_vector[j] - 4; //spin filp
                    temp[next_j] = base_vector[next_j] + 4;
                    setHamiltonianElem(k, K / 2., std::move(temp));
                }
                // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                    //spin up
                temp = base_vector;
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                    temp[next_j] -= 1; // anihilate spin-up electron
                    temp[j] += 1; // create spin-up electron
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        setHamiltonianElem(k, -t, std::move(temp));
                    }
                    else  setHamiltonianElem(k, +t, std::move(temp));
                }
                //spin down
                temp = base_vector;
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        temp[next_j] -= 2; // anihilate spin-down electron
                        temp[j] += 2; // create spin-down electron
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            setHamiltonianElem(k, -t, std::move(temp));
                        }
                        else  setHamiltonianElem(k, +t, std::move(temp));
                    }
                }
            }
            // electron repulsion
            if (base_vector[j] == 7 || base_vector[j] == 3)
                H(k, k) += U;
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            temp = base_vector;
            if (base_vector[j] == 5) {// S_i^+ s_i^-
                temp[j] = 2;
                setHamiltonianElem(k, -J_H, std::move(temp));
            }
            //Diagonal - z part
            if (base_vector[j] == 1 || base_vector[j] == 6)
                H(k, k) -= 2.0 * J_H * 0.25;
            if (base_vector[j] == 2 || base_vector[j] == 5)
                H(k, k) += 2.0 * J_H * 0.25;
        }
        //if (k % (u64)std::pow(10, 2) == 0) out << k << endl;
    }
    if (show_system_size_parameters) {
        out << "Hamiltonian complete" << endl;
        out << "size(H) = " << (sizeof(H) + H.n_nonzero * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
    }
}

void Lanczos::update_parameters(double t, double U, double K, double J_H, double Sz) {
        this->t = t; this->U = U, this->K = K; this->J_H = J_H;
        this->Sz = Sz;
        this->mapping.reset(new std::vector<u64>());
        //this->mapping = std::vector<u64>();
        generate_mapping();
        this->N = mapping->size();
        //this->H.reset(new sp_mat(N, N));
        this->H = sp_mat(N, N);
        if (!memory_over_performance)
            Hamiltonian();
}
void Lanczos::Hamil_vector_multiply_kernel(u64 start, u64 stop, vec& initial_vec, vec& result_vec_threaded) {
    std::vector<int> base_vector(L);
    int PBC = 0;
    int next_j;
    u64 idx;
    int s_i, s_j;
    for (int k = start; k < stop; k++) {
        int_to_binary(mapping->at(k), base_vector);
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;

                if (base_vector[j] < 4) s_i = 1;
                else s_i = 0;
                if (base_vector[next_j] < 4) s_j = 1;
                else s_j = 0;

                //localised spin diagonal part
                result_vec_threaded(k) += K * (s_i - 0.5) * (s_j - 0.5) * initial_vec(k);
                //Kinetic localised spin part: S+ S-
                if (s_i == 0 && s_j == 1) {
                    base_vector[j] -= 4; //spin filp
                    base_vector[next_j] += 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    result_vec_threaded(k) += K / 2. * initial_vec(idx); //   S_i^+ * S_i+1^-
                    base_vector[j] += 4; //spin filp
                    base_vector[next_j] -= 4;
                }
                if (s_i == 1 && s_j == 0) {
                    base_vector[j] += 4; //spin filp
                    base_vector[next_j] -= 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    result_vec_threaded(k) += K / 2. * initial_vec(idx); //   S_i^+ * S_i+1^-
                    base_vector[j] -= 4; //spin filp
                    base_vector[next_j] += 4;
                }
                //---------------------

                // electron hopping
                //spin up
                //only odd numbers have up-electrons  //even numbers lack one up-electron
                if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) { //  j+1 -> j
                    base_vector[next_j] -= 1; // anihilate spin-up electron
                    base_vector[j] += 1; // create spin-up electron
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    base_vector[next_j] += 1;
                    base_vector[j] -= 1;
                    if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                        result_vec_threaded(k) += -t * initial_vec(idx);
                    }
                    else
                        result_vec_threaded(k) += t * initial_vec(idx);
                }
                if (base_vector[j] % 2 == 1 && base_vector[next_j] % 2 == 0) { // j -> j+1
                    base_vector[next_j] += 1; // anihilate spin-up electron
                    base_vector[j] -= 1; // create spin-up electron
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                    base_vector[next_j] -= 1;
                    base_vector[j] += 1;
                    if (base_vector[j] % 2 == 1 && base_vector[next_j] % 4 == 2) { // sign change when up-electron hop over down-electron (if 3=du)
                        result_vec_threaded(k) += -t * initial_vec(idx);
                    }
                    else
                        result_vec_threaded(k) += t * initial_vec(idx);
                }
                // spin down
                // the if statement contains every possible down-electron hopping: next_j->j
                if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                    if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                        base_vector[next_j] -= 2; // anihilate spin-down electron
                        base_vector[j] += 2; // create spin-down electron
                        idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                        base_vector[next_j] += 2;
                        base_vector[j] -= 2;
                        if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                            result_vec_threaded(k) += -t * initial_vec(idx);
                        }
                        else
                            result_vec_threaded(k) += t * initial_vec(idx);
                    }
                }
                // the if statement contains every possible down-electron hopping: j->next_j
                if (base_vector[j] % 4 == 2 || base_vector[j] % 4 == 3) {
                    if (base_vector[next_j] % 4 == 0 || base_vector[next_j] % 4 == 1) {
                        base_vector[j] -= 2; // anihilate spin-down electron
                        base_vector[next_j] += 2; // create spin-down electron
                        idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                        base_vector[next_j] -= 2;
                        base_vector[j] += 2;
                        if ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 1) || ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 1)) || ((base_vector[j] % 4 == 3 && base_vector[next_j] % 4 == 0))) {
                            result_vec_threaded(k) += -t * initial_vec(idx);
                        }
                        else
                            result_vec_threaded(k) += t * initial_vec(idx);
                    }
                }
            }
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            if (base_vector[j] == 5) {// S_i^+ s_i^-
                base_vector[j] = 2;
                idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                base_vector[j] = 5;
                result_vec_threaded(k) += -J_H * initial_vec(idx);
            }
            if (base_vector[j] == 2) {// S_i^+ s_i^-
                base_vector[j] = 5;
                idx = binary_search(mapping, 0, N - 1, binary_to_int(base_vector));
                base_vector[j] = 2;
                result_vec_threaded(k) += -J_H * initial_vec(idx);
            }
            // electron repulsion
            if (base_vector[j] == 7 || base_vector[j] == 3)
                result_vec_threaded(k) += U * initial_vec(k);
            //Diagonal - z part
            if (base_vector[j] == 1 || base_vector[j] == 6)
                result_vec_threaded(k) -= 2.0 * J_H * 0.25 * initial_vec(k);
            if (base_vector[j] == 2 || base_vector[j] == 5)
                result_vec_threaded(k) += 2.0 * J_H * 0.25 * initial_vec(k);
        }
    }
}
void Lanczos::Hamil_vector_multiply(vec& initial_vec, vec& result_vec) {
    result_vec.zeros();
    std::vector<arma::vec> result_threaded(num_of_threads);
    std::vector<std::thread> threads;
    //Hamil_vector_multiply_kernel(0, N, initial_vec, result_vec);
    threads.reserve(num_of_threads);
    for (int t = 0; t < num_of_threads; t++) {
        u64 start = t * N / num_of_threads;
        u64 stop = ((t + 1) == num_of_threads ? N : N * (t + 1) / num_of_threads);
        result_threaded[t] = arma::vec(stop - start, fill::zeros);
        threads.emplace_back(&Lanczos::Hamil_vector_multiply_kernel, this, start, stop, ref(initial_vec), ref(result_vec));
    }for (auto& t : threads) t.join();
}


void Lanczos::Build_Lanczos_Hamil_wKrylov(vec& initial_vec) {
    this->Krylov_space = mat(N, 1);
    this->H_L = mat(1, 1, fill::zeros);

    Krylov_space.col(0) = initial_vec;

    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H * Krylov_space.col(0);

    double alfa = arma::dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        try {
            Krylov_space.resize(N, j + 1);
            H_L.resize(j + 1, j + 1);
        }
        catch (const bad_alloc& e) {
            std::cout << "Memory exceeded " << e.what() << "\n";
            out << "dim(H) = " << Krylov_space.size() * sizeof(Krylov_space(0, 0)) * std::pow(10, -9) << " GB" << "\n";
            assert(false);
        }

        beta = sqrt(dot(tmp, tmp));
        Krylov_space.col(j) = tmp / beta;
        vec tmp2 = Krylov_space.col(j);
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H * tmp2;

        alfa = arma::dot(Krylov_space.col(j), tmp);
        if (use_reorthonormalization) {
            vec temporary(N, fill::zeros);
            for (int k = 0; k <= j; k++)
                temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
            tmp = tmp - temporary;
            temporary.zeros();
            for (int k = 0; k <= j; k++)
                temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
            tmp = tmp - temporary;
        }
        else
            tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;
    }
}
void Lanczos::Build_Lanczos_Hamil(vec& initial_vec) {
    this->H_L = mat(lanczos_steps, lanczos_steps, fill::zeros);
    double beta = arma::dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)
    // already normalized input random vector
    this->randVec_inKrylovSpace = vec(lanczos_steps);
    randVec_inKrylovSpace(0) = arma::dot(initial_vec, initial_vec); // =1

    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H * initial_vec;

    double alfa = arma::dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;
    vec tmp2_prev = initial_vec;
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        double beta = sqrt(arma::dot(tmp, tmp));
        vec tmp2 = tmp / beta;
        randVec_inKrylovSpace(j) = dot(tmp2, initial_vec);

        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H * tmp2;

        alfa = arma::dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * tmp2_prev;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;

        tmp2_prev = tmp2;
    }
}

void Lanczos::Diagonalization() {
    //vec rand = randu<vec>(N);
    this->eigenvalues = vec(lanczos_steps, fill::zeros);
    this->eigenvectors = mat(lanczos_steps, lanczos_steps, fill::zeros);
    vec rand = Create_Random_vec(N);
    Build_Lanczos_Hamil(rand);
    eig_sym(eigenvalues, eigenvectors, H_L);
}
void Lanczos::Lanczos_GroundState() {
    vec initial_vec = Create_Random_vec(N);
    Build_Lanczos_Hamil(initial_vec);
    try {
        eig_sym(eigenvalues, eigenvectors, H_L);
    } catch (const bad_alloc& e) {
        out << "Memory error:  " << e.what() << "\n";
        //out << "dim(H) = " << H_L.size() * sizeof(H_L(0, 0)) << "\n";
        assert(false);
    }
    out << eigenvalues(0) << endl << endl;
    vec GS = eigenvectors.col(0);

    this->ground_state = vec(N, fill::zeros);

    double beta = dot(initial_vec, initial_vec);
    initial_vec = initial_vec / sqrt(beta); // normalized Krylov_space(j=0)

    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H * initial_vec;

    double alfa = dot(initial_vec, tmp);
    tmp = tmp - alfa * initial_vec;

    for (int j = 1; j < lanczos_steps; j++) {
        beta = sqrt(dot(tmp, tmp));
        vec tmp2 = tmp / beta;

        ground_state += GS(j) * tmp2;
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H * tmp2;
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        initial_vec = tmp2;
        //if (show_system_size_parameters) out << j << "lanczos" << endl;
    }
}
void Lanczos::Lanczos_convergence(vec& initial_vec) {
    this->Krylov_space = mat(N, 1);
    this->H_L = mat(1, 1, fill::zeros);
    vec e_prev(5, fill::zeros);
    ofstream convergence("Lanczos convergence.txt");
    ofstream energy("Energy Lanczos.txt");
    convergence << std::setprecision(20) << std::fixed;

    Krylov_space.col(0) = initial_vec;

    double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
    Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
    vec tmp(N, fill::zeros);
    if (memory_over_performance) Hamil_vector_multiply(initial_vec, tmp); // tmp = H * Krylov_space(0)
    else tmp = H * Krylov_space.col(0);

    double alfa = arma::dot(Krylov_space.col(0), tmp);
    tmp = tmp - alfa * Krylov_space.col(0);
    H_L(0, 0) = alfa;
    for (int j = 1; j < lanczos_steps; j++) {
        Krylov_space.resize(N, j + 1);
        H_L.resize(j + 1, j + 1);

        beta = sqrt(dot(tmp, tmp));
        Krylov_space.col(j) = tmp / beta;
        vec tmp2 = Krylov_space.col(j);
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H * tmp2;

        alfa = arma::dot(Krylov_space.col(j), tmp);
        //tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);
        vec temporary(N, fill::zeros);
        for (int k = 0; k <= j; k++)
            temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
        tmp = tmp - temporary;
        temporary.zeros();
        for (int k = 0; k <= j; k++)
            temporary += dot(tmp, Krylov_space.col(k)) * Krylov_space.col(k);
        tmp = tmp - temporary;

        H_L(j, j) = alfa;
        H_L(j, j - 1) = beta;
        H_L(j - 1, j) = beta;
        eig_sym(eigenvalues, eigenvectors, H_L);

        for (int k = 0; k < j; k++) {
            energy << j << "\t" << eigenvalues(k) << endl;
        }
        if (j >= 5) {
            vec e(e_prev.size());
            e(0) = eigenvalues(0);
            e(1) = eigenvalues(1);
            e(2) = eigenvalues(2);
            e(3) = eigenvalues(3);
            e(4) = eigenvalues(4);
            convergence << j << "\t";
            for (int k = 0; k < e.size(); k++) {
                convergence << fabs((e_prev(k) - e(k)) / e(k)) << "\t";
            }
            convergence << endl;
            e_prev = e;
        }
        out << j << endl;
    }
    energy.close();
    convergence.close();
}

vec Lanczos::Heat_Capacity_Lanczos(int random_steps) {
    vec E_av(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec E_av2(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    this->partition_function = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    vec Cv(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    auto temperature = prepare_parameterVec(dT, T_end, dT);
    for (int r = 0; r < random_steps; r++) {
        auto rand_vec = Create_Random_vec(N);
        Build_Lanczos_Hamil(rand_vec);
        eig_sym(eigenvalues, eigenvectors, H_L);
        //E_av.zeros(); E_av2.zeros(); Z.zeros();
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
                overlap *= overlap;
                Z(k) += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - E0) / T);
                E_av(k) += eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - E0) / T);
                E_av2(k) += eigenvalues(m) * eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - E0) / T);
                partition_function(k) += Z(k);
            }
        }
    }
    E_av = E_av / Z * (double)N / (double)random_steps;
    E_av2 = E_av2 / Z * (double)N / (double)random_steps;
    //vec heat_cap = (E_av2 - arma::square(E_av)) / (double)L;
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k];
        Cv(k) = (E_av2(k) - E_av(k) * E_av(k)) / T / T / (double)L;
        //heat_cap(k) = heat_cap(k) / T / T;
    }
    //Cv += heat_cap / random_steps;
    //Cv_2 += arma::square(heat_cap) / random_steps;
    return Cv;
}
vec Lanczos::static_spin_susceptibility(int random_steps) {
    this->partition_function = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    auto temperature = prepare_parameterVec(dT, T_end, dT);
    for (int r = 0; r < random_steps; r++) {
        auto rand_vec = Create_Random_vec(N);
        Build_Lanczos_Hamil(rand_vec);
        eig_sym(eigenvalues, eigenvectors, H_L);
        //vec chi_tmp(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
                overlap *= overlap;
                partition_function(k) += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - E0) / T);
                chi_0(k) += Sz * Sz * overlap * std::exp(-(eigenvalues(m) - E0) / T) * (double)N / (double)random_steps / T / (double)L;
            }
        }
        //chi_0 += chi_tmp;
    }
    return chi_0;
}
vec Lanczos::thermal_average_lanczos(vec&& quantity, int& random_steps) {
    assert(quantity.n_elem > 0 && "given object is empty");
    double T = dT;
    int t = 0;
    vec result(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    while (T <= T_end) {
        double Z = 0; //partition function
        double overlap = 0; //overlap of random vectopr and lanczos eigenvectors
        for (int m = 0; m < lanczos_steps; m++) {
            overlap = dot(randVec_inKrylovSpace, eigenvectors.col(m));
            overlap *= overlap;
            Z += (double)N / (double)random_steps * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            result(t) += quantity(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
        }
        result(t) = result(t) / Z * (double)N / (double)random_steps;
        T += dT; t++;
    }
    return result;
}
vec Lanczos::Sq_lanczos(int random_steps, double T) {
    vec Sq(L + 1, fill::zeros);
    double Z = 0;
    this->Z_constT = 0;
    vector<int> vect(L);
    for (int r = 0; r < random_steps; r++) {
        auto rand_vec = Create_Random_vec(N);
        Build_Lanczos_Hamil_wKrylov(rand_vec);
        mat V;
        eig_sym(eigenvalues, V, H_L);
        eigenvectors = Krylov_space * V;
        for (int l = 0; l <= L; l++) {
            double q = (double)l * pi / ((double)L + 1.0);
            sp_cx_mat Sq_mat(sp_mat(N, N), sp_mat(N, N));
#pragma omp parallel for num_threads(num_of_threads)
            for (int p = 0; p < N; p++) {
                int_to_binary(mapping->at(p), vect);
                Sq_mat(p, p) = 0;
                for (int m = 0; m < L; m++) {
                    double Szm = 0;
                    if (vect[m] < 4) Szm += 0.5;
                    else Szm -= 0.5;
                    if (vect[m] % 4 == 1) Szm += 0.5;
                    else if (vect[m] % 4 == 2) Szm -= 0.5;
                    Sq_mat(p, p) += std::exp(cpx(1i * q * (m + 0.0))) * Szm;
                }
            }
            cx_vec temp = Sq_mat * Sq_mat.t() * cx_vec(rand_vec, vec(N, fill::zeros));
            double Sq_tmp = 0;
            for (int m = 0; m < lanczos_steps; m++) {
                double overlap = dot(rand_vec, eigenvectors.col(m));
                Z += (double)N / (double)random_steps * overlap * overlap * std::exp(-(eigenvalues(m) - E0) / T);
                overlap = real(overlap * cdot(cx_vec(eigenvectors.col(m), vec(N, fill::zeros)), temp));
                Sq_tmp += overlap * std::exp(-(eigenvalues(m) - E0) / T);
            }
            Sq(l) += (double)N / (double)random_steps * 2.0 * Sq_tmp / pi / (L + 1.0);
        }
    }
    Z_constT = Z / (L + 1.0);
    return Sq;
}
