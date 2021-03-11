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
    if (this->lanczos_steps > this->N) this->lanczos_steps = this->N;
    if (show_system_size_parameters)
        out << "dim = " << N << endl;
    try {
        this->H = sp_mat(N, N); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    if (!memory_over_performance) {
        Hamiltonian();
        if (show_system_size_parameters) {
            out << "Hamiltonian complete" << endl;
            out << "size(H) = " << (sizeof(H) + H.n_nonzero * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
        }
    }
}

void Lanczos::setHamiltonianElem(u64& k, double value, std::vector<int>&& temp) {
    u64 idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
    H(idx, k) += value;
    H(k, idx) += value;
}
void Lanczos::Hamiltonian() {
    std::vector<int> base_vector(L);
    std::vector<int> temp(base_vector);
    //#pragma omp parallel for num_threads(num_of_threads)
    for (u64 k = 0; k < N; k++) {
        int_to_binary(mapping->at(k), base_vector);
        //print_base_vector(base_vector, out);
        temp = base_vector;
        int s_i, s_j; //i=j, j=j+1
        int next_j;
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;
                if (K != 0) {
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
                }
                if (t != 0 && num_of_electrons != 0) {
                    // electron hopping j+1 -> j  (j->j+1 is hermitian conjunagte)
                        //spin up
                    temp = base_vector;
                    //only odd numbers have up-electrons  //even numbers lack one up-electron
                    if (base_vector[next_j] % 2 == 1 && base_vector[j] % 2 == 0) {
                        temp[next_j] -= 1; // anihilate spin-up electron
                        temp[j] += 1; // create spin-up electron
                            if (base_vector[next_j] % 4 == 3 && base_vector[j] % 2 == 0) {
                                setHamiltonianElem(k, (PBC==1 && j==L-1)? +t : -t, std::move(temp));
                            }
                            else  setHamiltonianElem(k, (PBC == 1 && j == L - 1) ? -t : +t, std::move(temp));
                    }
                    //spin down
                    temp = base_vector;
                    // the if statement contains every possible down-electron hopping: next_j->j
                    if (base_vector[next_j] % 4 == 2 || base_vector[next_j] % 4 == 3) {
                        if (base_vector[j] % 4 == 0 || base_vector[j] % 4 == 1) {
                            temp[next_j] -= 2; // anihilate spin-down electron
                            temp[j] += 2; // create spin-down electron
                            if ((base_vector[next_j] % 4 == 3 && base_vector[j] % 4 == 1) || (base_vector[j] % 4 == 1 && base_vector[next_j] % 4 == 2)) {
                                setHamiltonianElem(k, (PBC == 1 && j == L - 1) ?  +t : -t, std::move(temp));
                            }
                            else  setHamiltonianElem(k, (PBC == 1 && j == L - 1) ?  -t : +t, std::move(temp));
                        }
                    }
                }
            }
            // electron repulsion
            if (num_of_electrons != 0) {
                if (base_vector[j] == 7 || base_vector[j] == 3)
                    H(k, k) += U;
            }
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            if (J_H != 0) {
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
    int next_j;
    u64 idx;
    int s_i, s_j;
    for (u64 k = start; k < stop; k++) {
        int_to_binary(mapping->at(k), base_vector);
        for (int j = 0; j <= L - 1; j++) {
            if (j < L - 1 || PBC == 1) {
                if (PBC == 1 && j == L - 1) next_j = 0;
                else next_j = j + 1;
                if (K != 0) {
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
                }
                // electron hopping
                if (t != 0 && num_of_electrons != 0) {
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
            }
            // electron repulsion
            if (num_of_electrons != 0) {
                if (base_vector[j] == 7 || base_vector[j] == 3)
                    result_vec_threaded(k) += U * initial_vec(k);
            }
            // electron-localised spin interaction ( interorbital electronic spin interaction)
            if (J_H != 0) {
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
                //Diagonal - z part
                if (base_vector[j] == 1 || base_vector[j] == 6)
                    result_vec_threaded(k) -= 2.0 * J_H * 0.25 * initial_vec(k);
                if (base_vector[j] == 2 || base_vector[j] == 5)
                    result_vec_threaded(k) += 2.0 * J_H * 0.25 * initial_vec(k);
            }
        }
    }
}
void Lanczos::Hamil_vector_multiply(vec& initial_vec, vec& result_vec) {
    result_vec.zeros();
    std::vector<arma::vec> result_threaded(num_of_threads);
    std::vector<std::thread> threads;
    threads.reserve(num_of_threads);
    for (int t = 0; t < num_of_threads; t++) {
        u64 start = (double)N / (double)num_of_threads * t;
        u64 stop = ((t + 1) == num_of_threads ? N : (double) N / (double)num_of_threads * (t + 1));
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
        //out << j << endl;
    }
    //out << "Watch out! It's Krylov!!" << endl;
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
        if (show_system_size_parameters) out << j << "lanczos" << endl;
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
void Lanczos::Lanczos_GroundState(vec& initial_vec) {
    this->ground_state = vec(N, fill::zeros);
    Build_Lanczos_Hamil(initial_vec);
    try {
        eig_sym(eigenvalues, eigenvectors, H_L);
    }
    catch (const bad_alloc& e) {
        out << "Memory error:  " << e.what() << "\n";
        //out << "dim(H) = " << H_L.size() * sizeof(H_L(0, 0)) << "\n";
        assert(false);
    }
    //out << eigenvalues(0) << endl << endl;
    vec GS = eigenvectors.col(0);


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

        this->ground_state += GS(j) * tmp2;
        if (memory_over_performance) Hamil_vector_multiply(tmp2, tmp); // tmp = H * tmp2
        else tmp = H * tmp2;
        alfa = dot(tmp2, tmp);
        tmp = tmp - alfa * tmp2 - beta * initial_vec;

        initial_vec = tmp2;
        if (show_system_size_parameters) out << j << "lanczos" << endl;
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
        vec rand_vec = Create_Random_vec(N);
        vec randum;
        if (use_reorthonormalization) {
            Build_Lanczos_Hamil_wKrylov(rand_vec);
            mat V;
            eig_sym(eigenvalues, V, H_L);
            eigenvectors = Krylov_space * V;
            randum = rand_vec;
        }
        else {
            Build_Lanczos_Hamil(rand_vec);
            eig_sym(eigenvalues, eigenvectors, H_L);
            randum = randVec_inKrylovSpace;
        }
        //E_av.zeros(); E_av2.zeros(); Z.zeros();
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randum, eigenvectors.col(m));
                overlap *= overlap;
                Z(k) += overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
                E_av(k) += eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
                E_av2(k) += eigenvalues(m) * eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            }
        }
    }
    Z = Z * (double)N / (double)random_steps;
    E_av = E_av / Z *(double)N / (double)random_steps;
    E_av2 = E_av2 / Z *(double)N / (double)random_steps;
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

vec Lanczos::entropy(int random_steps) {
    vec E_av(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec Z(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    vec S(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
    auto temperature = prepare_parameterVec(dT, T_end, dT);
    for (int r = 0; r < random_steps; r++) {
        vec rand_vec = Create_Random_vec(N);
        vec randum;
        if (use_reorthonormalization) {
            Build_Lanczos_Hamil_wKrylov(rand_vec);
            mat V;
            eig_sym(eigenvalues, V, H_L);
            eigenvectors = Krylov_space * V;
            randum = rand_vec;
        }
        else {
            Build_Lanczos_Hamil(rand_vec);
            eig_sym(eigenvalues, eigenvectors, H_L);
            randum = randVec_inKrylovSpace;
        }
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randum, eigenvectors.col(m));
                overlap *= overlap;
                Z(k) += overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
                E_av(k) += eigenvalues(m) * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            }
        }
    }
    Z = Z * (double)N / (double)random_steps;
    E_av = E_av / Z * (double)N / (double)random_steps;
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
    for (int k = 0; k < temperature.size(); k++) {
        double T = temperature[k];
        S(k) = 1 / (double)L * (std::log(Z(k)) + (E_av(k) - eigenvalues(0)) / T);
    }
    return S;
}

vec Lanczos::static_spin_susceptibility(int random_steps) {
    this->partition_function = vec(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    vec chi_0(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);

    auto temperature = prepare_parameterVec(dT, T_end, dT);
    for (int r = 0; r < random_steps; r++) {
        vec rand_vec = Create_Random_vec(N);
        vec randum;
        if (use_reorthonormalization) {
            Build_Lanczos_Hamil_wKrylov(rand_vec);
            mat V;
            eig_sym(eigenvalues, V, H_L);
            eigenvectors = Krylov_space * V;
            randum = rand_vec;
        }
        else {
            Build_Lanczos_Hamil(rand_vec);
            eig_sym(eigenvalues, eigenvectors, H_L);
            randum = randVec_inKrylovSpace;
        }
        //vec chi_tmp(static_cast<int>((T_end - dT) / dT + 1), fill::zeros);
#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
        for (int k = 0; k < temperature.size(); k++) {
            double overlap;
            double T = temperature[k];
            for (int m = 0; m < lanczos_steps; m++) {
                overlap = dot(randum, eigenvectors.col(m));
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
    double Z_constT;
    std::vector<sp_cx_mat> Sq_matrices(L + 1);
    int s = L + 1;
#pragma omp parallel for shared(Sq_matrices) num_threads(s)
    for (int l = 0; l <= L; l++) {
        std::vector<int> vect(L);
        double q = (double)l * pi / ((double)L + 1.0);
        sp_cx_mat Sq_mat(sp_mat(N, N), sp_mat(N, N));
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
        Sq_matrices[l] = Sq_mat;
    }
    double Sq0 = 0, Z = 0;
    for (int r = 0; r < random_steps; r++) {
        vec rand_vec = Create_Random_vec(N);
        Build_Lanczos_Hamil_wKrylov(rand_vec);
        mat V;
        eig_sym(eigenvalues, V, H_L);
        eigenvectors = Krylov_space * V;
        for (int l = 0; l <= L; l++) {
            cx_vec temp = Sq_matrices[l] * Sq_matrices[l].t() * cx_vec(rand_vec, vec(N, fill::zeros));
            double Sq_tmp = 0;
            for (int m = 0; m < lanczos_steps; m++) {
                double overlap = dot(rand_vec, eigenvectors.col(m));
                Z += (double)N / (double)random_steps * overlap * overlap * std::exp(-(eigenvalues(m) - eigenvalues(0))/ T);
                overlap = real(overlap * cdot(cx_vec(eigenvectors.col(m), vec(N, fill::zeros)), temp));
                Sq_tmp += overlap * std::exp(-(eigenvalues(m) - eigenvalues(0)) / T);
            }
            Sq0 = (double)N / (double)random_steps * 2.0 * Sq_tmp / pi / (L + 1.0);
            if (Sq0 != Sq0) Sq0 = 0; // catch NaN
            Sq(l) += Sq0;
            if (Sq0 > std::pow(10, 3)) out << "whaaaaat?" << endl;
        }
        out << "R = " << r << endl;
    }
    Z_constT = Z / ((double)L + 1.0);
    return Sq / Z_constT;
}

vec Lanczos::Sq_T0(int random_steps) {
    vec Sq(L + 1, fill::zeros);
    for (int r = 0; r < random_steps; r++) {
        vec randvec = Create_Random_vec(N);
        Lanczos_GroundState(randvec); // get the GS from lanczos procedure
        mat corr_mat = correlation_matrix(); // <GS| Sz(i)Sz(j) |GS> correlations
        int num_thr = (L + 1 > num_of_threads) ? num_of_threads : L + 1;
#pragma omp parallel for num_threads(num_thr)
        for (int l = 0; l <= L; l++) {
            double q = (double)l * pi / ((double)L + 1.0);
            cpx Sq0 = 0;
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < L; k++) {
                    Sq0 += std::exp(cpx(1i * q * double(j - k))) * corr_mat(j, k);
                }
            }
            Sq(l) += real(2.0 * Sq0 / pi / (L + 1.0));
        }
    }
    return Sq / (double)random_steps;
}

void Lanczos::SSF_T0() {
    stringstream Ustr, Jhstr, nstr;
    nstr << setprecision(2) << fixed << (double)this->num_of_electrons / (double)this->L;
    Ustr << setprecision(2) << fixed << this->U / W;
    Jhstr << setprecision(2) << fixed << this->J_H / ((this->U == 0) ? 1.0 : this->U);
    std::ofstream SpinFactorFile("SSF_L=" + std::to_string(L) + "_U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
    SpinFactorFile << std::setprecision(16) << std::fixed;
    int s = L + 1;
    double domega = 0.0005;
//#pragma omp parallel for num_threads(s)
    vec randvec = Create_Random_vec(N);
    Lanczos_GroundState(randvec);
    for (int l = 0; l <= L; l++) {
        double q;
        if(PBC) q = 2*(double)l * pi / (double)L;
        else q = (double)l * pi / ((double)L + 1.0);
        sp_cx_mat Sq(N, N);
//#pragma omp parallel for num_threads(num_of_threads)
        for (int p = 0; p < N; p++) {
            std::vector<int> vect(L);
            int_to_binary(mapping->at(p), vect);
            Sq(p, p) = 0;
            for (int m = 0; m < L; m++) {
                double Szm = 0;
                if (vect[m] < 4) Szm += 0.5;
                else Szm -= 0.5;
                //if (vect[m] % 4 == 1) Szm += 0.5;
                //else if (vect[m] % 4 == 2) Szm -= 0.5;
                Sq(p, p) += std::exp(cpx(1i * q * (m + 0.0))) * Szm;
                //Sq(p, p) += sin(q * (m + 0.0)) * Szm;
            }
        }
        cx_vec Sq_GS = Sq * this->ground_state;
        cx_double alfa = cdot(cx_vec(this->ground_state, vec(N, fill::zeros)), Sq.t()*Sq_GS);
        std::unique_ptr<Lanczos> Hamil(new Lanczos(L, num_of_electrons, t, U, K, J_H, (this->num_of_electrons % 2 == 0) ? 0 : 1, lanczos_steps));
        vec initial_vec = real(Sq_GS / sqrt(alfa));
        Hamil->Build_Lanczos_Hamil(initial_vec);
        double SpinFactor = 0;
        for (double omega = 0; omega < 1; omega += domega) { // Calculate S(q,w)
            cx_double z = omega + 1i * eta + eigenvalues(0);
            cx_double Continous_Fraction = z - Hamil->H_L(this->lanczos_steps - 1, this->lanczos_steps - 1);
            for (int m = this->lanczos_steps - 2; m >= 0; m--) {
                Continous_Fraction = z - Hamil->H_L(m, m) - Hamil->H(m, m + 1) * Hamil->H_L(m, m + 1) / Continous_Fraction;
            }
            SpinFactor = -1 / (L + 1.0) / pi * imag(alfa / Continous_Fraction);
            SpinFactorFile << q << "\t\t" << omega << "\t\t" << SpinFactor << endl;
        }
        out << q << endl;
    }
    SpinFactorFile.close();
}