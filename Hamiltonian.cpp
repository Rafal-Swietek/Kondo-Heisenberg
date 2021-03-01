#include "Hamiltonian.h"


double pi = 3.14159265358979323846;
double E0 = 0;

std::mutex my_mut;

HamiltonianKH::HamiltonianKH() {
    this->L = 8;
    this->num_of_electrons = 12;
    this->Sz = 0;
    this->U = W;
    this->K = 4 * 0.15 * 0.15 / U;
    this->t = 0.5;
    this->J_H = U / 4;
}
HamiltonianKH::HamiltonianKH(int L, int num_of_electrons, double t, double U, double K, double J_H, double Sz) {
	this->L = L; //number of sites
	this->num_of_electrons = num_of_electrons; //number of electrons in lower orbital
	this->t = t; this->U = U; 
	this->K = K;
	this->J_H = J_H;
    this->Sz = Sz;

    this->mapping = std::make_unique<std::vector<u64>>();
    //this->mapping = std::vector<u64>();
	generate_mapping();
	this->N = mapping->size();
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    Hamiltonian();

    if (show_system_size_parameters) {
        out << "Hamiltonian complete" << endl;
        out << "size(H) = " << (sizeof(H) + H.n_elem * sizeof(double)) / std::pow(10, 9) << " gb" << "\n";
    }
}
//-------------------------

void HamiltonianKH::update_parameters(double t, double U, double K, double J_H, double Sz) {
    this->t = t; this->U = U, this->K = K; this->J_H = J_H;
    this->Sz = Sz;
    this->mapping.reset(new vector<u64>());
    //this->mapping = std::vector<u64>();
    generate_mapping();
    this->N = mapping->size();
    Hamiltonian();
}

void HamiltonianKH::setHamiltonianElem(u64& k, double value, std::vector<int>&& temp){
    u64 idx = binary_search(mapping, 0, N - 1, binary_to_int(temp)); // findElement(std::move(mapping), binary_to_int(std::move(temp))); //
        H(idx, k) += value;
        H(k, idx) += value;
}
void HamiltonianKH::Hamiltonian() {
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
                            setHamiltonianElem(k, (PBC == 1 && j == L - 1) ? +t : -t, std::move(temp));
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
                                setHamiltonianElem(k, (PBC == 1 && j == L - 1) ? +t : -t, std::move(temp));
                            }
                            else  setHamiltonianElem(k, (PBC == 1 && j == L - 1) ? -t : +t, std::move(temp));
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
}

u64 findElement(std::vector<u64>* vector, u64 element) {
    std::vector<u64>::iterator it = find(vector->begin(), vector->end(), element);
    assert(it != vector->end() && "Element not present in the array");
    return (u64)std::distance(vector->begin(), it);
}


//generates the vector, which maps the base_vector index to the index in given subblock
std::tuple<int, int, int> calculateSpinElements(int L, u64& j, std::vector<int>& temp) {
    int bSz = 0; //bosonic total spin - spin of upper orbital locked to n=1 filling
    int fSz = 0; //fermionic total spin
    int N_e = 0; // numer of electrons in given state
    int_to_binary(j, temp);

    for (int k = 0; k < L; k++) {
        if (temp[k] < 4) bSz += 1;
        else bSz -= 1;
        if (temp[k] % 4 == 1) {
            fSz += 1;
            N_e += 1;
        }
        else if (temp[k] % 4 == 2) {
            fSz -= 1;
            N_e += 1;
        }
        else if (temp[k] % 4 == 3)
            N_e += 2;
    }

    return std::make_tuple(bSz, fSz, N_e);
}
void HamiltonianKH::mapping_kernel(u64 start, u64 stop, my_uniq_ptr& map_threaded, int _id) {
    int n = 1;
    //out << "A new thread joined tha party! from " << start << " to " << stop << endl;
    std::vector<int> temp(L);
    bool statement;
    for (u64 j = start; j < stop; j++) {
        int bSz = 0, fSz = 0, N_e = 0;
        std::tie(bSz, fSz, N_e) = calculateSpinElements(this->L, j, temp);
        switch (model) {
        case 0: 
            if(Sz_symmetry) statement = (bSz + fSz == this->Sz) && (N_e == this->num_of_electrons); //gKH 
            else statement = (N_e == this->num_of_electrons); //gKH - no spin-symmetry block
            break;
        case 1: statement = (bSz == this->Sz) && (N_e == 0); // Heisenberg
            break;
        case 2: 
            if(Sz_symmetry) statement = (bSz == -L) && (fSz == this->Sz) && (N_e == this->num_of_electrons); // Hubbard
            else statement = (bSz == -L) && (N_e == this->num_of_electrons); // Hubbard no spin symmetry
            break;
        }
        if(statement)
             map_threaded->push_back(j);
    }
}
void HamiltonianKH::generate_mapping() {
    u64 start = 0, stop = (u64)std::pow(8, L);
    if(num_of_threads == 1)
        mapping_kernel(start, stop, mapping, 0);
    else {
        //Threaded
        //std::vector<my_uniq_ptr> map_threaded(num_of_threads);
        std::vector<my_uniq_ptr> map_threaded(num_of_threads);
        std::vector<std::thread> threads;
        threads.reserve(num_of_threads);
        for (int t = 0; t < num_of_threads; t++) {
            start = std::pow(8, L) / (double)num_of_threads * t;
            stop = ((t + 1) == num_of_threads ? (u64)std::pow(8, L) : u64(std::pow(8, L) / (double)num_of_threads * (double)(t + 1)));
            //map_threaded[t] = my_uniq_ptr(new std::vector<u64>());
            map_threaded[t] = std::make_unique<std::vector<u64>>();
            //map_threaded[t]->reserve(std::pow(2, L)); ??
            threads.emplace_back(&HamiltonianKH::mapping_kernel, this, start, stop, ref(map_threaded[t]), t);
        }
        for (auto& t : threads) t.join();

        for (auto& t : map_threaded)
            mapping->insert(mapping->end(), std::make_move_iterator(t->begin()), std::make_move_iterator(t->end()));
        mapping->shrink_to_fit();
        //sort(mapping->begin(), mapping->end());
    }
    if (show_system_size_parameters) {
        out << "Mapping generated with  " << mapping->size() << "  elements" << endl;
        out << "Last element = " << mapping->at(mapping->size() - 1) << endl;
    }
    //out << mapping[0] << " " << mapping[mapping->size() - 1] << endl;
    assert(mapping->size() > 0 && "Not possible number of electrons - no. of states < 1");
}
//----------------------------------------------------

void HamiltonianKH::Diagonalization() {
	try {
		arma::eig_sym(eigenvalues, eigenvectors, H);
        this->ground_state = eigenvectors.col(0);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
        out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
        	assert(false);
    }
    //out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
}

vec HamiltonianKH::Total_Density_of_states(std::vector<double>& omega_vec) {
    vec resultDOS(omega_vec.size());

    double maximum = 0;
#pragma omp parallel for shared(omega_vec, resultDOS) num_threads(num_of_threads)
    for (int w = 0; w < omega_vec.size(); w++) {
        double omega = omega_vec[w];
        double DOS = 0;
    //#pragma omp parallel for shared(omega_vec, resultDOS) reduction(+: DOS)
        for (u64 n = 0; n < N; n++)
            DOS += -1. / (double)L / pi * cpx(1. / (omega + eta * 1i - eigenvalues(n))).imag();

        resultDOS(w) = DOS;
    }
    return resultDOS;
}

vec HamiltonianKH::static_structure_factor(double T) { /// here sth wrong?., check
    vec static_structure_factor(L + 1, fill::zeros);
#pragma omp parallel for num_threads(L+1)
    for (int l = 0; l <= L; l++) {
        double q = (double)l * pi / ((double)L + 1.0);
        cx_vec cpx_vec, a;
        std::vector<int> vect(L);
        sp_cx_mat Sq(sp_mat(N, N), sp_mat(N, N));
        for (int p = 0; p < N; p++) {
            int_to_binary(mapping->at(p), vect);
            Sq(p, p) = 0;
            for (int m = 0; m < L; m++) {
                double Szm = 0;
                if (vect[m] < 4) Szm += 0.5;
                else Szm -= 0.5;
                if (vect[m] % 4 == 1) Szm += 0.5;
                else if (vect[m] % 4 == 2) Szm -= 0.5;
                Sq(p, p) += std::exp(cpx(1i * q * (m + 0.0))) * Szm;
            }
        }
        cpx Sq0 = 0;
        if (T > 0) {
            for (int n = 0; n < N; n++) {
                cpx_vec = cx_vec(eigenvectors.col(n), vec(N, fill::zeros));
                a = Sq * cpx_vec;
                Sq0 += std::exp(-eigenvalues(n) / T) * cdot(a, a);
            }
        }
        else {
            cpx_vec = cx_vec(eigenvectors.col(0), vec(N, fill::zeros));
            a = Sq * cpx_vec;
            Sq0 = cdot(a,a);
        }
        /*Sq0 = 0;
        mat cor_mat = correlation_matrix();
        for (int m = 0; m < L; m++) {
            for (int k = 0; k < L; k++) {
                Sq0 += std::exp(1i * q * (m - k + 0.0)) * cor_mat(m,k);
            }
        }*/
        static_structure_factor(l) = real(2.0 * Sq0 / pi / (L + 1.0));
    }
    return static_structure_factor;
}
double HamiltonianKH::partition_function(double T) {
    double Z = 0;
    if (T > 0) {
        for (int n = 0; n < eigenvalues.size(); n++)
            Z += std::exp(-(eigenvalues(n) - E0) / T);
    }
    else
        Z = 1;
    return Z;
}

void HamiltonianKH::show_ground_state() {
    stringstream Ustr, Jhstr, nstr;
    Ustr << setprecision(2) << fixed << U / W;
    Jhstr << setprecision(2) << fixed << J_H / ((U == 0) ? 1.0 : U);
    nstr << setprecision(2) << fixed << (double)this->num_of_electrons / (double)this->L;
    ofstream GSfile;
    GSfile.open("GS_L=" + std::to_string(L) + "__U=" + Ustr.str() + "W__Jh=" + Jhstr.str() + "U_n=" + nstr.str() + "_PBC=" + std::to_string(PBC) + ".txt");
    vec GS((u64)std::pow(2, L), fill::zeros);
    std::vector<int> base_vector(L);
    //ostream& stream = GSfile;
    for (u64 k = 0; k < N; k++) {
        int_to_binary(mapping->at(k), base_vector);
        int val = 0;
        for (int j = 0; j < L; j++)
            val += (1 - int(base_vector[base_vector.size() - 1 - j] / 4)) * (int)std::pow(2, j);
        GS(val) += ground_state(k) * ground_state(k);
    }
    //GS = arma::abs(GS);
    GS = GS / arma::norm(GS); //normalizing to be sure
    GSfile << endl;
    double maximum = arma::max(GS);
    std::vector<int> vec(L);
    GSfile << "State\t\t probability\n\n";
    for (u64 k = 0; k < GS.size(); k++) {
        if (std::fabs(GS(k)) >= 0.05 * maximum) {
            u64 temp = k;
            for (int p = 0; p < L; p++) {
                vec[vec.size() - 1 - p] = temp % 2;
                temp = static_cast<int>((double)temp / 2.);
            }
            print_base_vector(vec, GSfile);
            GSfile << "\t\t" << GS(k) * GS(k) << endl;
        }
    }
    GSfile.close();
}
mat HamiltonianKH::correlation_matrix() {
    mat cor_mat(L, L, fill::zeros);
    vector<int> vect(L);
    for (int p = 0; p < N; p++) {
        int_to_binary(mapping->at(p), vect);
        for (int m = 0; m < L; m++) {
            double Szm = 0;
            if (vect[m] < 4) Szm += 0.5;
            else Szm -= 0.5;
            if (vect[m] % 4 == 1) Szm += 0.5;
            else if (vect[m] % 4 == 2) Szm -= 0.5;

            for (int k = 0; k < L; k++) {
                double Szk = 0;
                if (vect[k] < 4) Szk += 0.5;
                else Szk -= 0.5;
                if (vect[k] % 4 == 1) Szk += 0.5;
                else if (vect[k] % 4 == 2) Szk -= 0.5;
                cor_mat(m, k) += ground_state(p) * Szm * Szk * ground_state(p);
                //cor_mat(m, k) += ground_state(p) * Szk * ground_state(p);
            }
        }
    }
    return cor_mat;
}
double HamiltonianKH::total_spin_squared(vec&& state) {
    double S2 = 0;
    u64 idx;
    vector<int> vect(L), temp(L);
    for (int p = 0; p < N; p++) { //! 1 -localized, 2 - electrons
        int_to_binary(mapping->at(p), vect);
        for (int m = 0; m < L; m++) {
            double Sz1m = 0, Sz2m = 0;

            if (vect[m] < 4) Sz1m = 0.5;
            else Sz1m = -0.5;
            if (vect[m] % 4 == 1) Sz2m = 0.5;
            else if (vect[m] % 4 == 2) Sz2m = -0.5;

            S2 += state(p) * (Sz1m * Sz1m + Sz1m * Sz2m + Sz2m * Sz1m + Sz2m * Sz2m) * state(p); //  <Sz(i) Sz(i)>
            for (int k = m; k < L; k++) {
                double Sz1k = 0, Sz2k = 0;

                if (vect[k] < 4) Sz1k = 0.5;
                else Sz1k = -0.5;
                if (vect[k] % 4 == 1) Sz2k = 0.5;
                else if (vect[k] % 4 == 2) Sz2k = -0.5;

                S2 += 2 * state(p) * (Sz1m * Sz1k + Sz1m * Sz2k + Sz2m * Sz1k + Sz2m * Sz2k) * state(p); //  <Sz(i) Sz(j)> dla j>i
                // <S^+ S^-> + <S^- S^+>
                temp = vect;
                if ((vect[m] < 4) && (vect[k] >= 4)) {
                    temp[m] += 4;
                    temp[k] -= 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    S2 += 2*state(idx) * state(p);
                }
                // <S^+ s^-> + <S^- s^+>
                temp = vect;
                if (vect[m] >= 4 && vect[k] % 4 == 1) {
                    temp[m] -= 4;
                    temp[k] += 1;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    S2 += 2*state(idx) * state(p);
                }
                // <s^+ S^-> + <s^- S^+>
                temp = vect;
                if (vect[m] % 4 == 2 && vect[k] < 4) {
                    temp[m] -= 1;
                    temp[k] += 4;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    S2 += 2*state(idx) * state(p);
                }
                // <s^+ s^-> + <s^- s^+>
                temp = vect;
                if ((vect[m] % 4 == 1) && (vect[k] % 4 == 2)) {
                    temp[m] += 1;
                    temp[k] -= 1;
                    idx = binary_search(mapping, 0, N - 1, binary_to_int(temp));
                    S2 += 2*state(idx) * state(p);
                }
            }
        }
    }
    return S2;
}

void HamiltonianKH::printEnergy(double Ef) {
    ofstream Efile;
    stringstream Ustr, Nstr;
    Ustr << setprecision(1) << fixed << U;
    Nstr << setprecision(2) << fixed << (double)num_of_electrons / (double)L;
    Efile.open("E_L=" + std::to_string(L) + ",N_e="+std::to_string(num_of_electrons)+"_U=" + Ustr.str() + ".txt");

    for (int k = 0; k < 20; k++)
        Efile << U << "\t\t" << eigenvalues(k) - Ef << endl;

    Efile.close();

}

//----------------------------------------------------
//----------------------------------------------------

// Find Peaks
void diff(vector<float> in, vector<float>& output)
{
    output = vector<float>(in.size() - 1);

    for (int i = 1; i < in.size(); ++i)
        output[i - 1] = in[i] - in[i - 1];
}
void vectorProduct(vector<float> a, vector<float> b, vector<float>& output)
{
    output = vector<float>(a.size());

    for (int i = 0; i < a.size(); ++i)
        output[i] = a[i] * b[i];
}
void findIndicesLessThan(vector<float> in, float threshold, vector<int>& indices)
{
    for (int i = 0; i < in.size(); ++i)
        if (in[i] < threshold)
            indices.push_back(i + 1);
}
void selectElements(vector<float> in, vector<int> indices, vector<float>& output)
{
    for (int i = 0; i < indices.size(); ++i)
        output.push_back(in[indices[i]]);
}
void selectElements(vector<int> in, vector<int> indices, vector<int>& output)
{
    for (int i = 0; i < indices.size(); ++i)
        output.push_back(in[indices[i]]);
}
void signVector(vector<float> in, vector<int>& output)
{
    output = vector<int>(in.size());

    for (int i = 0; i < in.size(); ++i)
    {
        if (in[i] > 0)
            output[i] = 1;
        else if (in[i] < 0)
            output[i] = -1;
        else
            output[i] = 0;
    }
}


void Peaks::findPeaks(vector<float> x0, vector<int>& peakInds)
{
    int minIdx = distance(x0.begin(), min_element(x0.begin(), x0.end()));
    int maxIdx = distance(x0.begin(), max_element(x0.begin(), x0.end()));
    float sel = (x0[maxIdx] - x0[minIdx]) / 4.0;

    int len0 = x0.size();

    vector<float> dx;
    diff(x0, dx);
    replace(dx.begin(), dx.end(), 0.0f, -Peaks::EPS);
    vector<float> dx0(dx.begin(), dx.end() - 1);
    vector<float> dx1(dx.begin() + 1, dx.end());
    vector<float> dx2;

    vectorProduct(dx0, dx1, dx2);

    vector<int> ind;
    findIndicesLessThan(dx2, 0, ind); // Find where the derivative changes sign

    vector<float> x;

    vector<int> indAux(ind.begin(), ind.end());
    selectElements(x0, indAux, x);
    x.insert(x.begin(), x0[0]);
    x.insert(x.end(), x0[x0.size() - 1]);;


    ind.insert(ind.begin(), 0);
    ind.insert(ind.end(), len0);

    int minMagIdx = distance(x.begin(), min_element(x.begin(), x.end()));
    float minMag = x[minMagIdx];
    float leftMin = minMag;
    int len = x.size();

    if (len > 2)
    {
        float tempMag = minMag;
        bool foundPeak = false;
        int ii;

        // Deal with first point a little differently since tacked it on
        // Calculate the sign of the derivative since we tacked the first
        //  point on it does not neccessarily alternate like the rest.
        vector<float> xSub0(x.begin(), x.begin() + 3);//tener cuidado subvector
        vector<float> xDiff;//tener cuidado subvector
        diff(xSub0, xDiff);

        vector<int> signDx;
        signVector(xDiff, signDx);

        if (signDx[0] <= 0) // The first point is larger or equal to the second
        {
            if (signDx[0] == signDx[1]) // Want alternating signs
            {
                x.erase(x.begin() + 1);
                ind.erase(ind.begin() + 1);
                len = len - 1;
            }
        }
        else // First point is smaller than the second
        {
            if (signDx[0] == signDx[1]) // Want alternating signs
            {
                x.erase(x.begin());
                ind.erase(ind.begin());
                len = len - 1;
            }
        }

        if (x[0] >= x[1])
            ii = 0;
        else
            ii = 1;

        float maxPeaks = ceil((float)len / 2.0);
        vector<int> peakLoc(maxPeaks, 0);
        vector<float> peakMag(maxPeaks, 0.0);
        int cInd = 1;
        int tempLoc;

        while (ii < len)
        {
            ii = ii + 1;//This is a peak
            //Reset peak finding if we had a peak and the next peak is bigger
            //than the last or the left min was small enough to reset.
            if (foundPeak)
            {
                tempMag = minMag;
                foundPeak = false;
            }

            //Found new peak that was lager than temp mag and selectivity larger
            //than the minimum to its left.

            if (x[ii - 1] > tempMag && x[ii - 1] > leftMin + sel)
            {
                tempLoc = ii - 1;
                tempMag = x[ii - 1];
            }

            //Make sure we don't iterate past the length of our vector
            if (ii == len)
                break; //We assign the last point differently out of the loop

            ii = ii + 1; // Move onto the valley

            //Come down at least sel from peak
            if (!foundPeak && tempMag > sel + x[ii - 1])
            {
                foundPeak = true; //We have found a peak
                leftMin = x[ii - 1];
                peakLoc[cInd - 1] = tempLoc; // Add peak to index
                peakMag[cInd - 1] = tempMag;
                cInd = cInd + 1;
            }
            else if (x[ii - 1] < leftMin) // New left minima
                leftMin = x[ii - 1];

        }

        // Check end point
        if (x[x.size() - 1] > tempMag && x[x.size() - 1] > leftMin + sel)
        {
            peakLoc[cInd - 1] = len - 1;
            peakMag[cInd - 1] = x[x.size() - 1];
            cInd = cInd + 1;
        }
        else if (!foundPeak && tempMag > minMag)// Check if we still need to add the last point
        {
            peakLoc[cInd - 1] = tempLoc;
            peakMag[cInd - 1] = tempMag;
            cInd = cInd + 1;
        }

        //Create output
        if (cInd > 0)
        {
            vector<int> peakLocTmp(peakLoc.begin(), peakLoc.begin() + cInd - 1);
            selectElements(ind, peakLocTmp, peakInds);
            //peakMags = vector<float>(peakLoc.begin(), peakLoc.begin()+cInd-1);
        }



    }


}
