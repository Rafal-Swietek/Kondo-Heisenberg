#ifndef __MySpMat_CPP
#define __MySpMat_CPP
#include "sparse_mat.h"

template<typename eT>
MySpMat<eT>::MySpMat(ull_int& n) {
    num_of_elem = n;
    dim = n;
    elements = std::unique_ptr<std::vector<eT>>(new std::vector<eT>(n));
    col_shift = std::unique_ptr<std::vector<std::complex<ull_int>>>(new std::vector<std::complex<ull_int>>());
}

template<typename eT>
MySpMat<eT>::~MySpMat() {}


template<typename eT>
void MySpMat<eT>::resize(ull_int& new_size) {
    elements->resize(new_size);
    num_of_elem = new_size;
}

template<typename eT>
arma::Col<eT> MySpMat<eT>::operator*(arma::Col<eT>& vec) {
    Col<eT> result(vec.size(), fill::zeros);
    int sign = 1, tmp = 0;
    for (int k = 0; k < dim; k++) {
        result(k) = elements->at(k) * vec(k);
        /*while (sign > 0) {
            result(k) += elements->at(dim + tmp) * vec(abs(col_shift->at(tmp)));
            result(abs(col_shift->at(tmp))) += elements->at(dim + tmp) * vec(k);
            tmp++;
            if (tmp < col_shift->size() - 1) {
                sign = col_shift->at(tmp) / abs(col_shift->at(tmp)) * col_shift->at(++tmp) / abs(col_shift->at(++tmp));
            }
            else
                sign = -1;
        }
        if (tmp < col_shift->size() - 1) sign = 1;
        else sign = -1;*/
    }
    return result;
}

template<typename eT>
eT& MySpMat<eT>::operator()(const ull_int& row, const ull_int& col) {
    if (row != col) {
        //auto idx = std::find(col_shift->begin(), col_shift->end(), cpx(row, col));
            elements->emplace_back(0);
            col_shift->emplace_back(cpx(row, col));
            return elements->at(elements->size() - 1);
      
    }
    else
        return elements->at(row); // diagonal elements are stored first
}


template<typename eT>
std::ostream& operator<<(std::ostream& os, MySpMat<eT>& mat) {
    auto it = mat.elements->begin();
    os << "elem\t offset:" << std::endl;
    for (; it != mat.elements->end(); ++it) {
        /*if (it - mat.elements->begin() < mat.dim)
            os << (*it) << std::endl;
        else os << (*it) << "\t" << mat.col_shift->at(it - mat.elements->begin() - mat.dim) << std::endl;*/
        os << (*it) << std::endl;
    }
    return os;
}


#endif



// LOGIC OF SpMat CLASS::
// 1 2 0
// 2 5 1
// 0 1 7
// SpMat() = [1 5 7 2 1]  --  making use of diagonal mirror symmetry
