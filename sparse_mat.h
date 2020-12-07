#ifndef __MySpMat_H
#define __MySpMat_H
#include "headers.h"

//sparse matrix multithreaded class, especially for sparse hamiltonian, 2xN matrix where N is number of elements
template <typename eT> class  MySpMat {
public:
	ull_int num_of_elem;
	ull_int dim;
	std::unique_ptr<std::vector<eT>> elements; //offdiagonal elements start after diagonal (up to dim(H) is diag)
	std::unique_ptr<std::vector<std::complex<ull_int>>> col_shift;
public:
	MySpMat(ull_int& n);
	//MySpMat(const std::unique_ptr<MySpMat>& SpMat); //copy constructor
	~MySpMat();

	void resize(ull_int& new_size);
	//void herm_trans(); // hermitian conjungate
	//T& operator()(const ull_int& row, const ull_int& col);
	eT& operator()(const ull_int& row, const ull_int& col);
	arma::Col<eT> operator*(arma::Col<eT>& vec);

};

template<typename eT>
std::ostream& operator<<(std::ostream& os, MySpMat<eT>& mat);


template <typename eT> class SpMat_elem {
private:
	MySpMat<eT>& parent;
	const ull_int row;
	const ull_int col;
	
	inline SpMat_elem(MySpMat<eT>& parent, const ull_int& row, const ull_int& col);

	friend MySpMat<eT>;
public:
	inline SpMat_elem<eT>& operator=(const eT value);
	inline SpMat_elem<eT>& operator+=(const eT value);
	//inline SpMat_elem<eT>& operator-=(const eT value); TODO
	//inline SpMat_elem<eT>& operator*=(const eT value); TODO
	//inline SpMat_elem<eT>& operator/=(const eT value); TODO

};



#include "sparse_mat.cpp"
#endif 