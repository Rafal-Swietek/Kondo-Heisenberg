#pragma once
#include "headers.h"

// ________________QR DECOMPOSITION using Givens roatations (for tridiagonal matrix)__________________
inline mat GivensMatrix(unsigned int i, unsigned int j, double s, double c, size_t n) {
	mat G = eye(n, n);
	G(i, i) = c; G(j, j) = s;
	G(i, j) = s; G(j, i) = -s;
	return G;
}
inline std::tuple<mat, mat> QR_decompose_triDiag(mat& matrix) {
	size_t n = matrix.n_cols;
	mat R = matrix, Q = eye(n, n);
	for (int j = 0; j < n - 1; j++) {//cols
		int k = j + 1; // tridiagonal
		double r = std::sqrt(R(j, j) * R(j, j) + R(k, j) * R(k, j));
		double c = R(j, j) / r;
		double s = -R(k, j) / r;
		mat G = GivensMatrix(k, j, s, c, n);
		R = G * R;
		Q = Q * G.t(); // for real this is only transpose
	}
	return std::make_tuple(Q, R);
}
inline void eig_QR_Lanczos(vec& eigenVal, mat& eigenVec, mat& matrix) {
	size_t n = matrix.n_cols;
	eigenVal = vec(n, fill::zeros);
	eigenVec = eye(n, n);
	mat A = matrix; //copy
	while (true) {
		vec tmp = eigenVal;
		auto [Q, R] = QR_decompose_triDiag(A); // O(3n^2) memory -- O(n^2) complexzity
		A = R * Q;
		eigenVec = eigenVec * Q;
		eigenVal = arma::diagvec(A);
		out << min(eigenVal) << endl;
		//out << arma::norm(eigenVal - tmp) << endl;
		if (abs(min(eigenVal) - min(tmp)) < 1e-10) break;
	}
}

// algorithm in this paper:
//  "The LL^T and QR methods for symmetric tridiagonal matrices"
//  James M. Ortega and Henry F. Kaiser
inline void eig_QR_TriDiag(vec& eigenVal, mat& matrix) {
	size_t n = matrix.n_cols;
	double u = 0, s = 0, s_prev = 0;
	vec b = arma::square(arma::diagvec(matrix, 1));
	vec a = arma::diagvec(matrix);
	while (true) {
		vec tmp = a;
		for (int k = 0; k < n; k++) {
			double gamma = a(k) - u;
			double b_iprev = (k >= 1) ? b(k - 1) : 0;
			double p = (s != 1) ? gamma * gamma / (1 - s) : (1 - s_prev) * b_iprev * b_iprev;
			double b_i = (k == n - 1) ? 0 : b(k);
			if (k > 0) b(k - 1) = s * (p + b_i);
			s_prev = s;
			s = b_i / (p + b_i);
			double a_inext = (k == n) ? 0 : a(k);
			u = s * (gamma + a_inext);
			a(k) = gamma + u;
		}
		out << min(a) << endl;
		if (abs(min(a) - min(tmp)) < 1e-10) break;
	}
	eigenVal = a; // set eigenvalues
	sort(eigenVal.begin(), eigenVal.end());
	out << eigenVal(0);
	//here go to calculating eigenvectors
}




//___________________________________________________________________________

// sturm sequence for eigen_decomposition
inline void Thompson(mat& Mat, vec& result, vec& RHS) {
	size_t n = Mat.n_cols;
	mat copy_mat = Mat;
	vec RHS_copy = RHS;
	result = vec(n, fill::zeros);
	for (int k = 1; k < n; k++) {
		double w = Mat(k, k - 1) / Mat(k - 1, k - 1);
		copy_mat(k, k) = Mat(k, k) - w * Mat(k - 1, k);
		RHS_copy(k) = RHS(k) - w * RHS(k - 1);
	}
	result(n - 1) = RHS_copy(n - 1) / copy_mat(n - 1, n - 1);
	for (int k = n - 2; k >= 0; k--)
		result(k) = (RHS_copy(k) - copy_mat(k, k + 1) * result(k + 1)) / copy_mat(k, k);
}
inline double Sturm_sequence(mat& Mat, double x, unsigned int k) {
	if (k == 0) return 1;
	else if (k == 1) return Mat(0, 0) - x;
	else return (Mat(0, k - 1) - x) * Sturm_sequence(Mat, x, k - 1) - Mat(1, k - 1) * Mat(1, k - 1) * Sturm_sequence(Mat, x, k - 2);
}
inline std::tuple<int, double> SturmPolynomial(const mat& Mat, const double x) {
	size_t n = Mat.n_cols;
	double pn = Mat(0, 0) - x, pn1 = 1;
	int sign_change = ((pn1 * pn < 0) || pn == 0) ? 1 : 0;
	for (int j = 1; j < n; j++) {
		double pn2 = pn1; pn1 = pn;
		pn = (Mat(j, j) - x) * pn1 - Mat(j, j - 1) * Mat(j, j - 1) * pn2;
		if ((pn1 * pn < 0) || pn == 0) sign_change++;
	}
	return std::make_tuple(sign_change, pn);
}
inline std::tuple<double, double, double, double> SetContainEigenVal(mat& Mat, unsigned int k) {
	double min = Mat(0, 0) - abs(Mat(1, 0));
	double max = Mat(0, 0) + abs(Mat(1, 0));
	size_t n = Mat.n_cols;
	for (int j = 1; j < n - 1; j++) {
		double tmp = Mat(j, j) - abs(Mat(j + 1, j)) - abs(Mat(j, j - 1));
		if (min > tmp) min = tmp;
		tmp = Mat(j, j) + abs(Mat(j + 1, j)) + abs(Mat(j, j - 1));
		if (max < tmp) max = tmp;
	}
	min = std::min(min, Mat(n - 1, n - 1) - abs(Mat(n - 1, n - 2)));
	max = std::max(max, Mat(n - 1, n - 1) + abs(Mat(n - 1, n - 2)));
	auto [k_min, poly_min] = SturmPolynomial(Mat, min);
	auto [k_max, poly_max] = SturmPolynomial(Mat, max);
	while (k_max > k || k_min < k - 1) {
		double mid = (min + max) / 2.;
		auto [k_mid, poly_mid] = SturmPolynomial(Mat, mid);
		if (k_mid >= k) {
			k_max = k_mid;
			poly_max = poly_mid;
			max = mid;
		}
		else {
			k_min = k_mid;
			poly_min = poly_mid;
			min = mid;
		}
	}
	return std::make_tuple(min, poly_min, max, poly_max);
}
inline std::tuple<double, double> bisect_eig(mat& Mat, unsigned int k) {
	auto [min, poly_min, max, poly_max] = SetContainEigenVal(Mat, k);
	if (poly_min == 0)
		max = min;
	while (max - min > tolerance) {
		double mid = (min + max) / 2.;
		auto [k_mid, poly_mid] = SturmPolynomial(Mat, mid);
		if (poly_mid == 0) {
			max = mid; min = mid;
		}
		else {
			if (poly_mid * poly_max < 0) {
				min = mid; poly_min = poly_mid;
			}
			else {
				max = mid; poly_max = poly_mid;
			}
		}
	}
	return std::make_tuple(min, max);
}
inline void eig_sturm(vec& eigenVal, mat& eigenVec, mat& Mat) {
	vec result;
	size_t n = Mat.n_cols;
	eigenVal = vec(n);
	eigenVec = mat(n, n);
	for (unsigned int k = 0; k < n; k++) {
		auto [min, max] = bisect_eig(Mat, k);
		double mu = (max - min) / 2.;
		vec w = Create_Random_vec(n);
		mat tmp = Mat - mu * eye(n, n);
		Thompson(tmp, result, w);
		double mu_prev = mu;
		mu = mu + 1 / arma::dot(result, w);
		while (abs(mu - mu_prev) > 1e-4) {
			w = result / arma::norm(result);
			Thompson(tmp, result, w);
			mu_prev = mu;
			mu = mu + 1 / arma::dot(result, w);
		}
		eigenVal(k) = mu;
		eigenVec.col(k) = w;
	}
}