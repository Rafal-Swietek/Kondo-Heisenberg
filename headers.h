#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <algorithm>
// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
//#define ARMA_NO_DEBUG // only after checking if no out_of_bounds, or memory outusage
//-------
#include <armadillo>
#include <cassert> // assert terminates program
#include <omp.h>
#include <time.h>
#include <numeric>
#include <utility> // auto, etc. 
#include <memory> // smart ptr
#include <thread>
#include<queue>
#include<mutex>
#include<condition_variable>
#include<functional>
#include<future>
#include "persistence1d.hpp"

using namespace std;
using namespace arma;

typedef unsigned long long int u64;
typedef std::complex<double> cpx;
typedef std::unique_ptr<std::vector<u64>> my_uniq_ptr;

// User makros
#define im cpx(0.0,1.0)
#define out std::cout << std::setprecision(16) << std::fixed
#define num_of_threads 16

#define memory_over_performance false // optimized by size --true-- (memory usage shortage) or performance --false--
#define show_system_size_parameters false // this parameter defines whether to print data such as system size for each object conctructor call
#define calculate_stand_dev true // enables the calculation of the standard deviation in the FTLM randomized method
#define calculate_Sq false
#define calulate_X0 true 
#define calculate_Cv false
#define calculate_entropy false
#define Sz_symmetry false
#define N_symmetry true
#define use_reorthonormalization (memory_over_performance)? false : false // enables in lanczos procedure full reorthogonalization - needs full krylov_space access


extern double pi;
extern double T; // temperature for Sq calculations
extern double dT; // temperature increment
extern double T_end; // temperature range (dT, T_end)
extern double domega;
extern double eta; // control parameter of imaginary part in DOS
extern double E0; // ground state energy
constexpr double W = 2.1; // kinetic energy bandwidth (for t_00=0.5)
constexpr auto tolerance = 1e-2; // ? do not remember what this is for, better leave it untouched ??

extern bool PBC; // allow periodic boundary conditions?
extern int model;

//----------------------------------------------------------------------------------------------
//--------------------------------------------------TOOLS---------------------------------------
//--source: https://github.com/PaulRitaldato1/ThreadPool  --------------------------------------
inline void print_base_vector(std::vector<int>& base_vector, std::ostream& out_str) {
	out_str << " |";
	for (auto it = base_vector.begin(); it != base_vector.end(); ++it)
		out_str << *it;
	out_str << ">";
}
inline std::vector<double> prepare_parameterVec(double _min, double _max, double step) {
	std::vector<double> parameter_vec;
	double temp = _min;
	while (temp <= _max) {
		parameter_vec.emplace_back(temp);
		temp += step;
	}
	return parameter_vec;
}

inline u64 binary_search(my_uniq_ptr& arr, u64 l_point, u64 r_point, u64 element) {
	if (r_point >= l_point) {
		u64 middle = l_point + (r_point - l_point) / 2;
		if (arr->at(middle) == element) return middle;
		else if (arr->at(middle) > element) return binary_search(arr, l_point, middle - 1, element);
		else return binary_search(arr, middle + 1, r_point, element);
	}
	//u64 j = findElement(arr, element);
	out << "Element " << element << " not present in the array" << endl;
	assert(false);
	return -1;
}
inline void int_to_binary(u64 idx, std::vector<int>& vec) {
	u64 temp = idx;
	for (int k = 0; k < vec.size(); k++) {
		vec[vec.size() - 1 - k] = static_cast<int>(temp % 8);
		temp = static_cast<u64>((double)temp / 8.);
	}
}
inline u64 binary_to_int(vector<int>& vec) {
	u64 val = 0;
	for (int k = 0; k < vec.size(); k++) {
		val += vec[vec.size() - 1 - k] * (u64)std::pow(8, k);
	}
	return val;
}

inline vec Create_Random_vec(u64 N) {
	vec random_vec(N, fill::zeros);
	double norm = 0;
	for (u64 j = 0; j < N; j++) {
		random_vec(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
		norm += random_vec(j) * random_vec(j);
	}
	return random_vec / norm;
}

/**
 * Overriding the ostream operator for pretty printing vectors.
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
	if (vec.size() != 0) {
		std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(os, " "));
		os << vec.back();
	}
	return os;
}
namespace Peaks {
	const float EPS = 2.2204e-16f;
	void findPeaks(vector<float> x0, vector<int>& peakInds);
}

//----------------------------------------------------------------------------------------------
//---------------------------------------ThreadPool class---------------------------------------
//----------------------------------------------------------------------------------------------

/*
class ThreadPool
{
public:

	ThreadPool()
	{
		m_shutdown.store(false, std::memory_order_relaxed);
		//m_shutdown = false;
		createThreads(1);
	}

	ThreadPool(std::size_t numThreads)
	{
		m_shutdown.store(false, std::memory_order_relaxed);
		//m_shutdown = false;
		createThreads(numThreads);
	}

	~ThreadPool()
	{
		m_shutdown.store(true, std::memory_order_relaxed);
		//m_shutdown = true;
		m_notifier.notify_all();

		for (std::thread& th : m_threads)
		{
			th.join();
		}
	}

	//add any number of arguments of function to queue
	template <typename Func, typename... Args>
	auto enqueue(Func&& f, Args&&... args)
	{
		//get return type of the function
		using RetType = std::invoke_result_t<Func, Args...>;

		auto task = std::make_shared<std::packaged_task<RetType()>>([&f, &args...]() { return f(std::forward<Args>(args)...); });

		{
			// lock jobQueue mutex, add job to the job queue 
			std::scoped_lock<std::mutex> lock(m_jobMutex);

			//place the job into the queue
			m_jobQueue.emplace([task]() {
				(*task)();
				});
		}
		m_notifier.notify_one();

		return task->get_future();
	}

	// utility functions
	std::size_t getThreadCount() const {
		return m_threads.size();
	}

private:

	using Job = std::function<void()>;
	std::vector<std::thread> m_threads;
	std::queue<Job> m_jobQueue;
	std::condition_variable m_notifier;
	std::mutex m_jobMutex;
	std::atomic<bool> m_shutdown;

	void createThreads(std::size_t numThreads)
	{
		m_threads.reserve(numThreads);
		for (int i = 0; i != numThreads; ++i)
		{
			m_threads.emplace_back(std::thread([this]()
				{
					while (true)
					{
						Job job;

						{
							std::unique_lock<std::mutex> lock(m_jobMutex);
							m_notifier.wait(lock, [this] {return !m_jobQueue.empty() || m_shutdown.load(std::memory_order_relaxed); });

							if (m_shutdown.load(std::memory_order_relaxed))
							{
								break;
							}

							job = std::move(m_jobQueue.front());

							m_jobQueue.pop();
						}
						job();
					}
				}));
		}
	}
}; */
/* end ThreadPool Class */