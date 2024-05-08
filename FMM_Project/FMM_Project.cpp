#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Factory.h"
#include "Transformations.h"
using namespace Eigen;
void FindPotentials(size_t m, size_t nq, unsigned char P)
{
	double eps{ 1e-1 };
	Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	BigAdjacencyFactory adjfactory{ m, domain };
	Factory factory{ adjfactory, nq };

	Calculate_FMM::Solver fmm_solver{ factory, P };

	auto result = fmm_solver.SumAllPotential();
}

void OptimalSplit(unsigned char P)
{
	double c_step{ 1.0 / 3.0 };				//step for pow
	for (double c = 0; c < 0.5; c+=c_step)
	{
		std::string folder = "test_data/";
		std::string filename = folder + "t_sum_" + std::to_string(c) + ".txt";
		std::ofstream outFile(filename); // Open file for writing
		if (!outFile)
		{
			std::cerr << "Error opening file: " << filename << std::endl;
			continue;
		}

		size_t r{ 1 };						//number of tries
		//m - grid size, nq - number of q in cell, real_n - m^2 * nq
		for (size_t k = 1, N = 10, m, nq, real_N; k < 5; k++)
		{
			m = static_cast<size_t>(std::pow(N, c));
			nq = N / m;
			real_N = m * m * nq;
			double t_sum{ 0.0 };				//time sum for avg time
			for (size_t i = 0; i < r; i++)
			{
				auto start_time = std::chrono::high_resolution_clock::now();
				FindPotentials(m, nq, P);
				auto end_time = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
				t_sum += static_cast<double>(duration.count());
			}
			t_sum = t_sum / r;
			outFile << "N: " << N << ", t_sum: " << t_sum << "ms" << std::endl;
			N *= 10;
		}
		outFile.close(); // Close the file
	}
}



int main()
{
	for (size_t p = 3; p < 5; p++)
	{
		OptimalSplit(p);
	}
	return 0;
}






/**
 * 1. Реализовать алгоритм п. 5 статьи fmm.pdf
 * 2. Для двумерного случая
 * 3. Сравнить время работы программы при разных m = {1,N^(2/3)} как функция N
 * 4. Посмотреть для разных P
 *
 * \param x
 * \param y
 * \param q
 * \param width
 * \param height
 * \param Cx
 * \param Cy
 * \param m
 * \return
 */

//std::vector<double> calculate(
//	const std::vector<std::complex<double>>& w,
//	//	const std::vector<double>& x, 
//	//	const std::vector<double>& y, 
//	const std::vector<double>& q,
//	size_t P,
//
//	double width, double height,
//	double Cx, double Cy,
//	size_t m
//)
//{
//	std::vector<double> out;
//	out.reserve(q.size());
//	out.push_back(1.0);
//	out.emplace_back(1.0);
//
//	//	out.resize(q.size());
//	//	out[0] = 1.0;
//
//	return out;
//}
//
//std::vector<double> calculate(
//	const std::vector<double>& x,
//	const std::vector<double>& y,
//	const std::vector<double>& q,
//	size_t P
//)
//{
//	// https://en.cppreference.com/w/cpp/algorithm/max_element
//	double width = (*std::max_element(x.begin(), x.end())) - *std::min_element(x.begin(), x.end());
//	double height = *std::max_element(y.begin(), y.end()) - *std::min_element(y.begin(), y.end());
//	double Cx = (*std::max_element(x.begin(), x.end()) + *std::min_element(x.begin(), x.end())) / 2.0;
//	double Cy{ (*std::max_element(y.begin(), y.end()) + *std::min_element(y.begin(), y.end())) / 2.0 };
//
//
//	std::vector<std::complex<double>> w;
//	for (size_t i = 0; i < x.size(); ++i)
//		w.emplace_back(x[i], y[i]);
//
//	calculate(w, q, P, width, height, Cx, Cy,
//		static_cast<size_t>(std::round(std::pow(q.size(), 2.0 / 3.0))));
//
//
//}
//
//void T_ofs();
//void T_ifo();
//void T_tfi();

//template<typename T>
//class Vector
//{
//	// 3*8 byte
//	size_t its_size;
//	size_t its_capacity;
//
//	T* its_data = nullptr;
//public:
//	void reserve(size_t n)
//	{
//		its_data = new T[n];
//		its_size = 0;
//	}
//	void push_back(const T& obj)
//	{
//		T[its_size] = obj;
//		++its_size;
//	}
//
//	void resize(size_t n)
//	{
//		reserve(n);
//		for (size_t i = 0; i < n; ++i)
//			T[i] = T();
//		its_size = n;
//	}
//};