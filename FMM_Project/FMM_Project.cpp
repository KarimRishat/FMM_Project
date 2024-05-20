#pragma once
#include <vector>
#include <omp.h>
#include <fstream>
#include <string>
#include <chrono>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Factory.h"
#include "Transformations.h"
using namespace Eigen;
using point_t = Eigen::dcomplex;


size_t FindPotentialsDirect(size_t m, size_t nq, unsigned char P)
{
	double eps{ 1e-1 };
	Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	BigAdjacencyFactory adjfactory{ m, domain };
	Factory factory{ adjfactory, nq };

	Calculate_FMM::Solver fmm_solver{ factory, 1 };

	SortedData data{ factory.get_sources() };

	Eigen::VectorXcd expected(data.interval_ids.back());

	expected.setZero();

	auto start_time = std::chrono::high_resolution_clock::now();

	//for (size_t target = 0; target < data.q.size(); ++target)
	//{
	//	point_t potential{ 0.0 };

	//	for (size_t source = 0; source < data.q.size(); ++source)
	//	{
	//		if (target != source)
	//			potential += (std::log(data.point[target] - data.point[source])) * data.q[target];
	//	}
	//	expected[target] = potential;
	//}
	for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
	{
		size_t start_id{ factory.grid.far_factory.cell_intervals[cell_id] };

		for (size_t source_id = 0; source_id < factory.grid.far_factory.cell_count[cell_id]; ++source_id)
		{
			size_t far_id = factory.grid.far_factory.cell_ids[start_id + source_id];

			auto matrix = fmm_solver.CreateMatrix(cell_id, far_id);

			Map<VectorXd> q_source(data.q.data() + data.interval_ids[far_id],
				data.interval_count[far_id]);

			expected.segment(data.interval_ids[cell_id], data.interval_count[cell_id])
				+= matrix * q_source;

		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	return duration.count();

}



VectorXcd FindPotentials(size_t m, size_t nq, unsigned char P)
{
	double eps{ 1e-1 };
	Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	BigAdjacencyFactory adjfactory{ m, domain };
	Factory factory{ adjfactory, nq };

	Calculate_FMM::Solver fmm_solver{ factory, P };

	VectorXcd result = fmm_solver.SumAllPotential();

	return result;
}

void OptimalSplit(unsigned char P)
{
	double c_step{ 1.0 / 3.0 };				//step for pow
	for (double c = c_step; c < 0.5; c+=c_step)
	{
		std::string folder = "test_data/ChooseP/m=";
		std::string filename = folder + "N^{" + std::to_string(c) +"} P =" + std::to_string(P) + ".txt";
		std::ofstream outFile(filename); // Open file for writing
		if (!outFile)
		{
			std::cerr << "Error opening file: " << filename << std::endl;
			continue;
		}
		outFile << "N t(ms)" << std::endl;
		size_t r{ 1 };						//number of tries
		//m - grid size, nq - number of q in cell, real_n - m^2 * nq
		for (size_t k = 1, N = 10, m, nq, real_N; k < 12; k++)
		{
			//m = static_cast<size_t>(std::pow(N, c));
			//m = (size_t)(std::pow(N, c));
			m = 1;
			nq = N / (m*m);
			if (nq == 0) nq = 1;
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
			t_sum = (t_sum / r) /*/ (std::pow(N, 4.0/3.0))*/;
			outFile << real_N << " " << t_sum << std::endl;
			N *= 2;
		}
		outFile.close(); // Close the file
	}
}

void DirectTime(unsigned char P)
{
	
	double c_step{ 1.0 / 3.0 };				//step for pow
	std::string folder = "test_data/test_direct/";
	std::string filename = folder + "N^(" + std::to_string(c_step) + ") Eigen P =" + std::to_string(P) + ".txt";
	std::ofstream outFile(filename); // Open file for writing
	if (!outFile)
	{
		std::cerr << "Error opening file: " << filename << std::endl;
	}
	outFile << "N t(ms)" << std::endl;
	size_t r{ 1 };						//number of tries
	//m - grid size, nq - number of q in cell, real_n - m^2 * nq
	for (size_t k = 1, N = 10, m, nq, real_N; k < 15; k++)
	{
		//m = static_cast<size_t>(std::pow(N, c));
		m = (size_t)(std::pow(N, c_step));
		nq = N / (m * m);
		if (nq == 0) nq = 1;
		real_N = m * m * nq;
		size_t t_sum{ 0 };				//time sum for avg time
		for (size_t i = 0; i < r; i++)
		{
			t_sum += FindPotentialsDirect(m, nq, P);
		}
		t_sum = (t_sum / r);
		outFile << real_N << " " << t_sum << std::endl;
		N *= 2;
	}
	outFile.close(); // Close the file
}

VectorXcd RealPotential(size_t m, size_t nq )
{
	double eps{ 1e-1 };
	Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	BigAdjacencyFactory adjfactory{ m, domain };
	Factory factory{ adjfactory, nq };

	Calculate_FMM::Solver fmm_solver{ factory, 1 };

	SortedData data{ factory.get_sources() };

	Eigen::VectorXcd expected(data.interval_ids.back());

	expected.setZero();

	for (size_t target = 0; target < data.q.size(); ++target)
	{
		point_t potential{ 0.0 };

		for (size_t source = 0; source < data.q.size(); ++source)
		{
			if (target != source)
				potential += (std::log(data.point[target] - data.point[source])) * data.q[target];
		}
		expected[target] = potential;
	}
	return expected;
}

double MaxError(VectorXcd a, VectorXcd b)
{
	double max = 0.0;
	for (size_t i = 0; i < a.size(); ++i)
	{
		max = std::abs(a(i).real() - b(i).real()) > max ? std::abs(a(i).real() - b(i).real()) : max;
	}
	return max;
}


void FindError(unsigned char P)
{
	double c_step{ 1.0 / 3.0 };				//step for pow
	std::string folder = "test_data/Error/";
	std::string filename = folder +  " P =" + std::to_string(P) + ".txt";
	std::ofstream outFile(filename); // Open file for writing
	if (!outFile)
	{
		std::cerr << "Error opening file: " << filename << std::endl;
	}
	outFile << "N err" << std::endl;
	size_t r{ 1 };						//number of tries
	//m - grid size, nq - number of q in cell, real_n - m^2 * nq
	for (size_t k = 1, N = 10, m, nq, real_N; k < 12; k++)
	{
		//m = static_cast<size_t>(std::pow(N, c));
		m = (size_t)(std::pow(N, c_step));
		nq = N / (m * m);
		if (nq == 0) nq = 1;
		real_N = m * m * nq;
		double error{ 0.0 };				//time sum for avg time
		error = MaxError(FindPotentials(m, nq, P), RealPotential(m, nq));
		outFile << real_N << " " << error << std::endl;
		N *= 2;
	}
	outFile.close(); // Close the file
}


int main()
{	/*for (size_t p = 6; p < 7; p++)
	{
		DirectTime(p);
	}*/
	for (size_t p = 10; p < 20; p++)
	{
		OptimalSplit(p);
	}
	/*for (size_t P = 10; P < 15; P++)
	{
		FindError(P);
	}*/
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