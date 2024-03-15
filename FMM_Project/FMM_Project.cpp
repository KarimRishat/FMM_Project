#pragma once
#include <vector>

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Factory.h"
#include "Transformations.h"

int main()
{
	//Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	//AdjacencyFactory adjfactory{2ull, domain};

	//Factory factory{ adjfactory, 2ull };

	//SortedData data{ factory.get_sources() };

	//using namespace Eigen;

	////ArrayXd zrr{ ArrayXd::Zero(10) };

	//using namespace Calculate_FMM;

	//TranslateOperator Translate_op{ data };

	//VectorXcd v{ Translate_op.T_outgoing_from_source_cell(1ull, 3) };

	//std::cout << v;

	//Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	//AdjacencyFactory adjfactory{ 2ull, domain };
	//Factory factory{ adjfactory, 2ull };

	//SortedData data{ factory.get_sources() };

	//unsigned char P = 5;

	//Calculate_FMM::TranslateOperator tras_op{ data, P };

	//for (size_t cell_id = 0; cell_id < data.its_cell_center.size(); ++cell_id)
	//	for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
	//	{
	//		// check p == 0
	//	//	EXPECT_EQ(tras_op.T_ofs(cell_id)(0ull, source), 1.0);
	//		for (size_t p = 1; p < P; ++P)
	//			// per sources in a cell
	//	//		EXPECT_EQ(tras_op.T_ofs(cell_id)(p, source), -1.0 / p * std::pow(data.point[0] - data.cell_center(0), p));
	//		{
	//		}
	//	}

	using namespace Eigen;

	Domain domain{ -1.0, 1.0, -1.0, 1.0 };
	BigAdjacencyFactory adjfactory{ 1ull, domain };
	Factory factory{ adjfactory, 1ull };
	unsigned char P = 2;

	Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

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