#pragma once
#include <complex>
#include <cmath>
#include <algorithm>
//#include <Eigen/Core>
//#include <Eigen/Dense>

#include "Data.h"

//std::vector<double> T_outgoing_source(const SortedData& data, unsigned char P)
//{}

std::vector<double> T_outgoing_source(const SortedData& data, unsigned char P, size_t cell_id){

	size_t start_id{data.interval_ids[cell_id]};
	size_t end_id{data.interval_ids[cell_id+1ull]};

	
	std::vector<double> q_sigma;
	q_sigma.reserve(P);
	double sum{0};
	for (auto q_i : data.q) {
		sum += q_i;
	}
	q_sigma.push_back(sum);
	auto [x_c, y_c]{Find_centre(data)};
//	std::complex<double> c_sigma{Find_centre(data)};

//	Eigen::MatrixXd A{Eigen::MatrixXd::Zero(P, data.q.size())};
//	Eigen::VectorXd q{Eigen::VectorXd::Zero(data.q.size())};
//	A(i,j) = 1.0;
//	q(i) = data.q[i];
//	return A*q;
	
	for (size_t q_id = start_id; i < end_id; ++q_id) 
	{
		double sum{0.0};
		for(unsigned char p = 1; p < P; ++p)
		{
			//double a{1.0};
			//double k{z - c_sigma};
			//a = a*k/p*(p-1);

			//	std::complex<double> z{data.x[q_id], data.y{q_id]};
			
			std::pair<double, double> point{ data.x[q_id], data.y[q_id] };
			sum += std::pow(FindDistance(point, centre), p) * (-1.0 / p) * data.q[q_id];
		}
		q_sigma.push_back(sum);
	}
	return q_sigma;
}

std::pair<double,double> Find_centre(SortedData data) {
	double c_x = ((*std::max_element(data.x.begin(), data.x.end())) + *std::min_element(data.x.begin(), data.x.end()))*0.5;
	double c_y = ((*std::max_element(data.y.begin(), data.y.end())) + *std::min_element(data.y.begin(), data.y.end())) * 0.5;
	return std::pair<double, double>{c_x, c_y};
}

double FindDistance(std::pair<double, double> x, std::pair<double, double> y)
{
	return std::sqrt((x.first - y.first) * (x.first - y.first) + (x.second - y.second));
}
