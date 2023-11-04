#pragma once
#include "Data.h"
#include <complex>
#include <cmath>
#include <algorithm>




std::vector<double> T_outgoing_source(SortedData data, int P){
	std::vector<double> q_sigma;
	q_sigma.reserve(P);
	double sum = 0;
	for (auto& q_i : data.q) {
		sum += q_i;
	}
	q_sigma.push_back(sum);
	std::pair<double, double> centre{ Find_centre(data) };
	for(size_t p = 1;p < P; ++p)
	{
		sum = 0;
		for (size_t i = 0; i < data.q.size(); ++i) {
			std::pair<double, double> point{ data.x[i],data.y[i] };
			sum += std::pow(FindDistance(point, centre), p) * (-1 / p) * data.q[i];
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