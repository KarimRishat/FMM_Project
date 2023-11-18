#pragma once
#include <complex>
#include <cmath>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Data.h"


namespace Calculate_FMM {

	using namespace Eigen;
	class Transform_Class {
	private:
		SortedData data;	//sorted field of charges
		unsigned char P;	//the error incurred



		//std::pair<double, double> Find_centre(size_t cell_id) {
		//	size_t start_id{ data.interval_ids[cell_id] };
		//	size_t end_id{ data.interval_ids[cell_id + 1ull] };
		//	double c_x = ((*std::max_element(data.x.begin(), data.x.end())) + *std::min_element(data.x.begin(), data.x.end())) * 0.5;
		//	double c_y = ((*std::max_element(data.y.begin(), data.y.end())) + *std::min_element(data.y.begin(), data.y.end())) * 0.5;
		//	return std::pair<double, double>{c_x, c_y};
		//}


		//makes compact representation of the sources
		VectorXcd T_outgoing_from_source(const SortedData& data, size_t cell_id) {	

			size_t start_id{ data.interval_ids[cell_id] };
			size_t end_id{ data.interval_ids[cell_id + 1ull] };

			/*std::vector<std::complex<double>> q_sigma;
			q_sigma.reserve(P);
			std::complex<double> sum{ 0.0 };
			for (size_t i = start_id; i < end_id; i++)
			{
				sum += data.q[i];
			}
			q_sigma.push_back(sum);*/
			
			//auto [x_c, y_c] {Find_centre(cell_id)};
			//	std::complex<double> c_sigma{Find_centre(data)};
			/*Eigen::MatrixXd T_ofs{ Eigen::MatrixXd::Zero(P, data.q.size()) };
			Eigen::VectorXd q{ Eigen::VectorXd::Zero(data.q.size()) };*/

			std::complex<double> center_sigma = data.cell_center[cell_id];

			MatrixXcd T_ofs{ MatrixXd::Zero(P, end_id - start_id + 1) };

			Fill_Tofs(T_ofs, start_id, end_id, P, center_sigma);	

			VectorXcd q{ VectorXd::Zero(end_id - start_id + 1) };

			Fill_sources_q(q, start_id, end_id);
			
			return T_ofs * q;
		}


		//Makes the outgoing from sources translation operator
		void Fill_Tofs(MatrixXcd& T_ofs, size_t start_id, size_t end_id, unsigned char P, std::complex<double> center_sigma) {

			size_t count = end_id - start_id + 1;

			for (size_t Tofs_col = 0; Tofs_col < count; ++Tofs_col) {
				T_ofs(0, Tofs_col) = 1.0;
			}

			for (size_t Tofs_row = 1; Tofs_row < P; ++Tofs_row)
			{
				for (size_t Tofs_col = 0; Tofs_col < count; ++Tofs_col)
				{
					std::complex<double> complex_p{ (-1) * (Tofs_row - 1) / (Tofs_row) };
					T_ofs(Tofs_row, Tofs_col) = complex_p * (data.x[Tofs_col + start_id] - center_sigma);
				}
			}
		}


		void Fill_sources_q(VectorXcd& q_source, size_t start_id, size_t end_id) {
			size_t count = end_id - start_id + 1;
			for (size_t i = 0; i < count; i++)
			{
				q_source(i) = data.q[start_id + i];
			}
		}

	public:
		Transform_Class(const SortedData& data) : data{ data } {};


		//VectorXcd T_outgoing_source(const SortedData& data, unsigned char P)
		//{}
	};

}


