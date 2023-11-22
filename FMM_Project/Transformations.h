#pragma once
#include <complex>
#include <cmath>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Data.h"

namespace Calculate_FMM {

	using namespace Eigen;

	class TranslateOperator {
	private:
		SortedData data;	//sorted field of charges
		unsigned char P;	//the error

		Eigen::MatrixXcd T_ofs_container;

		auto T_ofs(size_t cell_id) const
		{
			return Eigen::MatrixXcd::Map(
				T_ofs_container.data() + P*data.interval_ids[cell_id], 
				P, data.interval_count[cell_id]);
		}

		//std::pair<double, double> Find_centre(size_t cell_id) {
		//	size_t start_id{ data.interval_ids[cell_id] };
		//	size_t end_id{ data.interval_ids[cell_id + 1ull] };
		//	double c_x = ((*std::max_element(data.x.begin(), data.x.end())) + *std::min_element(data.x.begin(), data.x.end())) * 0.5;
		//	double c_y = ((*std::max_element(data.y.begin(), data.y.end())) + *std::min_element(data.y.begin(), data.y.end())) * 0.5;
		//	return std::pair<double, double>{c_x, c_y};
		//}

		/**
		 * \brief Determine centers of cells
		 * 
		 * \param start_id
		 * \param end_id
		 * \return 
		 */
		std::complex<double> Find_center(size_t start_id, size_t end_id) {
			double c_x = ((*std::max_element(data.x.begin() + start_id, data.x.begin() + end_id))
				+ *std::min_element(data.x.begin() + start_id, data.x.begin() + end_id)) * 0.5;
			double c_y = ((*std::max_element(data.y.begin() + start_id, data.y.begin() + end_id))
				+ *std::min_element(data.y.begin() + start_id, data.y.begin() + end_id)) * 0.5;
			return std::complex<double>{c_x, c_y};

		//	return data.cell_center(start_id);
		}
		
		//Makes the outgoing from sources translation operator
		void Fill_Tofs(MatrixXcd& T_ofs, size_t cell_id)
		{
			Fill_Tofs(
				T_ofs, 
				data.interval_ids[cell_id], 
				data.interval_ids[cell_id + 1ull], 
				data.cell_center(cell_id));
		}

		void Fill_Tofs(MatrixXcd& T_ofs, size_t start_id, size_t end_id, std::complex<double> center_sigma) {

			for (size_t source_local = 0, source_global = start_id; source_global < end_id; ++source_local, ++source_global)
			{
				T_ofs(0ull, source_local) = 1.0;
				T_ofs(1ull, source_local) = -(data.point[source_global] - center_sigma);
				for (size_t p = 2; p < P; ++p)
					T_ofs(p, source_local) = T_ofs(p - 1, source_global - start_id) * (p - 1.0) / (double)p * (data.point[source_global] - center_sigma);
			}
		}

		void Fill_sources_q(VectorXd& q_source, size_t start_id, size_t end_id) {
		//	auto q_sources_temp{ data.cell_sources(cell_id) };

			size_t count = end_id - start_id + 1;
			for (size_t i = 0; i < count; i++)
			{
				q_source(i) = data.q[start_id + i];
			}
		}


	public:
		TranslateOperator(const SortedData& data, size_t P) : 
			data{ data }, P{ P } 
		{}


		/*VectorXcd T_outgoing_from_source_field(unsigned char P)
		{
			this->P = P;
			for (size_t cell_id = 0; cell_id < data.interval_ids.back(); cell_id++)
			{
				VectorXcd q_sigma{ T_outgoing_from_source_cell(cell_id)};
			}

			
		}*/

		//makes compact representation of the sources
		VectorXcd T_outgoing_from_source_cell(size_t cell_id) {

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


			UseDataType::point_t center_sigma{ Find_center(start_id,end_id) };

			//std::complex<double> center_sigma = data.cell_center[cell_id];

			MatrixXcd T_ofs{ MatrixXcd::Zero(P, end_id - start_id ) };

			Fill_Tofs(T_ofs, start_id, end_id, center_sigma);

			return T_ofs(cell_id) * data.cell_sources(cell_id);
		}
	};

}


