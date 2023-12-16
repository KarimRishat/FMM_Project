#pragma once
#include <complex>
#include <cmath>
#include <algorithm>
//#include <Eigen/Core>
//#include <Eigen/Dense>

#include "Data.h"

namespace Calculate_FMM
{

	using namespace Eigen;
	using namespace DataStructs;

	class TranslateOperator
	{
		using point_t = Eigen::dcomplex;
	private:
		SortedData data; // sorted field of charges
		unsigned char P; // the error

		MatrixXcd T_ofs_container;
		 
		
		// Makes the outgoing from sources translation operator for cell
		MatrixXcd Fill_Tofs( size_t cell_id)
		{
			MatrixXcd T_ofs(P, data.interval_count[cell_id]);
			size_t start_id{ data.interval_ids[cell_id] };
			size_t end_id{ data.interval_ids[cell_id + 1ull] };

			point_t center_sigma{ data.cell_center(cell_id) };
			for (size_t source_local = 0, source_global = start_id;
				 source_global < end_id; ++source_local, ++source_global)
			{
				T_ofs(0ull, source_local) = 1.0;
				T_ofs(1ull, source_local) = -(data.point[source_global] - center_sigma);
				for (size_t p = 2; p < P; ++p)
					T_ofs(p, source_local) = T_ofs(p - 1, source_global - start_id) * 
					(p - 1.0) / (double)p * (data.point[source_global] - center_sigma);
			}
			return T_ofs;
		}



		void Fill_T_ofs_container()
		{
			for (size_t cell_id = 0; cell_id < data.interval_ids.back(); cell_id++)
				{
					MatrixXcd temp_T_ofs = Fill_Tofs(cell_id);
					T_ofs_container.block
						(P, data.interval_count[cell_id], 0, data.interval_ids[cell_id]) = temp_T_ofs;
				}
		}


		// void Fill_sources_q(VectorXd& q_source, size_t start_id, size_t end_id)
		//{
		////	auto q_sources_temp{ data.cell_sources(cell_id) };

		//	size_t count = end_id - start_id + 1;
		//	for (size_t i = 0; i < count; i++)
		//	{
		//		q_source(i) = data.q[start_id + i];
		//	}
		//}

	public:

		TranslateOperator(const SortedData &data, unsigned char P) : data{data}, P{P}
		{
			T_ofs_container = MatrixXcd::Ones(P, data.interval_ids.size());
			Fill_T_ofs_container();
		}

		auto T_ofs(size_t cell_id) const
		{
			return MatrixXcd::Map(
				T_ofs_container.data() + P * data.interval_ids[cell_id],
				P, data.interval_count[cell_id]);
		}

		// VectorXcd T_outgoing_from_source_field()
		//{
		//	for (size_t cell_id = 0; cell_id < data.interval_ids.back(); cell_id++)
		//	{
		//		VectorXcd q_sigma{ T_outgoing_from_source_cell(cell_id)};
		//	}

		//
		//}

		// makes compact representation of the sources
		//VectorXcd T_outgoing_from_source_cell(size_t cell_id)
		//{

		//	size_t start_id{data.interval_ids[cell_id]};
		//	size_t end_id{data.interval_ids[cell_id + 1ull]};

		//	/*point_t center_sigma{data.cell_center(cell_id)};*/

		//	MatrixXcd T_ofs{MatrixXcd::Zero(P, end_id - start_id)};

		//	Fill_Tofs(T_ofs, cell_id);

		//	return T_ofs(cell_id) * data.cell_sources(cell_id);
		//}
	};

}
