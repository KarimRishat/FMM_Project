#pragma once
#include <complex>
#include <cmath>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

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

		//Matrix P*size(q)
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


		// Makes the T_ofs operators for all cells in field
		void Fill_T_ofs_container()
		{
			for (size_t cell_id = 0; cell_id < data.interval_count.size(); ++cell_id)
				{
					T_ofs_container.middleCols
						(data.interval_ids[cell_id], data.interval_count[cell_id]) = Fill_Tofs(cell_id);
				}
		}

	public:

		TranslateOperator(const SortedData &data, unsigned char P) : data{data}, P{P}
		{
			T_ofs_container = MatrixXcd::Ones(P, data.interval_ids.back());
			Fill_T_ofs_container();
		}

		auto T_ofs(size_t cell_id) const
		{
			return MatrixXcd::Map(
				T_ofs_container.data() + P * data.interval_ids[cell_id],
				P, data.interval_count[cell_id]);
		}



		//q^sigma - outgoing expansion of Omega


		VectorXcd Outgoing_expansion()
		{
			Map<VectorXd> sources(data.q.data(), data.q.size());

			VectorXcd result(data.interval_count.size() * P);

			for (size_t cell_id = 0; cell_id < data.interval_count.size(); ++cell_id)
			{
				size_t start_id{ data.interval_ids[cell_id] };
				size_t n = data.interval_count[cell_id];
				result.segment(cell_id * P, P) = T_ofs(cell_id) * sources.segment(start_id, n);
			}
			return result;
		}


	};
}
