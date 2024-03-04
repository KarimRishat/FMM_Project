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

	using point_t = Eigen::dcomplex;

	class TranslateOperator
	{
		
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
				size_t n{ data.interval_count[cell_id] };
				result.segment(cell_id * P, P) = T_ofs(cell_id) * sources.segment(start_id, n);
			}
			return result;
		}


		

	};


	class Incoming_translate_operator
	{

		//Matrix P*P*size(farfactory.cell_intervals)
		MatrixXcd T_ifo_container;

		Factory Grid_with_data;

		unsigned char P;

		// Makes the T_ifo operators for all cells in field
		void Fill_T_ifo_container()
		{
			//SortedData data{ Grid_with_data.get_sources() };

			FarFactory* far_cells = &Grid_with_data.grid.far_factory;

			size_t far_cell_id, count_far, start_id;

			for (size_t cell_id = 0; cell_id < Grid_with_data.grid.cell_centers.size(); ++cell_id)
			{
				start_id = cell_id * far_cells->cell_intervals[cell_id];
				count_far = far_cells->cell_intervals[cell_id + 1] - far_cells->cell_intervals[cell_id];
				for (size_t i = 0; i < count_far; ++i)
				{
					far_cell_id = far_cells->cell_ids[start_id + i];
					T_ifo_container.middleCols
					(P * (start_id + i), P) = Fill_Tifo(cell_id, far_cell_id);
				}
			}
		}

		// Makes the incoming from outgoing translation operator for cell
		MatrixXcd Fill_Tifo(size_t target_id, size_t source_id)
		{
			SortedData data{ Grid_with_data.get_sources() };
			//(c_sigma - c_tau)
			point_t c_delta = data.cell_center(source_id) - data.cell_center(target_id);

			MatrixXcd T_ifo(P, P);

			T_ifo(0, 0) = std::log(-c_delta);

			T_ifo(0, 1) = (-1.0) / c_delta;

			T_ifo(1, 0) = T_ifo(0, 1);

			T_ifo(1, 1) = T_ifo(0, 1);

			//First 2 rows
			for (size_t j = 2; j < P; j++)
			{
				T_ifo(0, j) = (-1.0) * T_ifo(0, j - 1) / c_delta;
				T_ifo(1, j) = (point_t)(j / (j - 1)) * (-1.0) * T_ifo(1, j - 1) / c_delta;
			}

			for (size_t p = 2; p < P; p++)
			{

				T_ifo(p, 0) = (point_t)((p - 1) / p) * T_ifo(p - 1, 0) / c_delta;
				T_ifo(p, 1) = T_ifo(p - 1, 1) / c_delta;

				for (size_t j = 2; j < P; j++)
				{
					T_ifo(p, j) = (-1.0) * (point_t)((p + j - 1) / (j - 1)) * T_ifo(p, j - 1) / c_delta;
				}
			}

			return T_ifo;

		}

		/*auto T_ifo(size_t target_id, size_t source_id) const
		{
			return MatrixXcd::Map(
				T_ifo_container.data() + P * data.interval_ids[target_id],
				P, data.interval_count[target_id]);

			return MatrixXcd::Map(
				T_ifo_container.data() + P * data.interval_ids[target_id],
				P, P);

		}*/

		/*VectorXcd Incoming_expansion()
		{
			Map<VectorXd> sources(data.q.data(), data.q.size());

			VectorXcd result(data.interval_count.size() * P);

			for (size_t cell_id = 0; cell_id < data.interval_count.size(); ++cell_id)
			{
				size_t start_id{ data.interval_ids[cell_id] };
				size_t n{ data.interval_count[cell_id] };
				result.segment(cell_id * P, P) = T_ofs(cell_id) * sources.segment(start_id, n);
			}
			return result;
		}*/
	};



}
