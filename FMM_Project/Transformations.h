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

		////Matrix P*P*size(farfactory.cell_intervals)
		//MatrixXcd T_ifo_container;

		const Factory* Grid_with_data;

		unsigned char P;

		// Makes the T_ifo operators for all cells in field
		void Fill_T_ifo_container()
		{
			const auto& far_cells = Grid_with_data->grid.far_factory;

			for (size_t cell_id = 0, far_cell_id, count_far, start_id; 
				cell_id < Grid_with_data->grid.cell_centers.size(); ++cell_id)
			{
				start_id = far_cells.cell_intervals[cell_id];
				count_far = far_cells.cell_count[cell_id];
				for (size_t i = 0; i < count_far; ++i)
				{
					far_cell_id = far_cells.cell_ids[start_id + i];
					T_ifo_container.middleCols
					(P * (start_id + i), P) = Fill_Tifo(cell_id, far_cell_id);
				}
			}
		}

		// Makes the incoming from outgoing translation operator for cell
		MatrixXcd Fill_Tifo(size_t target_id, size_t source_id)
		{
			const auto& cell_centers = Grid_with_data->grid.cell_centers;

			//(c_sigma - c_tau)
			point_t c_delta = cell_centers[source_id] - cell_centers[target_id];

			MatrixXcd T_ifo(P, P);

			T_ifo(0, 0) = std::log(-c_delta/*cell_centers[target_id] - cell_centers[source_id]*/);

			T_ifo(0, 1) = (-1.0) / c_delta;

			T_ifo(1, 0) = T_ifo(0, 1);

			T_ifo(1, 1) = T_ifo(0, 1)/c_delta;

			//First 2 rows
			for (size_t j = 2; j < P; j++)
			{
				T_ifo(0, j) = (-1.0) * T_ifo(0, j - 1) / c_delta;
				T_ifo(1, j) = ((point_t)j / (point_t)(j - 1)) * (-1.0) * T_ifo(1, j - 1) / c_delta;
			}

			for (size_t p = 2; p < P; p++)
			{
				T_ifo(p, 0) = ((point_t)(p - 1) / (point_t)p) * T_ifo(p - 1, 0) / c_delta;
				T_ifo(p, 1) = T_ifo(p - 1, 1) / c_delta;

				for (size_t j = 2; j < P; j++)
				{
					T_ifo(p, j) = (-1.0) * ((point_t)(p + j - 1) / (point_t)(j - 1)) * T_ifo(p, j - 1) / c_delta;
				}
			}

			return T_ifo;

		}
	public:

		MatrixXcd T_ifo_container;

		Incoming_translate_operator(const Factory &factory, unsigned char P)
			:Grid_with_data{ &factory }, P{ P }
		{
			/*T_ifo_container = MatrixXcd::Zero(P, 
				P * P * (Grid_with_data->grid.far_factory.cell_intervals.size() + 1));*/
			T_ifo_container = MatrixXcd::Zero(P,
				P * (Grid_with_data->grid.far_factory.cell_intervals.back()));
			Fill_T_ifo_container();
		}

		//Matrix P*far_cells_count of target cell
		auto T_ifo(size_t target_id) const
		{
			auto a = P * Grid_with_data->grid.far_factory.cell_intervals[target_id];
			auto b = P * Grid_with_data->grid.far_factory.cell_count[target_id];
			return T_ifo_container.block(0, a, P, b);
			/*return MatrixXcd::Map(
				T_ifo_container.data() + P * Grid_with_data->grid.far_factory.cell_intervals[target_id],
				P, P * Grid_with_data->grid.far_factory.cell_count[target_id]);*/
		}


		VectorXcd Incoming_expansion()
		{

			TranslateOperator tofs{ Grid_with_data->get_sources(),P };

			auto sources = tofs.Outgoing_expansion();

			size_t cells_count = Grid_with_data->grid.cell_centers.size();

			VectorXcd result(cells_count * P);

			for (size_t cell_id = 0, far_cells_count; cell_id < cells_count; ++cell_id)
			{

				far_cells_count = Grid_with_data->grid.far_factory.cell_count[cell_id];

				MatrixXcd t_ifo = T_ifo(cell_id);

				VectorXcd sum_vector = VectorXcd::Zero(P);

				size_t far_interval = Grid_with_data->grid.far_factory.cell_intervals[cell_id];

				for (size_t local_far_id = 0; local_far_id < far_cells_count; ++local_far_id)
				{
					//std::cout << "\nsum vector:\n" << sum_vector << std::endl;
					size_t far_id = Grid_with_data->grid.far_factory.cell_ids[far_interval + local_far_id];
					/*std::cout << std::endl << "tifo\n" << t_ifo.block(0, P * local_far_id, P, P) << std::endl;
					std::cout << std::endl << "sources\n" << sources.segment(P * far_id, P) << std::endl;*/
					sum_vector += t_ifo.block(0, P * local_far_id, P, P) * sources.segment(P * far_id, P);
				}
				//std::cout << std::endl<<"incoming\n" << sum_vector << std::endl;
				result.segment(cell_id * P, P) = sum_vector;
			}
			return result;
		}
	};


	//Translate targets from incoming operator
	class Target_translate_operator
	{

	private:
		SortedData data; // sorted field of charges
		VectorXcd incoming_expansion; //u^tau
		unsigned char P; // the error

		//Matrix size(q)*P
		MatrixXcd T_tfi_container;


		// Makes the targets from incoming operator for cell
		MatrixXcd Fill_Ttfi(size_t cell_id)
		{
			MatrixXcd T_tfi(data.interval_count[cell_id], P);

			size_t start_id{ data.interval_ids[cell_id] };

			size_t end_id{ data.interval_ids[cell_id + 1ull] };

			point_t center_tau{ data.cell_center(cell_id) };
			for (size_t source_local = 0, source_global = start_id;
				source_global < end_id; ++source_local, ++source_global)
			{
				T_tfi(source_local,0) = 1.0;
				auto delta = data.point[source_global] - center_tau;
				for (size_t p = 1; p < P; ++p)
					T_tfi(source_local,p) = T_tfi(source_local,p-1) * delta;
			}
			return T_tfi;
		}


		// Makes the T_ofs operators for all cells in field
		void Fill_T_tfi_container()
		{
			for (size_t cell_id = 0; cell_id < data.interval_count.size(); ++cell_id)
			{
				/*T_tfi_container.middleCols
				(data.interval_ids[cell_id], data.interval_count[cell_id]) = Fill_Ttfi(cell_id);*/
				T_tfi_container.block(data.interval_ids[cell_id], 0,
					data.interval_count[cell_id], P) = Fill_Ttfi(cell_id);
			}
		}

	public:
		Target_translate_operator(const SortedData& data, unsigned char P, const VectorXcd& v) : data{ data }, P{ P }, incoming_expansion{v}
		{
			T_tfi_container = MatrixXcd::Zero(data.interval_ids.back(),P);
			Fill_T_tfi_container();
		}

		auto T_tfi(size_t cell_id) const
		{
			/*return MatrixXcd::Map(
				T_tfi_container.data() + P * data.interval_ids[cell_id],
				P, data.interval_count[cell_id]);*/
			return T_tfi_container.block(data.interval_ids[cell_id], 0, data.interval_count[cell_id], P);
		}

		//u_1 - final potential from far cells

		VectorXcd Final_potential()
		{
			VectorXcd result(data.q.size());

			result.setZero();

			for (size_t cell_id = 0; cell_id < data.interval_count.size(); ++cell_id)
			{
				size_t start_id{ data.interval_ids[cell_id] };
				size_t n{ data.interval_count[cell_id] };
				MatrixXcd a = T_tfi(cell_id);
				VectorXcd b = incoming_expansion.segment(cell_id * P, P);
				VectorXcd c = a * b;
				result.segment(start_id, n) = c;
			}
			return result;
		}
	};

	class Solver
	{
		const Factory* Grid_with_data;

		SortedData data;

		unsigned char P;

		VectorXcd FindOwnPotential(size_t cell_id)
		{
			Map<VectorXd> q_target(data.q.data() + data.interval_ids[cell_id],
				data.interval_count[cell_id]);
			return CreateMatrix(cell_id, cell_id) * q_target;
		}

		VectorXcd FindNeighbourPotential(size_t cell_id)
		{
			VectorXcd result = VectorXcd::Zero(data.interval_count[cell_id]);

			size_t start_id{ Grid_with_data->grid.adjacency_factory.cell_intervals[cell_id] };

			size_t neighb_count{ Grid_with_data->grid.adjacency_factory.cell_count[cell_id] };

			for (size_t source_local = 0, source_global; source_local < neighb_count; ++source_local)
			{
				source_global = Grid_with_data->grid.adjacency_factory.cell_ids[start_id + source_local];

				auto matrix = CreateMatrix(cell_id, source_global);

				Map<VectorXd> q_source(data.q.data() + data.interval_ids[source_global], data.interval_count[source_global]);

				result += matrix * q_source;
			}

			return result;

		}

	public:
		Solver(const Factory& factory, unsigned char P) :
			Grid_with_data{ &factory }, P{ P }, data{ factory.get_sources() } {}

		MatrixXcd CreateMatrix(size_t target_id, size_t source_id)
		{
			size_t rows{ data.interval_count[target_id] }, cols{ data.interval_count[source_id] },
				start_id_target{ data.interval_ids[target_id] }, start_id_source{ data.interval_ids[source_id] };

			MatrixXcd A(rows, cols);

			for (size_t row_id = 0; row_id < rows; ++row_id)
			{
				for (size_t col_id = 0; col_id < cols; ++col_id)
				{
					if (target_id == source_id && col_id == row_id)
						A(row_id, col_id) = 0;
					else
					{
						A(row_id, col_id) = log(data.point[start_id_target + row_id] - data.point[start_id_source + col_id]);
					}
				}
			}

			return A;

		}


		VectorXcd SumAllPotential()
		{

			Incoming_translate_operator t_ifo{ *Grid_with_data, P };

			auto incoming_vector = t_ifo.Incoming_expansion();

			Target_translate_operator t_tfi{ data, P, incoming_vector };

			auto farPot = t_tfi.Final_potential();

			VectorXcd ownPot = VectorXcd::Zero(data.q.size());

			VectorXcd neighbPot = VectorXcd::Zero(data.q.size());

			size_t cell_count = Grid_with_data->grid.cell_centers.size();

			for (size_t cell_id = 0; cell_id < cell_count; ++cell_id)
			{
				ownPot.segment(data.interval_ids[cell_id], data.interval_count[cell_id])
					= FindOwnPotential(cell_id);
			}

			for (size_t cell_id = 0; cell_id < cell_count; ++cell_id)
			{
				neighbPot.segment(data.interval_ids[cell_id], data.interval_count[cell_id])
					= FindNeighbourPotential(cell_id);
			}

			return farPot + ownPot + neighbPot;
		}


	};

}
