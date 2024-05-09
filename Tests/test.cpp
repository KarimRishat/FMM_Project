#include "pch.h"
#include "../FMM_Project/Factory.h"
#include "../FMM_Project/Transformations.h"
#include <fstream>
#define EPS 1e-13

// https://learn.microsoft.com/ru-ru/visualstudio/test/how-to-use-google-test-for-cpp?view=vs-2022
// https://github.com/google/googletest/blob/main/docs/primer.md
// https://learn.microsoft.com/ru-ru/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022

using point_t = Eigen::dcomplex; 
using namespace Eigen;
using namespace DataStructs;

namespace DataGenerators
{
	TEST(DomainTest, DefaultState) {

		Domain domain{};

		EXPECT_EQ(domain.height, 1.0);
		EXPECT_EQ(domain.width, 1.0);
		EXPECT_EQ(domain.x_a, 0.0);
		EXPECT_EQ(domain.x_b, 1.0);
		EXPECT_EQ(domain.y_a, 0.0);
		EXPECT_EQ(domain.y_b, 1.0);
	}

	TEST(AdjacencyTest, NoDivision) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 1ull, domain };

		EXPECT_EQ(adjfactory.x_grid[0], -1.0);
		EXPECT_EQ(adjfactory.x_grid[1], 1.0);
		EXPECT_EQ(adjfactory.y_grid[0], -1.0);
		EXPECT_EQ(adjfactory.y_grid[1], 1.0);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids.size(), 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_intervals.front(), 0ull);

		EXPECT_EQ(adjfactory.adjacency_factory.cell_intervals.back(), adjfactory.adjacency_factory.cell_ids.size());
	}


	TEST(AdjacencyTest, DivideByTwo) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 2ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids.size(), 12ull);

		EXPECT_EQ(adjfactory.adjacency_factory.cell_intervals.back(), adjfactory.adjacency_factory.cell_ids.size());

		for (size_t id = 0ull; id < adjfactory.adjacency_factory.cell_intervals.size(); ++id)
			EXPECT_EQ(adjfactory.adjacency_factory.cell_intervals[id], id * 3ull);

		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[0], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[1], 2ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[2], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[3], 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[4], 2ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[5], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[6], 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[7], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[8], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[9], 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[10], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[11], 2ull);
	}

	TEST(AdjacencyTest, DivideByThree) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids.size(), 40ull);

		EXPECT_EQ(adjfactory.adjacency_factory.cell_intervals.back(), adjfactory.adjacency_factory.cell_ids.size());


		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[0], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[1], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[2], 4ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[3], 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[4], 2ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[5], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[6], 4ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[7], 5ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[8], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[9], 4ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[10], 5ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[11], 0ull);

		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[12], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[13], 4ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[14], 6ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[15], 7ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[16], 0ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[17], 1ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[18], 2ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[19], 3ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[20], 5ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[21], 6ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[22], 7ull);
		EXPECT_EQ(adjfactory.adjacency_factory.cell_ids[23], 8ull);
	}


	TEST(FarCellsTest, NoDivision) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 1ull, domain };

		EXPECT_EQ(adjfactory.x_grid[0], -1.0);
		EXPECT_EQ(adjfactory.x_grid[1], 1.0);
		EXPECT_EQ(adjfactory.y_grid[0], -1.0);
		EXPECT_EQ(adjfactory.y_grid[1], 1.0);
		EXPECT_EQ(adjfactory.far_factory.cell_ids.size(), 0ull);
		EXPECT_EQ(adjfactory.far_factory.cell_intervals.front(), 0ull);

		EXPECT_EQ(adjfactory.far_factory.cell_intervals.back(), adjfactory.far_factory.cell_ids.size());
	}


	TEST(FarCellsTest, DivideByTwo) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 2ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.far_factory.cell_ids.size(), 0ull);

		EXPECT_EQ(adjfactory.far_factory.cell_intervals.back(), adjfactory.far_factory.cell_ids.size());
	}


	TEST(FarCellsTest, DivideByThree) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.far_factory.cell_ids.size(), 32ull);

		EXPECT_EQ(adjfactory.far_factory.cell_intervals.back(), adjfactory.far_factory.cell_ids.size());

		EXPECT_EQ(adjfactory.far_factory.cell_ids[0], 2ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[1], 5ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[2], 6ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[3], 7ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[4], 8ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[5], 6ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[6], 7ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[7], 8ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[8], 0ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[9], 3ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[10], 6ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[11], 7ull);

		EXPECT_EQ(adjfactory.far_factory.cell_ids[12], 8ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[13], 2ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[14], 5ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[15], 8ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[16], 0ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[17], 3ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[18], 6ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[19], 0ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[20], 1ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[21], 2ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[22], 5ull);
		EXPECT_EQ(adjfactory.far_factory.cell_ids[23], 8ull);
	}

}


namespace TranslateOps
{
	size_t fact(size_t n)
	{
		if (n == 0)
			return 1;
		size_t res = 1;
		for (size_t i = 2; i <= n; i++)
			res = res * i;
		return res;
	}

	double nCr(size_t n, size_t r)
	{
		return (double)fact(n) / (fact(r) * fact(n - r));
	}


	MatrixXcd FindTifo(point_t cells_dif,size_t P, point_t cells_diff_neg)
	{
		MatrixXcd result(P, P);

		result(0, 0) = std::log(cells_diff_neg);

		for (size_t j = 1; j < P; ++j)
		{
			result(0, j) = std::pow(-1.0, j) / std::pow(cells_dif, j);
		}

		for (size_t i = 1; i < P; ++i)
		{
			result(i, 0) = -1.0 / ((point_t)i * std::pow(cells_dif, i));

			for (size_t j = 1; j < P; j++)
			{
				result(i, j) = std::pow(-1.0, j) * nCr(i+j-1,j-1) / std::pow(cells_dif, (i + j));
			}
		}
		return result;
	}



	TEST(Tofs, SingleCharge)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 1ull, domain };
		Factory factory{ adjfactory, 1ull };

		SortedData data{ factory.get_sources() };

		unsigned char P = 2;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		EXPECT_EQ(tras_op.T_ofs(0)(0, 0), point_t(1.0));
		EXPECT_EQ(tras_op.T_ofs(0)(1, 0), -point_t(1.0) * (data.point[0] - data.cell_center(0)));
	}

	TEST(Tofs, FourCellsTwoCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 2ull, domain };
		Factory factory{ adjfactory, 2ull };

		SortedData data{ factory.get_sources() };

		unsigned char P = 5;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		for (size_t cell_id = 0; cell_id < data.its_cell_center.size(); ++cell_id)
		{
			for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
			{
				size_t start_id{ data.interval_ids[cell_id] };
				// check p == 0
				EXPECT_EQ(tras_op.T_ofs(cell_id)(0ull, source), 1.0);
				for (size_t p = 1; p < P; ++p)
					// per sources in a cell
				{
					auto result = tras_op.T_ofs(cell_id)(p, source);
					auto expected = -1.0 / p * std::pow(data.point[start_id + source] - data.cell_center(cell_id), p);
					EXPECT_NEAR(result.real(), expected.real(), EPS);
					EXPECT_NEAR(result.imag(), expected.imag(), EPS);
				}
			}
		}
	}


	TEST(Tofs, NineCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };

		SortedData data{ factory.get_sources() };

		unsigned char P = 5;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		for (size_t cell_id = 0; cell_id < data.its_cell_center.size(); ++cell_id)
		{
			for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
			{
				size_t start_id{ data.interval_ids[cell_id] };
				EXPECT_EQ(tras_op.T_ofs(cell_id)(0ull, source), 1.0);
				for (size_t p = 1; p < P; ++p)
				{
					auto result = tras_op.T_ofs(cell_id)(p, source);
					auto expected = -1.0 / p * std::pow(data.point[start_id + source] - data.cell_center(cell_id), p);
					EXPECT_NEAR(result.real(), expected.real(), EPS);
					EXPECT_NEAR(result.imag(), expected.imag(), EPS);
				}
			}
		}
	}


	TEST(TofsMultiply, SingleCharge)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 1ull, domain };
		Factory factory{ adjfactory, 1ull };

		SortedData data{ factory.get_sources() };

		unsigned char P = 2;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		Eigen::VectorXcd result(tras_op.Outgoing_expansion());

		EXPECT_EQ(result(0), point_t(1.0));
		EXPECT_EQ(result(1), -point_t(1.0) * (data.point[0] - data.cell_center(0)) * 1.0);
	}


	TEST(TofsMultiply, NineCellsThreeChargesFiveP)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };

		SortedData data{ factory.get_sources() };

		unsigned char P = 5;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		Eigen::VectorXcd result(tras_op.Outgoing_expansion());

		Eigen::VectorXcd expected(P * data.its_cell_center.size());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());

		for (size_t cell_id = 0; cell_id < data.its_cell_center.size(); ++cell_id)
		{
			size_t start_id{ data.interval_ids[cell_id] };
			for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
			{
				expected(cell_id * P) += data.q[start_id + source];
			}
			for (size_t p = 1; p < P; ++p)
			{
				for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
				{
					expected(cell_id * P + p) += -1.0 / p *
						std::pow(data.point[start_id + source] - data.cell_center(cell_id), p)
							* data.q[start_id + source];
				}
			}
		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), EPS);
			EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}

	}


	TEST(TofsMultiply, FiveChargesInOneCellTwoP)
	{
		const std::vector<double>& x{ 4,1,-1,-1,-2 };
		const std::vector<double>& y{ 2,1,1,-2,2 };
		const std::vector<double>& q{ 1,3,2,-2,-1 };
		SortedData data{ x,y,q,std::vector<size_t>{ 0,5},std::vector <point_t>{} };

		unsigned char P = 2;

		Calculate_FMM::TranslateOperator tras_op{ data, P };

		Eigen::VectorXcd result(tras_op.Outgoing_expansion());

		EXPECT_EQ(result(0), point_t(3.0));
		EXPECT_EQ(result(1), -point_t(6.0, 9.0));
	}


	TEST(Tifo, NineCellsThreeChargesMultiply)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		auto result = t_ifo.Incoming_expansion();

		Eigen::VectorXcd expected = VectorXcd::Zero(P * factory.grid.cell_centers.size());

		EXPECT_EQ(result.size(), expected.size());

		Calculate_FMM::TranslateOperator tofs{ factory.get_sources(),P };

		auto sources = tofs.Outgoing_expansion();

		for (size_t cell_id = 0, far_cell_id, count_far, start_id;
			cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			VectorXcd sum_vector = VectorXcd::Zero(P);
			count_far = factory.grid.far_factory.cell_count[cell_id];
			start_id = factory.grid.far_factory.cell_intervals[cell_id];
			for (size_t source_id = 0; source_id < factory.grid.far_factory.cell_count[cell_id]; ++source_id)
			{
				far_cell_id = factory.grid.far_factory.cell_ids[start_id + source_id];
				point_t cells_diff = factory.grid.cell_centers[far_cell_id] - factory.grid.cell_centers[cell_id];
				point_t cells_diff_neg = factory.grid.cell_centers[cell_id] - factory.grid.cell_centers[far_cell_id];
				auto temp = FindTifo(cells_diff, P, cells_diff_neg);
				sum_vector += temp * sources.segment(P * far_cell_id, P);
			}
			expected.segment(cell_id * P, P) = sum_vector;

		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), EPS);
			EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}


	TEST(Tifo, SixteenCellsTwentyFiveChargesMultiply)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 4ull, domain };
		Factory factory{ adjfactory, 25ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		auto result = t_ifo.Incoming_expansion();

		Eigen::VectorXcd expected = VectorXcd::Zero(P * factory.grid.cell_centers.size());

		EXPECT_EQ(result.size(), expected.size());

		Calculate_FMM::TranslateOperator tofs{ factory.get_sources(),P };

		auto sources = tofs.Outgoing_expansion();

		for (size_t cell_id = 0, far_cell_id, count_far, start_id;
			cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			VectorXcd sum_vector = VectorXcd::Zero(P);
			count_far = factory.grid.far_factory.cell_count[cell_id];
			start_id = factory.grid.far_factory.cell_intervals[cell_id];
			for (size_t source_id = 0; source_id < factory.grid.far_factory.cell_count[cell_id]; ++source_id)
			{
				far_cell_id = factory.grid.far_factory.cell_ids[start_id + source_id];
				point_t cells_diff = factory.grid.cell_centers[far_cell_id] - factory.grid.cell_centers[cell_id];
				point_t cells_diff_neg = factory.grid.cell_centers[cell_id] - factory.grid.cell_centers[far_cell_id];
				auto temp = FindTifo(cells_diff, P, cells_diff_neg);
				sum_vector += temp * sources.segment(P * far_cell_id, P);
			}
			expected.segment(cell_id * P, P) = sum_vector;

		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), EPS);
			EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}


	TEST(Tifo, NineCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		Calculate_FMM::TranslateOperator tofs{ factory.get_sources(),P };

		auto sources = tofs.Outgoing_expansion();

		for (size_t cell_id = 0, far_cell_id, count_far, start_id;
			cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			auto tifo_matrix = t_ifo.T_ifo(cell_id);
			count_far = factory.grid.far_factory.cell_count[cell_id];
			start_id = factory.grid.far_factory.cell_intervals[cell_id];
			for (size_t source_id = 0; source_id < count_far; ++source_id)
			{
				far_cell_id = factory.grid.far_factory.cell_ids[start_id + source_id];
				point_t cells_diff = factory.grid.cell_centers[far_cell_id] - factory.grid.cell_centers[cell_id];
				auto result = tifo_matrix.block(0, P * source_id, P, P);
				point_t cells_diff_neg = factory.grid.cell_centers[cell_id] - factory.grid.cell_centers[far_cell_id];
				auto expected = FindTifo(cells_diff, P, cells_diff_neg);
				for (size_t i = 0; i < P; ++i)
				{
					for (size_t j = 0; j < P; ++j)
					{
						EXPECT_NEAR(result(i,j).real(), expected(i,j).real(), EPS);
						EXPECT_NEAR(result(i,j).imag(), expected(i,j).imag(), EPS);
					}
				}

			}

		}
	}
	TEST(Ttfi, NineCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		auto incoming_vector = t_ifo.Incoming_expansion();

		Calculate_FMM::Target_translate_operator t_tfi{ factory.get_sources(),P,incoming_vector };

		SortedData data{ factory.get_sources() };


		for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{

			size_t start_id{ data.interval_ids[cell_id] };
			for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
			{
				EXPECT_EQ(t_tfi.T_tfi(cell_id)(source, 0), 1.0);
				for (size_t p = 1; p < P; ++p)
				{
					auto result = t_tfi.T_tfi(cell_id)(source, p);
					auto expected = std::pow(data.point[start_id + source] - data.cell_center(cell_id), p);
					EXPECT_NEAR(result.real(), expected.real(), EPS);
					EXPECT_NEAR(result.imag(), expected.imag(), EPS);
				}
			}
		}
	}

	TEST(TtfiMultiply, NineCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		auto incoming_vector = t_ifo.Incoming_expansion();

		Calculate_FMM::Target_translate_operator t_tfi{ factory.get_sources(),P,incoming_vector };

		auto result = t_tfi.Final_potential();

		SortedData data{ factory.get_sources() };

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());

		for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			size_t start_id{ data.interval_ids[cell_id] };
			
			for (size_t source = 0; source < data.interval_count[cell_id]; ++source)
			{
				auto delta = data.point[start_id + source] - data.cell_center(cell_id);
				for (size_t p = 0; p < P; ++p)
				{
					expected(start_id + source) += std::pow(delta, p) * incoming_vector(cell_id * P + p);
				}
			}
			
		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), EPS);
			EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}


	TEST(FullResult, FourCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 2, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 3;

		Calculate_FMM::Solver fmm_solver{ factory, P };

		SortedData data{ factory.get_sources() };

		auto result = fmm_solver.SumAllPotential();

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());

		for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			//size_t start_id{ data.interval_ids[cell_id] };

			for (size_t source_id = 0; source_id < factory.grid.cell_centers.size(); ++source_id)
			{
				auto matrix = fmm_solver.CreateMatrix(cell_id, source_id);

				Map<VectorXd> q_source(data.q.data() + data.interval_ids[source_id],
					data.interval_count[source_id]);

				expected.segment(data.interval_ids[cell_id], data.interval_count[cell_id])
					+= matrix * q_source;

			}

		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), EPS);
			EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}

	TEST(FullResult, NineCellsOneCharge)
	{
		double eps{ 1e-1 };
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3, domain };
		Factory factory{ adjfactory, 1ull };
		unsigned char P = 3;

		Calculate_FMM::Solver fmm_solver{ factory, P };

		SortedData data{ factory.get_sources() };

		auto result = fmm_solver.SumAllPotential();

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());

		for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{

			for (size_t source_id = 0; source_id < factory.grid.cell_centers.size(); ++source_id)
			{
				auto matrix = fmm_solver.CreateMatrix(cell_id, source_id);

				Map<VectorXd> q_source(data.q.data() + data.interval_ids[source_id], 
					data.interval_count[source_id]);

				expected.segment(data.interval_ids[cell_id], data.interval_count[cell_id]) 
					+= matrix * q_source;

			}

		}

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), eps);
			//EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}

	TEST(FullResult, NineCellsOneChargeOnlyVectors)
	{
		double eps{ 1e-1 };
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3, domain };
		Factory factory{ adjfactory, 1ull };
		unsigned char P = 3;

		Calculate_FMM::Solver fmm_solver{ factory, P };

		SortedData data{ factory.get_sources() };

		auto result = fmm_solver.SumAllPotential();

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());


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

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), eps);
			//EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}

	TEST(FullResult, twentyOneCellsTwentytwoChargeOnlyVectors)
	{
		double eps{ 1e-1 };
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 21, domain };
		Factory factory{ adjfactory, 22ull };
		unsigned char P = 5;

		Calculate_FMM::Solver fmm_solver{ factory, P };

		SortedData data{ factory.get_sources() };

		auto result = fmm_solver.SumAllPotential();

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());


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

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), eps);
			//EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}

	TEST(FullResult, FarPotential)
	{
		double eps{ 1e-10 };
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 10ull };
		unsigned char P = 100;

		Calculate_FMM::Solver fmm_solver{ factory, P };

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		VectorXcd incoming_vector = t_ifo.Incoming_expansion();

		Calculate_FMM::Target_translate_operator t_tfi{ factory.get_sources(),P,incoming_vector };

		VectorXcd result{ t_tfi.Final_potential() };

		SortedData data{ factory.get_sources() };

		Eigen::VectorXcd expected(data.interval_ids.back());

		expected.setZero();

		EXPECT_EQ(result.size(), expected.size());


		for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			size_t start_id{ factory.grid.far_factory.cell_intervals[cell_id]};

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

		for (size_t i = 0; i < result.size(); i++)
		{
			EXPECT_NEAR(result(i).real(), expected(i).real(), eps);
			//EXPECT_NEAR(result(i).imag(), expected(i).imag(), EPS);
		}
	}


}