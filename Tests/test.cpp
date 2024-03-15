#include "pch.h"
#include "../FMM_Project/Factory.h"
#include "../FMM_Project/Transformations.h"
#define EPS 1e-15

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


	TEST(Tifo, NineCellsThreeCharges)
	{
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		BigAdjacencyFactory adjfactory{ 3ull, domain };
		Factory factory{ adjfactory, 3ull };
		unsigned char P = 5;

		Calculate_FMM::Incoming_translate_operator t_ifo{ factory, P };

		auto result = t_ifo.Incoming_expansion();

		std::cout << result;


		/*for (size_t cell_id = 0; cell_id < factory.grid.cell_centers.size(); ++cell_id)
		{
			MatrixXcd t_ifo_matrix = t_ifo.T_ifo(cell_id);



		}*/

		/*EXPECT_EQ(t_ifo.T_ifo(0)(0, 0), point_t(1.0));
		EXPECT_EQ(t_ifo.T_ifo(0)(1, 0), -point_t(1.0) * (data.point[0] - data.cell_center(0)));*/
	}

}