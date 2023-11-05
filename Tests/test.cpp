#include "pch.h"
#include "../FMM_Project/Factory.h"

// https://learn.microsoft.com/ru-ru/visualstudio/test/how-to-use-google-test-for-cpp?view=vs-2022
// https://github.com/google/googletest/blob/main/docs/primer.md
// https://learn.microsoft.com/ru-ru/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022

namespace DataGenerators
{
	TEST(DomainTest, DefaualtState) {

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
		AdjacencyFactory adjfactory{ 1ull, domain };

		EXPECT_EQ(adjfactory.x_grid[0], -1.0);
		EXPECT_EQ(adjfactory.x_grid[1], 1.0);
		EXPECT_EQ(adjfactory.y_grid[0], -1.0);
		EXPECT_EQ(adjfactory.y_grid[1], 1.0);
		EXPECT_EQ(adjfactory.cell_ids.size(), 0ull);
		EXPECT_EQ(adjfactory.cell_intervals.front(), 0ull);

		EXPECT_EQ(adjfactory.cell_intervals.back(), adjfactory.cell_ids.size());
	}

	TEST(AdjacencyTest, DivideByTwo) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		AdjacencyFactory adjfactory{ 2ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.cell_ids.size(), 12ull);

		EXPECT_EQ(adjfactory.cell_intervals.back(), adjfactory.cell_ids.size());

		for (size_t id = 0ull; id < adjfactory.cell_intervals.size(); ++id)
			EXPECT_EQ(adjfactory.cell_intervals[id], id * 3ull);

		EXPECT_EQ(adjfactory.cell_ids[0], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[1], 2ull);
		EXPECT_EQ(adjfactory.cell_ids[2], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[3], 0ull);
		EXPECT_EQ(adjfactory.cell_ids[4], 2ull);
		EXPECT_EQ(adjfactory.cell_ids[5], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[6], 0ull);
		EXPECT_EQ(adjfactory.cell_ids[7], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[8], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[9], 0ull);
		EXPECT_EQ(adjfactory.cell_ids[10], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[11], 2ull);
	}

	TEST(AdjacencyTest, DivideByThree) {
		Domain domain{ -1.0, 1.0, -1.0, 1.0 };
		AdjacencyFactory adjfactory{ 3ull, domain };

		EXPECT_EQ(adjfactory.x_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.x_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.y_grid.front(), -1.0);
		EXPECT_EQ(adjfactory.y_grid.back(), 1.0);
		EXPECT_EQ(adjfactory.cell_ids.size(), 40ull);

		EXPECT_EQ(adjfactory.cell_intervals.back(), adjfactory.cell_ids.size());


		EXPECT_EQ(adjfactory.cell_ids[0], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[1], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[2], 4ull);
		EXPECT_EQ(adjfactory.cell_ids[3], 0ull);
		EXPECT_EQ(adjfactory.cell_ids[4], 2ull);
		EXPECT_EQ(adjfactory.cell_ids[5], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[6], 4ull);
		EXPECT_EQ(adjfactory.cell_ids[7], 5ull);
		EXPECT_EQ(adjfactory.cell_ids[8], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[9], 4ull);
		EXPECT_EQ(adjfactory.cell_ids[10], 5ull);
		EXPECT_EQ(adjfactory.cell_ids[11], 0ull);

		EXPECT_EQ(adjfactory.cell_ids[12], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[13], 4ull);
		EXPECT_EQ(adjfactory.cell_ids[14], 6ull);
		EXPECT_EQ(adjfactory.cell_ids[15], 7ull);
		EXPECT_EQ(adjfactory.cell_ids[16], 0ull);
		EXPECT_EQ(adjfactory.cell_ids[17], 1ull);
		EXPECT_EQ(adjfactory.cell_ids[18], 2ull);
		EXPECT_EQ(adjfactory.cell_ids[19], 3ull);
		EXPECT_EQ(adjfactory.cell_ids[20], 5ull);
		EXPECT_EQ(adjfactory.cell_ids[21], 6ull);
		EXPECT_EQ(adjfactory.cell_ids[22], 7ull);
		EXPECT_EQ(adjfactory.cell_ids[23], 8ull);
	}
}