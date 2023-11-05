#include "pch.h"
#include "../FMM_Project/Factory.h"

// https://learn.microsoft.com/ru-ru/visualstudio/test/how-to-use-google-test-for-cpp?view=vs-2022
// https://github.com/google/googletest/blob/main/docs/primer.md
// https://learn.microsoft.com/ru-ru/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022

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
	Domain domain{-1.0, 1.0, -1.0, 1.0};
	AdjacencyFactory adjfactory{ 1ull, domain };

	EXPECT_EQ(adjfactory.x_grid[0], -1.0);
	EXPECT_EQ(adjfactory.x_grid[1], 1.0);
	EXPECT_EQ(adjfactory.y_grid[0], -1.0);
	EXPECT_EQ(adjfactory.y_grid[1], 1.0);
	EXPECT_EQ(adjfactory.cell_ids.size(), 0ull);
	EXPECT_EQ(adjfactory.cell_intervals[0], 0ull);
	EXPECT_EQ(adjfactory.cell_intervals[1], 0ull);
}