#pragma once
#include <vector>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <random>


#include "Data.h"

using namespace DataStructs;

/**
 * @brief Descriptor of entire rectangular domain with point sources.
 */
struct Domain
{
	double x_a, x_b, y_a, y_b;
	double width, height;
	Domain(double x_a = 0.0, double x_b = 1.0, double y_a = 0.0, double y_b = 1.0) :
		x_a{ x_a }, x_b{ x_b }, y_a{ y_a }, y_b{ y_b },
		width{ x_b - x_a }, height{ y_b - y_a }
	{}
};

// DRY principle
struct SubDomain :
	public Domain
{
	using Domain::Domain;
};

struct GridFactory
{
	Domain domain;
	size_t m;
	// m x m grid within the domain
	std::vector<double> x_grid, y_grid;

	explicit GridFactory(
		size_t m,
		const Domain& domain = Domain{}) :
		domain{ domain }, m{ m },
		x_grid(m + 1, 0.0),
		y_grid(m + 1, 0.0)
	{
		double hx{ domain.width / m }, hy{ domain.height / m };
		for (size_t i = 0; i < m + 1; ++i)
		{
			x_grid[i] = domain.x_a + i * hx;
			y_grid[i] = domain.y_a + i * hy;
		}
	}

	SubDomain sub_domain(size_t i, size_t j) const
	{
		assert(i < m);
		assert(j < m);

		return SubDomain{ x_grid[i] , x_grid[i + 1] , y_grid[j] , y_grid[j + 1] };
	}

};


struct AdjacencyFactory
{
	std::vector<size_t> cell_ids;
	std::vector<size_t> cell_intervals;
	size_t m;

	explicit AdjacencyFactory(size_t m = 1) : m{ m }
	{
		set_adjacent_cells();
	}
protected:
	void set_adjacent_cells()
	{
		cell_intervals.reserve(m * m + 1);
		cell_ids.reserve(8 * m * m - 4 * (m + 1));

		cell_intervals.push_back(0ull);
		for (size_t j = 0, l = 0; j < m; ++j)
		{
			for (size_t i = 0; i < m; ++i, ++l)
			{
				size_t count{ 0ull };
				if (i > 0)
				{
					++count;
					cell_ids.push_back(l - 1ul);
				}
				if (i < m - 1ull)
				{
					++count;
					cell_ids.push_back(l + 1ul);
				}
				if (j > 0)
				{
					++count;
					cell_ids.push_back(l - m);
				}
				if (j < m - 1ull)
				{
					cell_ids.push_back(l + m);
					++count;
				}

				if (i > 0 && j > 0)
				{
					++count;
					cell_ids.push_back(l - 1ul - m);
				}
				if (i > 0 && (j < m - 1ull))
				{
					++count;
					cell_ids.push_back(l - 1ul + m);
				}
				if ((i < m - 1ull) && j > 0)
				{
					++count;
					cell_ids.push_back(l + 1ul - m);
				}
				if ((i < m - 1ull) && (j < m - 1ull))
				{
					++count;
					cell_ids.push_back(l + 1ul + m);
				}

				assert(count < 9);
				assert((m < 2) || (count > 2));
				std::sort(cell_ids.begin() + cell_intervals.back(), cell_ids.end());
				cell_intervals.push_back(cell_intervals.back() + count);
			}
		}
	}

};

struct FarFactory
{
	std::vector<size_t> cell_ids;
	std::vector<size_t> cell_intervals;
	size_t m;
	explicit FarFactory(size_t m) : m{ m }
	{
		set_far_cells();
	}

protected:
	void set_far_cells()
	{
		cell_intervals.reserve(m * m + 1);
		size_t area = m * m;
		cell_ids.reserve((area - 9) * (m - 2) * (m - 2)
			+ 4 * (area - 4) + (area - 6) * (m - 2) * 4);

		cell_intervals.push_back(0ull);

		size_t count{ 0ull };

		for (size_t row = 0; row < m; ++row)
		{
			for (size_t col = 0; col < m; ++col)
			{
				count = 0ull;
				//Все дальние клетки справа на этой линии
				if (m > 2 && col < m - 2)
				{
					for (size_t line = col + 2; line < m; ++line)
					{
						cell_ids.push_back(m * row + line);
					}
					count += (m - col - 2);
				}

				//Все дальние клетки справа на линии ниже
				if (m > 2 && row < m - 1 && col < m - 2)
				{
					for (size_t line = col + 2; line < m; ++line)
					{
						cell_ids.push_back(m * (row + 1) + line);
					}
					count += (m - col - 2);
				}

				//Все дальние клетки справа на линии выше
				if (m > 2 && row > 0 && col < m - 2)
				{
					for (size_t line = col + 2; line < m; ++line)
					{
						cell_ids.push_back(m * (row - 1) + line);
					}
					count += (m - col - 2);
				}


				//Все дальние клетки слева на этой линии 
				if (col > 1)
				{
					for (size_t line = 0; line < col - 1; ++line)
					{
						cell_ids.push_back(m * row + line);
					}
					count += (col - 1);
				}

				//Все дальние клетки слева на нижней линии
				if (row < m - 1 && col > 1)
				{
					for (size_t line = 0; line < col - 1; ++line)
					{
						cell_ids.push_back(m * (row + 1) + line);
					}
					count += (col - 1);
				}

				//Все дальние клетки слева на верхней линии
				if (row > 0 && col > 1)
				{
					for (size_t line = 0; line < col - 1; ++line)
					{
						cell_ids.push_back(m * (row - 1) + line);
					}
					count += (col - 1);
				}

				//Все клетки ниже
				if (row > 1)
				{
					for (size_t local_row = 0; local_row < row - 1; ++local_row)
					{
						for (size_t i = 0; i < m; i++)
						{
							cell_ids.push_back(local_row * m + i);
						}
					}
					count += (m * (row - 1));
				}

				//Все клетки выше
				if (m > 2 && row < m - 2)
				{
					for (size_t local_row = row + 2; local_row < m; ++local_row)
					{
						for (size_t i = 0; i < m; i++)
						{
							cell_ids.push_back(local_row * m + i);
						}
					}
					count += (m * (m - row - 2));
				}

				std::sort(cell_ids.begin() + cell_intervals.back(), cell_ids.end());
				cell_intervals.push_back(cell_intervals.back() + count);

			}
		}
	}
};


struct BigAdjacencyFactory:
	public GridFactory
{
	using point_t = Eigen::dcomplex;
	std::vector<point_t> cell_centers;
	struct AdjacencyFactory adjacency_factory;
	struct FarFactory far_factory;

	explicit BigAdjacencyFactory(
		size_t m,
		const Domain& domain = Domain{}) :
		GridFactory(m, domain),
		adjacency_factory{m},
		far_factory{m}
	{
		double hx{ domain.width / m }, hy{ domain.height / m };
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < m; ++j)
				cell_centers.emplace_back(x_grid[j] + hx / 2, y_grid[i] + hy / 2);
	}

};



/**
 * @brief Factory that introduces a grid on the domain
 * and generates the adjacency list of adjoining cells.
 */

//adjFactory

//struct AdjacencyFactory
//{
//	using point_t = Eigen::dcomplex;
//	Domain domain;
//	size_t m;
//	// m x m grid within the domain
//	std::vector<double> x_grid, y_grid;
//	// descriptor of adjacent cells
//	// cell_ids --- contiguous cell ids adjacent to a cell I
//	// [cell_intervals[I]; cell_intervals[I+1]] --- contains the I-adjacent cells
//	std::vector<size_t> cell_ids, cell_intervals;
//	std::vector<point_t> cell_centers;
//
//	explicit AdjacencyFactory(
//		size_t m,
//		const Domain& domain = Domain{}) :
//		domain{ domain }, m{ m },
//		x_grid(m + 1, 0.0),
//		y_grid(m + 1, 0.0)
//	{
//		double hx{ domain.width / m }, hy{ domain.height / m };
//		for (size_t i = 0; i < m + 1; ++i)
//		{
//			x_grid[i] = domain.x_a + i * hx;
//			y_grid[i] = domain.y_a + i * hy;
//		}
//		for (size_t i = 0; i < m; ++i)
//			for (size_t j = 0; j < m; ++j)
//				cell_centers.emplace_back(x_grid[j] + hx / 2, y_grid[i] + hy / 2);
//		set_adjacent_cells();
//	}
//
//	SubDomain sub_domain(size_t i, size_t j) const
//	{
//		assert(i < m);
//		assert(j < m);
//
//		return SubDomain{ x_grid[i] , x_grid[i+1] , y_grid[j] , y_grid[j+1] };
//	}
//
//protected:
//	void set_adjacent_cells()
//	{
//		cell_intervals.reserve(m * m + 1);
//		cell_ids.reserve(8 * m * m - 4 * (m + 1));
//
//		cell_intervals.push_back(0ull);
//		for (size_t j = 0, l = 0; j < m; ++j)
//		{
//			for (size_t i = 0; i < m; ++i, ++l)
//			{
//				size_t count{ 0ull };
//				if (i > 0)
//				{
//					++count;
//					cell_ids.push_back(l - 1ul);
//				}
//				if (i < m - 1ull)
//				{
//					++count;
//					cell_ids.push_back(l + 1ul);
//				}
//				if (j > 0)
//				{
//					++count;
//					cell_ids.push_back(l - m);
//				}
//				if (j < m - 1ull)
//				{
//					cell_ids.push_back(l + m);
//					++count;
//				}
//
//				if (i > 0 && j > 0)
//				{
//					++count;
//					cell_ids.push_back(l - 1ul - m);
//				}
//				if (i > 0 && (j < m - 1ull))
//				{
//					++count;
//					cell_ids.push_back(l - 1ul + m);
//				}
//				if ((i < m - 1ull) && j > 0)
//				{
//					++count;
//					cell_ids.push_back(l + 1ul - m);
//				}
//				if ((i < m - 1ull) && (j < m - 1ull))
//				{
//					++count;
//					cell_ids.push_back(l + 1ul + m);
//				}
//
//				assert(count < 9);
//				assert((m < 2) || (count > 2));
//				std::sort(cell_ids.begin() + cell_intervals.back(), cell_ids.end());
//				cell_intervals.push_back(cell_intervals.back() + count);
//			}
//		}
//	}
//};

// https://en.cppreference.com/w/cpp/numeric/random
// https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c
// https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
/**
 * @brief Generator of point charges within every subcell 
 * of the main domain
 */
struct Factory
	//:public GridFactory
{
	size_t q_per_block;
	BigAdjacencyFactory grid;
	//AdjCellsFactory
	//FarCellsFactory

	explicit Factory(
		const BigAdjacencyFactory& grid,
		size_t q_per_block = 1ull) :
		q_per_block{ q_per_block },
		grid{grid},
		rd{}
	{
		size_t count{ grid.m * grid.m * q_per_block };
		x.reserve(count);
		y.reserve(count);
		q.reserve(count);
		interval_ids.reserve(grid.m* grid.m + 1ull);

		std::mt19937 gen(rd());
		interval_ids.push_back(0ull);
		for (size_t j = 0, l =0; j < grid.m; ++j)
		{
			for (size_t i = 0; i < grid.m; ++i, ++l)
			{
				// for every subdomain of the main domain
				// generate a set of point charges
				auto domain{ grid.sub_domain(i,j) };
				std::uniform_real_distribution<> uniform_x(domain.x_a, domain.x_b);
				std::uniform_real_distribution<> uniform_y(domain.y_a, domain.y_b);
				for (size_t q_id = 0; q_id < q_per_block; ++q_id)
				{
					x.push_back(uniform_x(gen));
					y.push_back(uniform_y(gen));
					q.push_back(1.0);
				}
				interval_ids.push_back(interval_ids.back() + q_per_block);
			}
		}
	}


	virtual void set_data()
	{
		size_t count{ grid.m * grid.m * q_per_block };
		x.reserve(count);
		y.reserve(count);
		q.reserve(count);
		interval_ids.reserve(grid.m * grid.m + 1ull);

		std::mt19937 gen(rd());
		interval_ids.push_back(0ull);
		for (size_t j = 0, l = 0; j < grid.m; ++j)
		{
			for (size_t i = 0; i < grid.m; ++i, ++l)
			{
				// for every subdomain of the main domain
				// generate a set of point charges
				auto domain{ grid.sub_domain(i,j) };
				std::uniform_real_distribution<> uniform_x(domain.x_a, domain.x_b);
				std::uniform_real_distribution<> uniform_y(domain.y_a, domain.y_b);
				for (size_t q_id = 0; q_id < q_per_block; ++q_id)
				{
					x.push_back(uniform_x(gen));
					y.push_back(uniform_y(gen));
					q.push_back(1.0);
				}
				interval_ids.push_back(interval_ids.back() + q_per_block);
			}
		}
	}

	auto get_sources() const
	{
		return SortedData{ x,y,q,interval_ids, grid.cell_centers };
	}

	std::vector<double> x, y, q;
	std::vector<size_t> interval_ids;

private:
	// Seed with a real random value
	std::random_device rd;
};
