#pragma once
#include <vector>
#include "Data.h"


struct Factory
{
	double x_a, x_b, y_a, y_b;
	size_t m, q_per_block;

	std::vector<double> x_grid, y_grid;

	Factory(
		size_t m, 
		size_t q_per_block = 1ull, 
		double x_a = 0.0, double x_b = 1.0, double y_a = 0.0, double y_b = 1.0) :
		x_a{ x_a }, x_b{ x_b }, y_a{ y_a }, y_b{ y_b }, m{ m },
		q_per_block{ q_per_block },
		x_grid(m+1, 0.0),
		y_grid(m+1, 0.0)
	{
		double h_x{ (x_b - x_a) / m }, h_y{ (y_b - y_a) / m };
		for (size_t i = 0; i < m + 1; ++i)
		{
			x_grid[i] = x_a + i * h_x;
			y_grid[i] = y_a + i * h_y;
		}

		interval_ids.push_back(0ull);
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < m; ++j)
			{
				for (size_t q_id = 0; q_id < q_per_block; ++q_id)
				{
				//	auto [x0, y0, q0] = generate_source(i, j);

					//x.push_back(x0 /*x_grid[i] + h_x / 2.0*/);
					//y.push_back(y0 /*y_grid[j] + h_y / 2.0*/);
					//q.push_back(q0 /*1.0*/);

					x.push_back(x_grid[i] + h_x / 2.0);
					y.push_back(y_grid[j] + h_y / 2.0);
					q.push_back(1.0);
				}
				interval_ids.push_back(interval_ids.back() + q_per_block);
			}
		}
	}

	auto get_sources() const
	{
		return SortedData{ x,y,q,interval_ids };
	}

	std::vector<double> x, y, q;
	std::vector<size_t> interval_ids;

private:
	/*auto generate_source(size_t i, size_t j)
	{
		return Data{ random(x), random(y), random(q) };
	}*/

};
