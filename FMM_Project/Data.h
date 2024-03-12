#pragma once
#include <vector>
#include <complex>
#include <cassert>

#include <Eigen/Core>

namespace DataStructs
{
	struct Data
	{
		Data(double x, double y, double q) :
			x{ x }, y{ y }, q{ q }
		{}

		Data(const Data&) = default;

		double x;
		double y;
		double q;

	};

	struct y_sort
	{
		bool operator() (const Data& lhs, const Data& rhs) const
		{
			return lhs.y < rhs.y;
		}
	};

	struct x_sort
	{
		bool operator() (const Data& lhs, const Data& rhs) const
		{
			return lhs.x < rhs.x;
		}
	};


	struct InputDataType
	{
		std::vector<Data> data;

		size_t size() const
		{
			return data.size();
		}
		const Data& operator[](size_t i) const
		{
			return data[i];
		}

		void Sort_InputDataType()
		{
			std::sort(data.begin(), data.end(), y_sort());

			for (size_t i = 0, start_idx = 0, idx = 0; idx < size(); i++)
			{
				if (data[idx].y < y_grid[i+1])
				{
					continue;
				}

				std::sort(data.begin() + start_idx, data.begin() + idx, x_sort());
				++i;
				start_idx = idx;
			}

			std::sort(data.begin() + start)

		}

	};

	struct UseDataType
	{
		using point_t = Eigen::dcomplex;
		std::vector<double> x, y, q;
		std::vector<point_t> point;

		UseDataType(
			const std::vector<double>& x,
			const std::vector<double>& y,
			const std::vector<double>& q) :
			x{ x }, y{ y }, q{ q }
		{
			//point.reserve(x.size());

			for (size_t i = 0; i < x.size(); ++i)
			{
				point.push_back(point_t{ x[i], y[i] });
			}
		}

		UseDataType(const InputDataType& data)
		{
			x.reserve(data.size());
			y.reserve(data.size());
			q.reserve(data.size());
			point.reserve(data.size());

			for (size_t i = 0; i < data.size(); ++i)
			{
				x.push_back(data[i].x);
				y.push_back(data[i].y);
				q.push_back(data[i].q);

				point[i] = point_t{ x[i], y[i] };
			}
		}
	};

	struct SortedData :
		public UseDataType
	{
		std::vector<size_t> interval_ids;
		std::vector<size_t> interval_count;	
		std::vector<point_t> its_cell_center;

		SortedData(
			const UseDataType& data,
			const std::vector<size_t>& interval_ids,
			const std::vector<point_t>& cell_center) :
			UseDataType{ data },
			interval_ids{ interval_ids },
			its_cell_center{ cell_center }
		{
			if (interval_ids[0] != 0ull)
				throw std::exception();

			interval_count.reserve(interval_ids.size() - 1);
			for (size_t i = 0; i < interval_ids.size() - 1; ++i)
				interval_count.push_back(interval_ids[i + 1] - interval_ids[i]);
			if (its_cell_center.size() == 0)
			{
				its_cell_center.reserve(interval_ids.size() - 1);
				for (size_t i = 0; i < interval_ids.size() - 1; ++i)
				{
					its_cell_center.push_back(Find_center(i));
				}
			}
		}


		void Sort()
		{
			std::sort(data)
		}


		SortedData(
			const std::vector<double>& x,
			const std::vector<double>& y,
			const std::vector<double>& q,
			const std::vector<size_t>& interval_ids,
			const std::vector<point_t>& cell_center) :
			SortedData{ UseDataType{ x,y,q },  interval_ids, cell_center }
		{	}

		auto cell_center(size_t cell_id) const
		{
			assert(cell_id < its_cell_center.size());
			return its_cell_center[cell_id];
		}

		auto cell_sources(size_t cell_id) const
		{
			assert(cell_id < its_cell_center.size());
			return Eigen::VectorXd::Map(q.data() + interval_ids[cell_id],
				interval_ids[cell_id + 1ull] - interval_ids[cell_id]);
		}

		point_t Find_center(size_t cell_id)
		{
			size_t start_id{ interval_ids[cell_id] };
			size_t end_id{ interval_ids[cell_id + 1ull] };
			double c_x = ((*std::max_element(x.begin() + start_id, x.begin() + end_id))
				+ *std::min_element(x.begin() + start_id, x.begin() + end_id)) * 0.5;
			double c_y = ((*std::max_element(y.begin() + start_id, y.begin() + end_id))
				+ *std::min_element(y.begin() + start_id, y.begin() + end_id)) * 0.5;
			return point_t{ c_x, c_y };

			//	return data.cell_center(start_id);
		}
	};
}




