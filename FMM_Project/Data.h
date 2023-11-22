#pragma once
#include <vector>
#include <complex>
#include <cassert>

struct Data
{
	Data(double x, double y, double q):
		x{ x }, y{ y }, q{ q }
	{}

	Data(const Data&) = default;

	double x;
	double y;
	double q;
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
	{}

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

	}

	SortedData(
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& q,
		const std::vector<size_t>& interval_ids,
		const std::vector<point_t>& cell_center) :
		SortedData{UseDataType{ x,y,q },  interval_ids, cell_center }
	{	}

	auto cell_center(size_t cell_id) const
	{
		assert(cell_id < its_cell_center.size());
		return its_cell_center[cell_id];
	}

	auto cell_sources(size_t cell_id) const
	{
		assert(cell_id < its_cell_center.size());
		return Eigen::VectorXd::Map(q.data() + interval_ids[cell_id], interval_ids[cell_id + 1ull] - interval_ids[cell_id]);
	}
};




