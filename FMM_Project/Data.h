#pragma once
#include <vector>


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
	std::vector<double> x, y, q;

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

		for (size_t i = 0; i < data.size(); ++i)
		{
			x.push_back(data[i].x);
			y.push_back(data[i].y);
			q.push_back(data[i].q);
		}
	}
};

struct SortedData :
	public UseDataType
{
	std::vector<size_t> interval_ids;

	SortedData(
		const UseDataType& data,
		const std::vector<size_t>& interval_ids) :
		UseDataType{ data },
		interval_ids{ interval_ids }
	{	
		if (interval_ids[0] != 0ull)
			throw std::exception();	
	}

	SortedData(
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& q,
		const std::vector<size_t>& interval_ids) :
		SortedData{UseDataType{ x,y,q },  interval_ids }
	{	}
};




